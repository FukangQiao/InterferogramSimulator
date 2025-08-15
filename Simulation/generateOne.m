% ===========================================================
% Filename:     generateOne.m
% Date:   	 	2025-08-14
% Author:    	Fukang Qiao
% Description:  Main function for simulation one
% ===========================================================

function outputs = generateOne(terrainRelated,params) %[50,50]
addpath('WZPUtil');

[rows,cols] = size(terrainRelated);
deformBbox = [];

%% Phase components %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phaseComponent = zeros(size(terrainRelated));
%%%%%% Randomly shifted phase
if 1
    ampRange = [-pi,pi]; % Amplitude range
    comp = randR(ampRange);
    phaseComponent = phaseComponent+comp;
end
%%%%%% Slope phase
if rand<params.probSlop
    ampRange = [0,cols/4.*pi]; % Amplitude range
    comp = getSlop(rows,cols,ampRange);
    phaseComponent = phaseComponent+comp;
end
%%%%%% Building phase
if rand<params.probBuilding
    numRange = [0,20]; % Number of buildings (integer)
    ampRange = [pi/2,4*pi]; % Amplitude range
    comp = getBuilding(rows,cols,numRange,ampRange);
    phaseComponent = phaseComponent+comp;
end
%%%%%% Atmospheric turbulence phase, i.e., fractal Perlin noises
if rand<params.probTurbulence
    ampRange = [0,pi]; % Amplitude range
    zmat = fractalPerlinNoise(rows,cols);
    comp=randR(ampRange)*zmat;
    phaseComponent = phaseComponent+comp;
end
%% Deformation components %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
deformComponent = zeros(size(terrainRelated));% 形变相位
%%%%%% Distorted two-dimensional Gaussian surface
if rand<params.probDeform
    ampRange = [-20*pi,20*pi]; % Amplitude range +-20pai
    deformNorm = getOriginNorm(rows,cols);
    deformAmp = randR(ampRange); %生成强度范围内的随机数
    deformation = deformAmp.*deformNorm; % Phase range
    deformComponent = deformComponent+deformation;

    if params.out.deformBbox && abs(deformAmp)>(pi/2) % If the maximum deformation amplitude is too small, it is ignored.
        % 由于二维高斯模拟 导致的形变为 Label=1
        [bbox,columnLabel]=getBoundingBox(abs(deformation)>(pi/3),1,0);
        deformBbox = [deformBbox;bbox];
%         showBbox(wrapToPi(deformation),bbox);colormap jet;colorbar;
    end
end

%%%%%% Deformation caused by earthquakes
if rand<params.probEarthquake
    faltParms.rows = rows;
    faltParms.cols = cols;
    faltParms.length=randR([4, faltParms.cols/3]);
    faltParms.width=randR([0.8, 1.2])*faltParms.length;
    faltParms.depth=randR([1, 20]);
    faltParms.dip=randR([-140, -40]);  % -90:1:0
    faltParms.strike=randR([0, 360]);
    faltParms.xloc=randR([0.1,0.9])*faltParms.cols;
    faltParms.yloc=randR([0.1,0.9])*faltParms.rows;
    faltParms.sslip=randR([-1.5,1.5]);
    faltParms.dslip=randR([-1.5,1.5]);
    faltParms.opening=(randR>0.95)*randR([0,10]);
    if rand>0.5; faltParms.ATorDT='AT'; else ; faltParms.ATorDT='DT'; end % AT or DT
    deformation=getFalts(faltParms);%是否为disloc参数
    deformComponent = deformComponent+deformation;

    if params.out.deformBbox
        % 由于Earthquake 导致的形变为 Label=2
        [bbox,columnLabel]=getBoundingBox(abs(deformation)>(pi/2),2,0);
        deformBbox = [deformBbox;bbox];
%         showBbox(wrapToPi(deformation),bbox);colormap jet;colorbar;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Phase combination %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
origin = terrainRelated + phaseComponent+deformComponent;
origin(isnan(terrainRelated)|isinf(terrainRelated))=0;
originWrapped = wrapToPi(origin);
if params.out.origin; outputs.origin = origin; end
if params.out.deformBbox; outputs.deformBbox = deformBbox; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 增加复数随机噪声
%%% Complex-valued noise with random amplitude %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 噪声与梯度有关
if params.noiseType==0 
    [waterMask,interf,coherence] = addComplexNoise(...
        imresize(origin,[rows*params.multilook(1),cols*params.multilook(2)]), ...
        imresize(terrainRelated,[rows*params.multilook(1),cols*params.multilook(2)]),...
        imresize(deformComponent,[rows*params.multilook(1),cols*params.multilook(2)]),...
        params);
    originWrapped(waterMask)=interf(waterMask); % Reserved waters
    
elseif params.noiseType==1
    snrV = randR(params.noiseSNRRange);
    interf = awgn(complex(exp(1i*imresize(origin,[rows*params.multilook(1),cols*params.multilook(2)]))),snrV,'measured'); %加入信噪比为*db的噪声，加入前预估信号的功率（强度）。
    interf = angle(multilook(interf,params.multilook(1),params.multilook(2)));
    coherence = nan;
else
    error(['Unsupported noise type: ' num2str(params.noiseType)]);
end

% interf = goldstein_filt(interf,32,0.5);
if params.out.interf; outputs.interf = interf; end
if params.out.coherence; outputs.coherence = coherence; end
if params.out.originWrapped; outputs.originWrapped = originWrapped; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Generate data related to branch-cuts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if params.out.VDRI  % Horizontal phase gradient + residual map
    % + vertical phase gradient + interferogram (4 channels)
    residues = getResidues(interf);
    FxIntf = interf(:,2:end)-interf(:,1:end-1); FxIntf(:,end+1)=0;
    FxIntf = wrapToPi(FxIntf);
    FyIntf = interf(2:end,:)-interf(1:end-1,:); FyIntf(end+1,:)=0;
    FyIntf = wrapToPi(FyIntf);
    VDRI = cat(3,FxIntf,residues/2,FyIntf,interf); % Vertical+Residual+Horizontal
    
    outputs.VDRI = VDRI;
end
% label
if params.out.branchCut
    difference = angle(exp(1i*(interf-origin)));
    od = origin+difference; % True value, used to generate label
    
    FxLabel = od(:,2:end)-od(:,1:end-1); FxLabel(:,end+1)=0;
    FxLabel = abs(FxLabel)>pi;
    FyLabel = od(2:end,:)-od(1:end-1,:); FyLabel(end+1,:)=0;
    FyLabel = abs(FyLabel)>pi;
    branchCut = cat(3,FxLabel,FyLabel);
    
    outputs.branchCut = branchCut;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Customized Output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% You can customize other types of output and implement them here.
% Note that the variable names are consistent with the names in the
% parameter settings in "main.m".


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end


%% Add complex-valued deformation-related noise
function [waterMask,originNoise,coh] = addComplexNoise(origin,terrainRelated,deformation,params)
% 设置最小值为0
tamp = terrainRelated-min2(terrainRelated);
% 计算局部标准差和均值
% 在8*8窗口计算形变相位的标准差和均值
[dataStd, srcDataMean] = getLocalVariance(deformation, 8);
dsum = dataStd*1+srcDataMean*0.15;

noiseAmp = tamp+dsum;
if max(max(noiseAmp))~=0
    noiseAmp = noiseAmp/max(max(noiseAmp)); % 归一化[0,1]
    noiseAmp = noiseAmp.^randR([0.1,0.7]); % 幂函数随机分布
end

noiseFunnelScale = randR([0,0.5]);
noiseFunnelNorm = noiseAmp*(1-noiseFunnelScale)+noiseFunnelScale;

noiseAmp = randR([0,1.2])*pi.*(noiseFunnelNorm);

%%%%%% Completely decorrelated area (water area)
[m,n]=size(terrainRelated);
% 生成水域、低相干性区域
if rand<params.probWater
    maxx=randsrc(1,1,[2, 3, 4, 5, 6]);
    maxy=randsrc(1,1,[2, 3, 4, 5, 6]);
    roriginNorm = fractalPerlinNoise(n,m,1,1,maxx,maxy);
    thres = randR([0,1]);
    waterMaskGen = imbinarize(roriginNorm,thres);
else
    waterMaskGen = zeros(m,n);
end

waterMask = isnan(terrainRelated); % The area with DEM = -32768 is considered as water area
waterMask = waterMask|waterMaskGen;

noiseAmp(waterMask) = randR([5,10])*pi; % 设置水域为大幅度噪声
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

coh=calCohFromNoiseStd(noiseAmp,params.multilook,1);% 相干性

noiseR = noiseAmp .* randn(size(origin));
noiseI = noiseAmp .* randn(size(origin));

origin(waterMask)=0;
originNoise = cos(origin)+noiseR+1i*(sin(origin)+noiseI);
originNoise = multilook(originNoise,params.multilook(1),params.multilook(2));
originNoise = angle(originNoise);

waterMask = multilook(waterMask,params.multilook(1),params.multilook(2))>0.5;
end