% ===========================================================
% Filename:     generate.m
% Date:   	 	2025-08-14
% Author:    	Fukang Qiao
% Description:  Main function for simulation multiple numbers
% ===========================================================

function generate(params)
saveFolderNames = fieldnames(params.out);
for idx = length(saveFolderNames):-1:1
    if ~params.out.(saveFolderNames{idx}); saveFolderNames(idx)=[]; end
end
createSubfolder(params.savePath);
for idx=1:length(saveFolderNames)
    createSubfolder(params.savePath,saveFolderNames{idx});
    if params.savePNGFlag && ~ismember(saveFolderNames{idx},{'deformBbox'})
        createSubfolder(params.savePath,[saveFolderNames{idx} 'Show']);
    end
end
% DEM文件
fileNames = listdir(params.demFolder,'/*.tif');
fileNum = length(fileNames);

if ~fileNum
    cycleNum = 1;
    warning('DEM file not found, DEM not used!');
else
    cycleNum = fileNum;
end
singleNum = ceil(params.totalNum/cycleNum);%如果多个DEM，则根据DEM数量分组处理

tic
for fileId=1:cycleNum
    if fileNum
        fileName = fileNames{fileId};
        % 读取DEM高程数据
        dem = imread(fullfile(params.demFolder, fileName));
        dem = double(dem);

        % Incidence angle range	18.3° - 46.8°
        incidenceAngle = randR([18.3,46.8]); % 32.55; 入射角
        wavelength = 0.0560; % unit: meter 波长
        LosPhase = 4*pi*dem*cos(incidenceAngle)/wavelength;
        LosPhase(dem==-32768|dem==inf)=nan;
    else
        LosPhase=nan;
    end
    % Total Numbers
    ibegin = (fileId-1)*singleNum;
    iend = fileId*singleNum-1; if iend>params.totalNum-1; iend=params.totalNum-1; end
    % 并行parallel
    if params.Parallel
        parfor i = ibegin:iend
            randomClipping(params,fileNum,LosPhase,fileId,i,saveFolderNames);
        end
    else
        for i = ibegin:iend
            randomClipping(params,fileNum,LosPhase,fileId,i,saveFolderNames);
        end
    end
end
toc
disp(['average time: ' num2str(toc/params.totalNum)]);
disp('over!');

end

%% Sub-functions for parallel processing
function randomClipping(params,fileNum,LosPhase,fileId,i,saveFolderNames)
if fileNum
    [m,n] = size(LosPhase);
    r1 = randi(m-params.sampleSize+1);
    c1 = randi(n-params.sampleSize+1);
    img = LosPhase(r1:r1+params.sampleSize-1,c1:c1+params.sampleSize-1);

    maxV = max2(img); minV=min2(img);
    if maxV~=minV
        img = (img-min2(img))/(max2(img)-min2(img));
        img = img.*randR([-6*pi,6*pi]); %amp=6*pi;;
    end
else
    img = zeros(params.sampleSize,params.sampleSize);
end

outputs = generateOne(img,params);

name = num2str(i,'%05d');
for idx=1:length(saveFolderNames)
    folderName = saveFolderNames{idx};
    if ~isfield(outputs,folderName);error([folderName ' is not included in the generated data, please check!']);end
    if isnan(outputs.(folderName)); warning(['outputs.' folderName ' is nan, not saved.']);continue;end
    if ismember(folderName,{'deformBbox'})
        if~isempty(outputs.(folderName))
            writematrix(outputs.(folderName),fullfile(params.savePath,folderName,[name,'.txt']),"Delimiter"," ");
        else
            %输出为空矩阵txt文件，便于YOLO样本训练
            writematrix([ ],fullfile(params.savePath,folderName,[name,'.txt']),"Delimiter"," ");
        end
    else
        % 修改输出文件后缀格式
        imwritebin3(outputs.(folderName),fullfile(params.savePath,folderName,[name,'.dat']));
        if params.savePNGFlag; imwrite(matToRGB(outputs.(folderName)),fullfile(params.savePath,[folderName 'Show'],[name,'.png'])); end
    end
    
end

disp(['fileId:' num2str(fileId) '    ' 'i:' num2str(i+1) '/' num2str(params.totalNum)]);
end