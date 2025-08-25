% ===========================================================
% Filename:     getOriginNorm.m
% Date:   	 	2025-08-14
% Author:    	Fukang Qiao
% Description:  Function for simulate deformation of points
% ===========================================================

%% getOriginNorm
function result = getOriginNorm(m,n)
addpath('WZPUtil');
if nargin < 1
    m=128;
    n=128;
end

%%% params
controlPoints = 5;% Number of control points in each direction, each point corresponds to a random vector
shiftParam = 20; % Strength of distortion文章里alpha
strengthen = 5;% Emphasize the greater impact of proximity 文章里p
% 位置随机 归一化[0,1]
rowLocationRange = [0.3,0.7]; % Location of deformation center (0-1)
colLocationRange = [0.3,0.7]; % Location of deformation center (0-1)
mu = [randR(rowLocationRange)*m, randR(colLocationRange)*n]; % u1 u2

% Random Sigma
N = 2;
D = diag(rand(N,1)); % 二维随机对角阵
U = orth(rand(N,N)); % 另一个二维随机矩阵正交基
Sigma = U' * D * U; % 互相关矩阵：表示形变区域的形状和大小
Sigma = Sigma*((m+n)*2); % s为缩放因子 此处为(m+n)*2

[X1,X2] = meshgrid(1:m,1:n); % x1,x1
X = [X1(:) X2(:)];
p = mvnpdf(X, mu, Sigma); % 生成f(x)
origin = reshape(p,m,n);
originNorm = origin/max(max(origin)); % Normalization [0,1]
% 控制点随机向量场
[uxmat,uymat] = randxymat(controlPoints,controlPoints);%文章里Vx Vy

[m,n] = size(originNorm);

My = repmat(round(linspace(1,m,controlPoints))',1,controlPoints);
Mx = repmat(round(linspace(1,n,controlPoints)),controlPoints,1);

result = zeros(m,n);
maxDistance = sqrt((m-1)^2+(n-1)^2);

for i=1:m
    for j=1:n
        distance = sqrt((My-i).^2+(Mx-j).^2);% 控制点与每个点的距离
        affect = influence(distance,maxDistance,strengthen);%公式中M
        affect=affect*shiftParam;
        
        % Pixel offset
        xShift = sum(sum(uxmat.*affect));
        yShift = sum(sum(uymat.*affect));
        % 由于XY的偏移量不是整数，所以需要重采样
        result(i,j) = myInterp2(originNorm,i+yShift,j+xShift);
    end
end

end

%% Subfunction
function zz = myInterp2(A,y,x) % Interpolation
y1 = [fix(y) ceil(y)];
x1 = [fix(x) ceil(x)];

[m,n]=size(A);
y1 = limIndex(y1,m);
x1 = limIndex(x1,n);

v1 = A(y1,x1);

if(x1(1)==x1(2))
    xv1 = v1(1,1);
    xv2 = v1(2,1);
else
    xv1 = (v1(1,2)-v1(1,1))/(x1(2)-x1(1))*(x-x1(1))+v1(1,1);
    xv2 = (v1(2,2)-v1(2,1))/(x1(2)-x1(1))*(x-x1(1))+v1(2,1);
end

if(y1(1)==y1(2))
    zz = xv1;
else
    zz = (xv2-xv1)/(y1(2)-y1(1))*(y-y1(1))+xv1;
end
end

% Calculate the effect of the random vector from the distance
function u=influence(d,maxD,strengthen)
t = 1-d./maxD;
t = t.^strengthen; % Used to emphasize that the closer the distance, the greater the impact, and the farther the distance, the less the impact
u = -2*t.^3+3*t.^2; % The smoothing function is designed empirically
% and its derivative is zero at t = 0 and t = 1,
% which avoids the discontinuity of the analog deformation signal.
end

% Prevents exceeding the original image when indexing
function input = limIndex(input,maxV,minV)
if nargin<3
    minV=1;
end
input(input>maxV)=maxV;
input(input<minV)=minV;
end

% Generate random vectors with different methods
function [uxmat,uymat]=randxymat(numx,numy,method)
if nargin<3
    method=1; % method 1 is recommended
end
if method==1
    theta = randR([0,2*pi],numy,numx);
    r = rand(numy,numx);
    uxmat = r.*cos(theta);
    uymat = r.*sin(theta);
elseif method==2
    uxmat = randR([-1,1],numy,numx);
    uymat = randR([-1,1],numy,numx);
elseif method==3
    num=numx*numy;
    uxmat=zeros(numy,numx);
    uymat=zeros(numy,numx);
    % Generate random vectors and determine if they are within the unit circle. Less efficient.
    for j=1:num
        k=0;
        while k==0
            randxy=rand(1,2);
            if (randxy(1)-0.5)^2+(randxy(2)-0.5)^2<=0.25
                k=k+1;
                uxmat(j)=randxy(1);
                uymat(j)=randxy(2);
            end
        end
    end
    uxmat=(uxmat-0.5)*2; % Adjusted to [-1,1]
    uymat=(uymat-0.5)*2;
end
end


