% ===========================================================
% Filename:     getSlop.m
% Date:   	 	2025-08-14
% Author:    	Fukang Qiao
% Description:  Simulate Slope phase
% ===========================================================

function img = getSlop(ro,co,ampRange)

%%% params
if nargin<3
    ampRange = [0,co/4.*pi]; % Amplitude range
end
rows = 1.5*ro;
cols = 1.5*co;

cx = 1:cols;
cy = func(cx);

rx = 1:rows;
ry = func(rx);

img = ry'*cy;
img = img./(max2(img)-min2(img));
img = img.*randR(ampRange);

img = imrotate(img,randR([0,360]),'bilinear');
% subplot(121),imagesc(wrapToPi(img));colormap jet;
[rows,cols]=size(img);
r1 = round((rows-ro)/2);
c1 = round((cols-co)/2);
img = img(r1:r1+ro-1,c1:c1+co-1);

% subplot(122),imagesc(wrapToPi(img));colormap jet;
end

function res=func(x)
if rand<0.2 % Probability of nonlinear transformation
    res= x.^2.*randn;
elseif rand<0.4 % Probability of nonlinear transformation
    res= sqrt(x).*randn;
else
    res = x*randn;
end
end