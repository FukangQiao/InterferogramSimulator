% ===========================================================
% Filename:     getLocalVariance.m
% Date:   	 	2025-08-14
% Author:    	Fukang Qiao
% Description:  function for calculation of variance and mean
% ===========================================================


function [dataStd, absDataMean] = getLocalVariance(srcData, r)
% 计算局部标准差和均值
[nHeight, nWidth] = size(srcData);
N               = boxfilter(ones(nHeight, nWidth), r);
absDataMean     = boxfilter(abs(srcData), r*2) ./ N;
srcDataMean     = boxfilter(srcData, r) ./ N;
srcDataDataMean = boxfilter(srcData.*srcData, r) ./ N;
dataVariance    = srcDataDataMean - srcDataMean .* srcDataMean;
dataVariance(dataVariance<0)=0;
dataStd         = sqrt(dataVariance);
end
