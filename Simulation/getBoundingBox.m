% ===========================================================
% Filename:     getBoundingBox.m
% Date:   	 	2025-08-14
% Author:    	Fukang Qiao
% Description:  Main function for drawing rectangle phase label
% ===========================================================

function [bbox,columnLabel]=getBoundingBox(BW,classId,showFlag)
if nargin<2
    classId=1;
end
if nargin<3
    showFlag=0;
end
if ~islogical(BW)%判断大梯度形变相位是否大于2pai
    error('Input must be a 2-D binary image.');
end
% 与YOLO.yaml格式一样
columnLabel = {'classId', 'xCenter', 'yCenter', 'width', 'height'};
bbox=[];

[m,n]=size(BW);
[L,Num] = bwlabel(BW);
if Num<1    
    return;
end
% 归一化xywh 0-1
for idx=1:Num
    [y,x]=find(L==idx);
    x1 = min(x); x2 = max(x);
    y1 = min(y); y2 = max(y);

    xCenter = (x1+x2)/2/n;
    yCenter = (y1+y2)/2/m;
    w = (x2-x1)/n;
    h = (y2-y1)/m;

    bbox = [bbox; classId xCenter yCenter w h];
end
if showFlag
    showBbox(BW,bbox);
end

end


