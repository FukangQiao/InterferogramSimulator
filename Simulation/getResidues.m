% ===========================================================
% Filename:     getResidues.m
% Date:   	 	2025-08-14
% Author:    	Fukang Qiao
% Description:  Calculate residual points
% ===========================================================

%getResidues Calculate residual points
% [a,b] = getResidues(w, plotFlag) 
% Returns the [row, column, value] of negative (a) and positive (b) residual points
% 
% residues = getResidues(w, plotFlag) 
% Returns the residual points in a matrix, which size is the same as w

function [a,b] = getResidues(w, plotFlag)
if nargin<1
    help getResidues
    return
end
if nargin<2
    plotFlag = 0;
end

Infinitesimal = 1e-9;

[M,N]=size(w);
m=M-1;n=N-1;

tmp11 = w(1:m,1:n);
tmp12 = w(1:m,2:N);
tmp21 = w(2:M,1:n);
tmp22 = w(2:M,2:N);

errorP=wrapToPi(tmp12-tmp11)+...
    wrapToPi(tmp22-tmp12)+...
    wrapToPi(tmp21-tmp22)+...
    wrapToPi(tmp11-tmp21);

a = [[errorP zeros(m,1)]; zeros(1,N)];

if nargout==2
    [x,y] = find(errorP<-Infinitesimal);
    a = [x+0.5,y+0.5,errorP(sub2ind(size(errorP),x,y))];
    [x,y] = find(errorP>Infinitesimal);
    b = [x+0.5,y+0.5,errorP(sub2ind(size(errorP),x,y))];
end

if plotFlag
    if ~exist('b','var')
        [x,y] = find(errorP<-Infinitesimal);
        a = [x+0.5,y+0.5,errorP(sub2ind(size(errorP),x,y))];
        [x,y] = find(errorP>Infinitesimal);
        b = [x+0.5,y+0.5,errorP(sub2ind(size(errorP),x,y))];
    end
    hold on;scatter(a(:,2),a(:,1),20,'b','filled','Marker','o','MarkerEdgeColor','k');
    hold on;scatter(b(:,2),b(:,1),20,'r','filled','Marker','o','MarkerEdgeColor','k');
end