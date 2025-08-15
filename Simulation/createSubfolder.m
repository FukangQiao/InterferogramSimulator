% ===========================================================
% Filename:     createSubfolder.m
% Date:   	 	2025-08-14
% Author:    	Fukang Qiao
% Description:  Main function for creating sub folder
% ===========================================================

function folderPath=createSubfolder(rootPath,subfolder,createFlag)
if nargin<2
    subfolder='';
end
if nargin<3
    createFlag=1;
end
if ~createFlag
    folderPath = '';
    return;
end
folderPath = fullfile(rootPath,subfolder);
if ~exist(folderPath,'dir')
    mkdir(folderPath); 
end
end