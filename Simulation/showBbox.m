% ===========================================================
% Filename:     showBbox.m
% Date:   	 	2025-08-14
% Author:    	Fukang Qiao
% Description:  Show rectangle labels information 
% ===========================================================

function showBbox(img,contents)
[m,n]=size(img);
imagesc(img);
for i=1:size(contents,1)
    content = contents(i,:);
    labelId = content(1);
    w = content(4);
    h = content(5);
    x = content(2)-w/2;
    y = content(3)-h/2;

    rectangle('Position',[x*n,y*m,w*n,h*m],'EdgeColor','r','LineWidth',2);
    text(x*n,y*m,num2str(labelId),'BackgroundColor','w');
end

end