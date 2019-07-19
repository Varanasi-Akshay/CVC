% We need to show labels
% I = imread('26.jpg');
% figure
% imshow(I)
% h = drawrectangle('Position',[250,250,300,300],'StripeColor','r');

%%
% NumberArrays=[20 20 120 120];
% %Show the Image
% figure(1), imshow('26.jpg');
% h = imrect(gca, NumberArrays);
% addNewPositionCallback(h,@(p) title(mat2str(p,3)));
% %recall function to draw the rectangle. 
% fcn = makeConstrainToRectFcn('imrect',get(gca,'XLim'),get(gca,'YLim'));
% setPositionConstraintFcn(h,fcn); 

%%
clc; clear all;close all;
img = imread('26.jpg');
fh = figure;
imshow( img, 'border', 'tight' ); %//show your image
hold on;
rectangle('Position', [900 750 50 50], 'EdgeColor', 'r'); %// draw rectangle on image
frm = getframe( fh ); %// get the image+rectangle
imwrite( frm.cdata, 'savedFileName.png' ); %// save to file
