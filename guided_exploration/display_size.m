clc;clear all;close all;
% c = checkerboard;

% c = imread('26.jpg');
c = imread('/home/akshay/Desktop/Data/Biopsy/stain_thumbnail.png');
% figure
% imshow(c,'InitialMagnification','fit')
% truesize
figure
imshow(c,'InitialMagnification','fit')
truesize([1000 1000]);
y = c(2000:2500,2000:2500,:);
figure
imshow(y,'InitialMagnification','fit')
truesize([1000 1000]);

z = y(200:250,200:250,:);
figure
imshow(z,'InitialMagnification','fit')
truesize([1000 1000]);


a = imread('26.jpg');
figure
imshow(a,'InitialMagnification','fit')
truesize
%truesize([1000 1000]);