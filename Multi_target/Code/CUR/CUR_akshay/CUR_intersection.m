clc;close all;clear all;
X  = imread('10-12813-01-2.bmp');
X = double(rgb2gray(X));
figure
imshow(uint8(X))
s = 100;
r = 100;
k= 50;
XX = CUR_v2(X,k,s,r);
figure
imshow(uint8(XX))