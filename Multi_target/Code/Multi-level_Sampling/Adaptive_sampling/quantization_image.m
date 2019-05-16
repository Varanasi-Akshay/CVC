clc; clear all;

% same image only just the name is different
img = imread('tubule_041.bmp');

% annotation of that
annotation = imread('tubule_041_anno.bmp');
annotation_double = double(annotation);

%% Quantization

img_double = double(rgb2gray(img));


% Since you want 32 bit 256/32
img_double = floor(img_double/64);

% to view it we divide by 32
img_final = (img_double)*64;

imtool(uint8(img_final))
%imtool(img_final)