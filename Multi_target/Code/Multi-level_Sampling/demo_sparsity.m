% The demo below is for wavelets
clc;clear;close all
%load('Shepp-Logan.mat') % test image, 1024 by 1024
addpath('/home/akshay/Desktop/CVC/Multi_target/Code/Multi-level_Sampling/spgl1-1.9')

img = imread('10-12813-01-2.bmp');
img_gray = rgb2gray(img);
I_full = double(img_gray);
I = I_full(1:512,1:512)/255;


%% Step 1: through trail and error, find a sparse approximation of the image in the wavelet domain
level = 3;
per = 0.03;
[II,S] = sparse_approx(I,level,per);
PSNR21 = psnr(II,I); % model error
wave_name='haar';
X=I;
for p=0:0.1:1
    [im_sparse, wave_coef,S]=wavelet_sparsify2d(X, p, wave_name)
end    