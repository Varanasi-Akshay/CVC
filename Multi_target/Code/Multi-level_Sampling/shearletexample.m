%% Add path
clear all;clc;
addpath('/home/akshay/Downloads/Multitarget/Codes/Multi-level_Sampling/ShearLab3Dv11/Util')
addpath('/home/akshay/Downloads/Multitarget/Codes/Multi-level_Sampling/ShearLab3Dv11/2D')


img = imread('10-12813-01-2.bmp');

%sizeX = 500;
%sizeY = 300;
useGPU = 0;

sigma = 30;
scales = 1;
thresholdingFactors = [0 3 3 4 4];

%data = randn(sizeX, sizeY);
data = rgb2gray(img);
data = double(data);
[sizeX, sizeY] = size(data);

% create shearlet system 
system = SLgetShearletSystem2D(useGPU,sizeX,sizeY,scales);

% co efficients
coeff = SLsheardec2D(data,system);

% thresholding
for j =1:system.nShearlets
   shearletIdx = system.shearletIdxs(j,:);
   %disp(shearletIdx);
   coeff(:,:,j) = coeff(:,:,j).*(abs(coeff(:,:,j)) > ...
   thresholdingFactors(shearletIdx(2)+1)*system.RMS(j)*sigma);
end
% Reconstruction
reconstruction = SLshearrec2D(coeff,system);
PSNR = SLcomputePSNR(data, reconstruction);