% The demo below is for wavelets
clc;clear;close all;imtool close all;
%load('Shepp-Logan.mat') % test image, 1024 by 1024
addpath('/home/akshay/Desktop/CVC/Multi_target/Code/Multi-level_Sampling/spgl1-1.9')

img = imread('10-12813-01-2.bmp');
img_gray = rgb2gray(img);
I_full = double(img_gray);
I = I_full(1:512,1:512)/255;

%% Step 1: through trail and error, find a sparse approximation of the image in the wavelet domain
level = 8;
per = 0.03;
[II,S] = sparse_approx(I,level,per);
PSNR21 = psnr(II,I);


% MulitWavelet plot
[C,S] = wavedec2(I,level,'haar');

%% Step 2: perform sampling
[n1,n2] = size(I);
N = n1*n2;
m = floor(per*N);
const = 3;
M = const*m;

c = 3; % 1 - uniform, 2 - Fourier+wavelet, 3 - Hadamard+wavelet
scheme = uniform_rect_samp_scheme(n1,n2,M);

temp_scheme=zeros(size(scheme));
temp_scheme(1:9)=scheme(1:9);
%temp_scheme(8)= temp_scheme(8) + scheme(9)+10000;
pos = random_rect_subsamp(temp_scheme);
IND = sub2ind([n1,n2],pos(:,1),pos(:,2));
IND = sort(IND);

Sample=zeros(size(I));
Sample(IND)=1;




sample_level=sample_level_count(I,IND,level);
b = fwht2(I);
b = b(IND)*sqrt(n1*n2);
bb = fwht2(II);
bb = bb(IND)*sqrt(n1*n2);
sigma = norm(bb-b);
opA = @(x,mode) partial_Hadamard(x,S,IND,n1,n2,level,mode);
C = spg_bpdn(opA,b,sigma);
III = waverec2(C,S,'haar'); % Changed from haar to db4
III = reshape(III,n1,n2);
psnr_values = psnr(III,I);    
