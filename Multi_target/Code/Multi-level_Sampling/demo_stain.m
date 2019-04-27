% The demo below is for wavelets
clc;clear;close all
%load('Shepp-Logan.mat') % test image, 1024 by 1024
addpath('/home/akshay/Downloads/Multitarget/Codes/Multi-level_Sampling/spgl1-1.9')

img = imread('10-12813-01-2.bmp');
img_gray = rgb2gray(img);
I_full = double(img_gray);
I = I_full(1:512,1:512)/255;
%I = I_512;


%for i=1:8
   
%% Step 1: through trail and error, find a sparse approximation of the image in the wavelet domain
level = 2;
per = 0.03;
[II,S] = sparse_approx(I,level,per);
PSNR21 = psnr(II,I); % model error

%% Step 2: perform sampling
[n1,n2] = size(I);
N = n1*n2;
m = floor(per*N);
const = 6;
M = const*m;

c = 3; % 1 - uniform, 2 - Fourier+wavelet, 3 - Hadamard+wavelet

switch c
    case 1
        % uniform sampling
        IND = randsample(n1*n2,M); % uniform sampling
        IND = sort(IND);
    case 2
        % multi-level sampling (Fourier+wavelet)
        nu = round(log2(min(n1,n2)));
        bounds = [2^(nu-4),2^(nu-2),2^(nu-1)];
        pos = random_circ_subsamp(n1,n2,bounds,M);
        IND = sub2ind([n1,n2],pos(:,1),pos(:,2));
        IND = sort(IND);
    case 3
        % multi-level sampling (Hadamard+wavelet)
        scheme = uniform_rect_samp_scheme(n1,n2,M);
        pos = random_rect_subsamp(scheme);
        IND = sub2ind([n1,n2],pos(:,1),pos(:,2));
        IND = sort(IND);
end



%% Step 3: recovery from samples

switch c
    case 1
        b = I(IND);
        bb = II(IND);
        sigma = norm(II(IND)-I(IND));
        opA = @(x,mode) partial_Uniform(x,S,IND,n1,n2,level,mode);
    case 2 
        b = fftshift(fft2(I))/sqrt(n1*n2);
        b = b(IND);
        bb = fftshift(fft2(II))/sqrt(n1*n2);
        bb = bb(IND);
        sigma = norm(bb-b);
        opA = @(x,mode) partial_Fourier(x,S,IND,n1,n2,level,mode);
    case 3
        b = fwht2(I);
        b = b(IND)*sqrt(n1*n2);
        bb = fwht2(II);
        bb = bb(IND)*sqrt(n1*n2);
        sigma = norm(bb-b);
        opA = @(x,mode) partial_Hadamard(x,S,IND,n1,n2,level,mode);
end

C = spg_bpdn(opA,b,sigma);
III = waverec2(C,S,'haar'); % Changed from haar to db4
III = reshape(III,n1,n2);

PSNR32 = psnr(III,II); % recovery error
PSNR31 = psnr(III,I); 

% %%% 2-norm error
% error21 = norm(I-II);
% error23 = norm(II-III);
error13 = norm(I-III)/norm(I)*100;

%%% Frobenius-norm error
error13_fro = norm(I-III,'Fro')/norm(I,'Fro')*100;

%%% inf-norm error
error13_inf = norm(I-III,'inf')/norm(I,'inf')*100;


%%% If value is scaled between 0-255 
% figure(1)
% imshow(uint8(I))
% figure(2)
% imshow(uint8(II))
% figure(3)
% imshow(uint8(III))

%%% If value is scaled between 0-1
imtool(I)
imtool(II)
imtool(III)

