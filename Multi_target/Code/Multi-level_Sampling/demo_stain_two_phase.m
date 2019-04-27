% The demo below is for wavelets
clc;clear;close all
%load('Shepp-Logan.mat') % test image, 1024 by 1024
addpath('/home/akshay/Downloads/Multitarget/Codes/Multi-level_Sampling/spgl1-1.9')
addpath('/home/akshay/Downloads/Multitarget/Codes/local_coherences')
img = imread('10-12813-01-2.bmp');
img_gray = rgb2gray(img);
I_full = double(img_gray);
I = I_full(1:512,1:512)/255;
%I = I_512;

%r = 100;

%for i=1:8
   
%% Step 1: through trail and error, find a sparse approximation of the image in the wavelet domain
level = 5;
per = 0.03;
[II,S] = sparse_approx(I,level,per);
PSNR21 = psnr(II,I); % model error

%% Step 2.a: perform sampling first multi-level
tic;
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
        IND1 = sort(IND);
    case 2
        % multi-level sampling (Fourier+wavelet)
        nu = round(log2(min(n1,n2)));
        bounds = [2^(nu-4),2^(nu-2),2^(nu-1)];
        pos = random_circ_subsamp(n1,n2,bounds,M);
        IND = sub2ind([n1,n2],pos(:,1),pos(:,2));
        IND1 = sort(IND);
        Z =fftshift(fft2(I))/sqrt(n1*n2);
    case 3
        % multi-level sampling (Hadamard+wavelet)
        scheme = uniform_rect_samp_scheme(n1,n2,M);
        pos = random_rect_subsamp(scheme);
        IND = sub2ind([n1,n2],pos(:,1),pos(:,2));
        IND1 = sort(IND);
        Z=fwht2(I)*sqrt(n1*n2);
end

X=zeros(n1,n2);
X(IND1)=Z(IND1);

imtool(X)


%% Step 2.b: Perform sampling using coherency sampling.
% beta = 0.6;
% 


r = 60;
%Step 2: sampling based on local coherences, needs r value

[U,Sing,V] = svd(X);
s = diag(Sing);


U = U(:,1:r);
V = V(:,1:r);

row = zeros(n1,1);
for i = 1:n1
    row(i) = norm(U(i,1:r))^2*n1/r;
end

col = zeros(n2,1);
for j = 1:n2
    col(j) = norm(V(j,1:r))^2*n2/r;
end
% 
%M = round(0.05*n1*n2); % Step 3 needs your input of sample size
% 
% % Step3: local coherences sampling
[IND2,~] = local_coherences_v2(row,col,M);
IND2 = unique(IND2);
% 
Y = zeros(n1,n2);
Y(IND2)=Z(IND2);
imtool(Z)

IND=[IND1 ; IND2];
% actual_sample_size = length(IND);
fprintf("Time taken for sampling %d \n",toc)


    

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
imtool(III)

%%% If value is scaled between 0-255 
% figure(1)
% imshow(uint8(I))
% figure(2)
% imshow(uint8(II))
% figure(3)
% imshow(uint8(III))

%%% If value is scaled between 0-1
% imtool(I)
% imtool(II)
%imtool(III)

