% The demo below is for wavelets
clc;clear;close all;imtool close all;
%load('Shepp-Logan.mat') % test image, 1024 by 1024
addpath('/home/akshay/Desktop/CVC/Multi_target/Code/Multi-level_Sampling/spgl1-1.9')

img = imread('/home/akshay/Desktop/Data/KIMIAPath24RGB_zip/KIMIAPath24 RGB/Training/s0/2.jpg');
img_reshaped = imresize(img, [512, 512]);
img_double = double(img_reshaped);
I = rgb2gray(img_double);

%% Step 1: perform sampling
[n1,n2] = size(I);
level = 7;
per = 0.30;
N = n1*n2;
m = floor(per*N);
const = 1;
M = const*m;

c = 3; % 1 - uniform, 2 - Fourier+wavelet, 3 - Hadamard+wavelet
scheme = uniform_rect_samp_scheme(n1,n2,M);
sample_indices={};
% using samples of that level only
for i=2:length(scheme)
    temp_scheme=zeros(size(scheme));
    %temp_scheme(1)=4;
    temp_scheme(1:i)=scheme(1:i);
    pos = random_rect_subsamp(temp_scheme);
    IND = sub2ind([n1,n2],pos(:,1),pos(:,2));
    IND = sort(IND);
    sample_indices{i-1} = IND;
end
%% Step 1: through trail and error, find a sparse approximation of the image in the wavelet domain
    
% 
% [II,S] = sparse_approx(I,level,per);
%     PSNR21 = psnr(II,I);
% 
%     MulitWavelet plot
%     [C,S] = wavedec2(I,level,'haimg_double(:,:,1)ar');
% 
%     dec=plotwavelet2(C,S,level,'haar',255,'square');
% 
%     coeff=dec{1}/255;
% 
%     showing the sparsity
%     Finding no.of coeff at each level
%     max_level=log2(size(I,1));
%     coeff_level=zeros(level,1);
%     for i=1:level
%         coeff_count=nnz(coeff(1:2^(max_level-level+i),1:2^(max_level-level+i)));
%         coeff_count = coeff_count-sum(coeff_level);
%         coeff_level(i) = coeff_count;
%     end

img_reconstructed=zeros(size(img_double));
for i=1:3
    
    I = img_double(:,:,i)/255;

    [II,S] = sparse_approx(I,level,per);
    PSNR21 = psnr(II,I);
    IND = sample_indices{level};
    %[sample_level{i-1}]=sample_level_count(I,IND,level);
    b = fwht2(I);
    b = b(IND)*sqrt(n1*n2);
    bb = fwht2(II);
    bb = bb(IND)*sqrt(n1*n2);
    sigma = norm(bb-b);
    opA = @(x,mode) partial_Hadamard(x,S,IND,n1,n2,level,mode);
    C = spg_bpdn(opA,b,sigma);
    III = waverec2(C,S,'haar'); % Changed from haar to db4
    III = reshape(III,n1,n2);
    %III = rescale(III);
    img_reconstructed(:,:,i)=rescale(III*255,0,255);    
   
end    

image(uint8(img_reconstructed))
% %imtool(III)
%         des=num2str(i-1);
%         des=strcat('Kimia_upto ',des);
%         %imshow(III)
%         hFigure = imtool(III,[]);
%         saveas(hFigure,des,'png')
%         psnr_values(i-1) = psnr(III,I);
    %     set(hFigure,'NumberTitle','off','Name',des);    
%plot(psnr_values)
imtool close all


