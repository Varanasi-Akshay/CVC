% The demo below is for wavelets
clc;clear;close all;imtool close all;


addpath('/home/akshay/Desktop/CVC/Multi_target/Code/Multi-level_Sampling/spgl1-1.9')

img = imread('/home/akshay/Desktop/Data/KIMIAPath24RGB_zip/KIMIAPath24 RGB/Training/s1/95.jpg');
img_reshaped = imresize(img, [512, 512]);
img_double = double(img_reshaped);
I = rgb2gray(img_double);

%% Step 1: perform sampling
[n1,n2] = size(I);
level = 6;
per = 0.30;
N = n1*n2;
m = floor(per*N);
const = 1;
M = const*m;

% Hadamard + Haar
scheme = uniform_rect_samp_scheme(n1,n2,M);
sample_indices={};
% using samples of that level only
for i=2:length(scheme) % since 1st value in schema should be 4 always
    temp_scheme=zeros(size(scheme));
    temp_scheme(1:i)=scheme(1:i);
    pos = random_rect_subsamp(temp_scheme);
    IND = sub2ind([n1,n2],pos(:,1),pos(:,2));
    IND = sort(IND);
    sample_indices{i-1} = IND; % saving the indices
end

% Testing for all the samples till level 8 
img_reconstructed=zeros(size(img_double));
for i=1:3
    
    I = img_double(:,:,i)/255;

    [II,S] = sparse_approx(I,level,per);
    PSNR21 = psnr(II,I);
    IND = sample_indices{level}; % samples till that level
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

% original image
figure()
image(uint8(img_double))


% reconstructed image
figure()
image(uint8(img_reconstructed))


