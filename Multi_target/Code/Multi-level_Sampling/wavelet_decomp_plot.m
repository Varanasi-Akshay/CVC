% n = 256;
% M = load_image('lena');
% % extract a sub-part of it
% M = rescale(crop(M,n));

addpath('/home/akshay/Desktop/CVC/Multi_target/Code/Multi-level_Sampling/spgl1-1.9')

img = imread('10-12813-01-2.bmp');
img_gray = rgb2gray(img);
I_full = double(img_gray);
I = I_full(1:512,1:512)/255;

M=I;

% options.wavelet_type = 'daubechies'; % kind of wavelet
% options.wavelet_vm = 4; % number of vanishing moments
% Jmin = 3; %  minimum scale


options.wavelet_type = 'daubechies'; % kind of wavelet
options.wavelet_vm = 1; % number of vanishing moments
Jmin = 3; %  minimum scale

% Wavelet
MW = perform_wavelet_transform(M, Jmin, +1, options);

% Recovered
M1 = perform_wavelet_transform(MW, Jmin, -1, options);
disp(['Error of recovery (should be close to 0): ' num2str(norm(M-M1, 'fro')/norm(M, 'fro')) '.'])

subplot(1,2,1);
imageplot(M); title('Image');
subplot(1,2,2);
plot_wavelet(MW,Jmin); title('Wavelet coefficients');