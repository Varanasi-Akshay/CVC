%% For making movie
workingDir = pwd;

%% Tensor toolbox path
addpath(genpath('~/tensor_toolbox-master/'))
addpath(genpath('/home/akshay/Downloads/Multitarget/Codes/functions'))



%% Reading the image

Img= imread('10-12813-01-2.bmp');


% % Downsampling/resizing
%Img=imresize(Img,0.2);

% Making it double for further processing
Img=double(Img);

% Reference For metrics
ref= Img;


% Size
[x,y,z]=size(Img);

% Reshaping
Img_reshape=reshape(Img,x*y,z);

% % Normalizing
% Img_reshape= Img_reshape/255;


%% Finding the low rank using cp and tucker

% for i=1:20
%     rank=i;
%         
%      % cp
%     M = cp_als(tensor(Img),rank,'maxiters',999);
%     A = M.U{1};
%     B = M.U{2};
%     C = M.U{3};
%     v = M.lambda;
%     I2 = ktensor(v,{A,B,C},[1 2 3]);
%     I2 = double(I2);
% %     filename=sprintf('cp_%d.jpg',rank);
% %     imwrite(uint8(I2),filename) ; 
%     [peaksnr_cp(i), snr_cp(i)] = psnr(I2, ref);
%     
%     
%     % tucker
%     M = tucker_als(tensor(Img),[rank rank 3],'maxiters',999);
%     G = M.core;
%     A = M.U{1};
%     B = M.U{2};
%     C = M.U{3};    
%     I3 = ttm(G,{A,B,C},[1 2 3]);
%     I3 = double(I3);
% %     filename=sprintf('tucker_%d.jpg',rank);
% %     imwrite(uint8(I3),filename) ; 
%     [peaksnr_tucker(i), snr_tucker(i)] = psnr(I3, ref);
% end
% 
% close all;
% clear;
% clc;


%% Adding the paths

% addpath(genpath('~/Desktop/tensor_toolbox-master/'))
% addpath(genpath('~/Desktop/Denoising_MNF/'))
% addpath(genpath('~/Desktop/Denoising_MNF/functions'))
% addpath(genpath('~/Desktop/Denoising_MNF/better_faster_mnf'))
% addpath(genpath('~/Desktop/Denoising_MNF/better_noise_estimation'))



%% CP decomposition
start = 20; 
step = 20;
ending = 660;
iter = 0;
CP_results = zeros(1+(ending-start)/step,13);
for cp_rank=start:step:ending
       T=cp_als(tensor(C_denoised),cp_rank);
       U1 = T.U{1};
       U2 = T.U{2};
       U3 = T.U{3};
       d= T.lambda;
       A=ttm(tendiag(d.',[cp_rank cp_rank cp_rank]),{U1 U2 U3},[1 2 3]);
       A=double(A);
       [psnr1,Rmse1, Ergas1, Sam1, Uiqi1,Ssim1,DD01,CC01,pfe1,snr1,mb] = quality_assessment(C_denoised,A, 0, 1.0);
       iter = 1+iter;%(cp_rank-start)/step;
       CP_results(iter,1:13)=[iter,cp_rank,psnr1,Rmse1, Ergas1, Sam1, Uiqi1,Ssim1,DD01,CC01,pfe1,snr1,mb];
end
save('CP','CP_results')



%% Tucker decoposition
start1 = 10;
start2 = 10; 
step = 10;
ending1 = 10;
ending2 = 40; 

iter=0;
%%Spectral Rank= 9
Tucker_results_9 = zeros((1+(ending1-start1)/step)*(1+(ending2-start2)/step),14);
for tucker_rank1=start1:step:ending1
    for tucker_rank2=start2:step:ending2    
       T=tucker_als(tensor(C_denoised),[tucker_rank1,tucker_rank2,9]);
       C = T.core;
       W = T.U{1};
       H = T.U{2};
       S = T.U{3};
       A=ttm(C,{W H S},[1 2 3]);
       A=double(A);
       [psnr1,Rmse1, Ergas1, Sam1, Uiqi1,Ssim1,DD01,CC01,pfe1,snr1,mb] = quality_assessment(C_denoised,A, 0, 1.0);
      
       iter = iter + 1;
       Tucker_results_9(iter,1:14)=[iter,tucker_rank1,tucker_rank2,psnr1,Rmse1, Ergas1, Sam1, Uiqi1,Ssim1,DD01,CC01,pfe1,snr1,mb];
       
    end   
end


save('Tucker','Tucker_results_9')




