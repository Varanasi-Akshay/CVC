%% For making movie
workingDir = pwd;

%% Tensor toolbox path
addpath(genpath('~/tensor_toolbox-master/'))
addpath(genpath('/home/akshay/Downloads/Multitarget/Codes/functions/'))

%% Reading the image

Img= imread('stains_775x522/10-12813-01-2.bmp');


% % Downsampling/resizing
%Img=imresize(Img,0.2);

% Making it double for further processing
img=(Img);
color_dim = 1024;
img_ind = rgb2ind(img, color_dim);
M1 = double(img_ind); 
mw = size(M1,1);
mh = size(M1,2);
    
%block_size = round(linspace(1,min(mh,mw)/2,100));
block_size = 25;    
blocks_reshaped = matrix2tensor(M1,block_size);

% Reference For metrics
ref= blocks_reshaped;

%% CP decomposition
start = 10; 
step = 5;
ending = block_size;
iter = 0;
CP_results = zeros(1+(ending-start)/step,14);
for cp_rank=start:step:ending
       T=cp_als(tensor(ref),cp_rank);
       U1 = T.U{1};
       U2 = T.U{2};
       U3 = T.U{3};
       d= T.lambda;
       A=ttm(tendiag(d.',[cp_rank cp_rank cp_rank]),{U1 U2 U3},[1 2 3]);
       A=double(A);
       [psnr1,Rmse1, Ergas1, Sam1, Uiqi1,Ssim1,DD01,CC01,pfe1,snr1,mb] = quality_assessment(ref,A, 0, 1.0);
       iter = 1+iter;%(cp_rank-start)/step;
       CP_results(iter,1:13)=[iter,cp_rank,psnr1,Rmse1, Ergas1, Sam1, Uiqi1,Ssim1,DD01,CC01,pfe1,snr1,mb];
end
save('CP','CP_results')



% %% Tucker decoposition
% start1 = 10;
% start2 = 10; 
% step = 10;
% ending1 = block_size;
% ending2 = block_size; 
% 
% iter=0;
% %%Spectral Rank= 9
% Tucker_results_9 = zeros((1+(ending1-start1)/step)*(1+(ending2-start2)/step),15);
% for tucker_rank1=start1:step:ending1
%     for tucker_rank2=start2:step:ending2    
%        T=tucker_als(tensor(ref),[tucker_rank1,tucker_rank2,9]);
%        C = T.core;
%        W = T.U{1};
%        H = T.U{2};
%        S = T.U{3};
%        A=ttm(C,{W H S},[1 2 3]);
%        A=double(A);
%        [psnr1,Rmse1, Ergas1, Sam1, Uiqi1,Ssim1,DD01,CC01,pfe1,snr1,mb] = quality_assessment(ref,A, 0, 1.0);
%       
%        iter = iter + 1;
%        Tucker_results_9(iter,1:14)=[iter,tucker_rank1,tucker_rank2,psnr1,Rmse1, Ergas1, Sam1, Uiqi1,Ssim1,DD01,CC01,pfe1,snr1,mb];
%        
%     end   
% end
% 
% 
% save('Tucker_Pavia','Tucker_results_9')




