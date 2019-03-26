%% For making movie
workingDir = pwd;
%name='Semal_rgb.avi';
% name='Tingqua_rgb.avi';
% name='Lidar_rank.avi';
% 
% % name='Tingqua_rank_difference.avi';
% %name='Semal_rank_difference.avi';
% 
% 
% 

%% Tensor toolbox path
addpath(genpath('~/tensor_toolbox-master/'))

%% Reading the image

% Semal
%Img=imread('Semal_(Bombax_ceiba)_flowers_in_Kolkata_W_IMG_4132.jpg');

% LIDAR
%Img=imread('laser_technologies_totalsight_denver_1200x399-100K.jpg');

% % Tingqua
 Img=imread('Shop_of_Tingqua,_the_painter.jpg');
% 
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

% SVD of that
[u,s,v]=svd(Img_reshape,'econ');


% Img_recon=u(:,1:2)*s(1:2,1:2)*v(:,1:2)';
% Img_recon_shape=reshape(Img_recon,x,y,z);
% image(uint8(Img_recon_shape))



%% For finding the rank of individual bands

% [U1,S1,V1]=svd(Img(:,:,1));
% [U2,S2,V2]=svd(Img(:,:,2));
% [U3,S3,V3]=svd(Img(:,:,3));

%% Finding the low rank of each band using matrix,cp and tucker

% for i=1:20
%     rank=i;
%     % Matrix rank    
%     Img_matrix_rank(:,:,1)=U1(:,1:rank)*S1(1:rank,1:rank)*V1(:,1:rank)';    
%     Img_matrix_rank(:,:,2)=U2(:,1:rank)*S2(1:rank,1:rank)*V2(:,1:rank)';    
%     Img_matrix_rank(:,:,3)=U3(:,1:rank)*S3(1:rank,1:rank)*V3(:,1:rank)';
%     filename=sprintf('Semal_matrix_rank_%d.jpg',rank);
%     imwrite(uint8(Img_matrix_rank),filename) ; 
%     [peaksnr_matrix(i), snr_matrix(i)] = psnr(Img_matrix_rank, ref);
%     
%      % cp
%     M = cp_als(tensor(Img),rank,'maxiters',999);
%     A = M.U{1};
%     B = M.U{2};
%     C = M.U{3};
%     v = M.lambda;
%     I2 = ktensor(v,{A,B,C},[1 2 3]);
%     I2 = double(I2);
%     filename=sprintf('Semal_cp_rank_%d.jpg',rank);
%     imwrite(uint8(I2),filename) ; 
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
%     filename=sprintf('Semal_tucker_rank_%d.jpg',rank);
%     imwrite(uint8(I3),filename) ; 
%     [peaksnr_tucker(i), snr_tucker(i)] = psnr(I3, ref);
% end

        

