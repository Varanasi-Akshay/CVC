clc;close all;clear all;imtool close all;

addpath('/home/akshay/Desktop/CVC/Multi_target/Code/Multi-level_Sampling');
addpath('/home/akshay/Desktop/Data/Biopsy');
addpath('/home/akshay/Desktop/Data');

% X = imread('yolo004.png');
% X = imread('stain_thumbnail.png');
% X  = imread('26.jpg');
%X  = imread('45.jpg');
X = imread('10264_00056.tif');

% 
% % if grayscale
% X = double(rgb2gray(X));
% % figure
% % imshow(uint8(X))
% 
% X=X(3001:5000,3001:6000);
% figure
% imshow(uint8(X))

% if rgb
X = double(X);
% X=X(3001:5000,3001:6000,:);
%X = X(1:500,1:500,:);
figure
imshow(uint8(X))

%% Parameters
% patch_size = 100;
patch_size = 50;

s = 30;
r = 30;
k= 10;

% no.of patches
m = size(X)/patch_size; 
disp(m)


% patch_data = zeros([patch_size*patch_size,m(1)*m(2)]);
patch_data = zeros([patch_size*patch_size*3,m(1)*m(2)]);

%% Reconstructing each patch using CUR
count = 0;
Img_reconst = zeros(size(X));
for i=1:m(1)
    for j=1:m(2)
        patch = X(1+(i-1)*patch_size:(i)*patch_size,1+(j-1)*patch_size:(j)*patch_size,:);
        for l=1:3
            patch_cur(:,:,l) = CUR_v2(patch(:,:,l),k,s,r);
        end    
%         patch_cur = CUR_v2(patch,k,s,r);
      %  patch_cur = patch;
        XX(1+(i-1)*patch_size:(i)*patch_size,1+(j-1)*patch_size:(j)*patch_size,:) = patch_cur;
        count = count + 1;
        patch_data(:,count) = reshape(patch_cur,[patch_size*patch_size*3,1]);
        
    end    
end

figure
imshow(uint8(XX))
figure
v = scree_plot_new(patch(:,:,1));

%% SVD
%[U, S, V] = svd(patch_data);

figure
[U,Sig,V,v] = scree_plot_new(patch_data);
rank = 40;
Coeff = Sig(1:rank,1:rank)*V(:,1:rank)';

patch_data_svd = U(:,1:rank)*Coeff;

%%% For checking if SVD rank is correct or not

% count = 0;
% Img_reconst = zeros(size(X));
% for i=1:m(1)
%     for j=1:m(2)
%         count = count + 1;
%         patch_img = patch_data_svd(:,count);
%         Img_reconst(1+(i-1)*patch_size:(i)*patch_size,1+(j-1)*patch_size:(j)*patch_size,:) = reshape(patch_img,[patch_size,patch_size,3]);
%         
%     end    
% end
% figure
% imshow(uint8(Img_reconst))

%% Clustering
n = 10;
%klist= 3;
klist=2:n;%the number of clusters you want to try
% myfunc = @(X,K)(kmeans(X, K));
% eva = evalclusters(net.IW{1},myfunc,'CalinskiHarabasz','klist',klist)
% classes=kmeans(net.IW{1},eva.OptimalK);

%rng('default');  % For reproducibility 'CalinskiHarabasz' 'gap',
num_trials = 10;
temp = inf;
for i=1:num_trials
    state = rng(i);
%     seed = rnd
%     rng(seed);


    eva = evalclusters(patch_data','kmeans','silhouette','KList',klist);
    [classes_temp, centers_temp, sumd_temp, D_temp ]=kmeans(patch_data',eva.OptimalK);
    loss = sum(sumd_temp);
    if loss < temp
        temp = loss;
        state_save = state;
        classes = classes_temp; 
        centers = centers_temp;
        sumd = sumd_temp;
        D = D_temp;
        eva_save = eva;
    end    
end
figure
plot(eva_save)

% How many in each cluster
% A = accumarray(classes,1)

% To visualize the patches
count=0;
vis_patch = zeros(size(X,1),size(X,2));
for i=1:m(1)
    for j=1:m(2)
        count=count+1;
        vis_patch(1+(i-1)*patch_size:(i)*patch_size,1+(j-1)*patch_size:(j)*patch_size)=classes(count);
    end
end    
figure
imagesc(vis_patch)

% figure
% pcolor(vis_patch);

% take the patch closest to class centers
[M,I] = min(D);
data_close_center = zeros(patch_size*patch_size*3,eva.OptimalK);
for i=1:eva.OptimalK
    patch_img = patch_data(:,I(i));
    data_close_center(:,i) = patch_img;
    patch_img = reshape(patch_img,[patch_size, patch_size,3]); 
    figure
    imshow(uint8(patch_img))
end    

for i=1:eva.OptimalK
    img = zeros(size(X));
    Rimg = zeros(size(X,1),size(X,2));
    Gimg = zeros(size(X,1),size(X,2));
    Bimg = zeros(size(X,1),size(X,2));
    ind=find(vis_patch==i);
    R = X(:,:,1);
    G = X(:,:,2);
    B = X(:,:,3); 
    
    Rimg(ind) = R(ind);
    Gimg(ind) = G(ind);
    Bimg(ind) = B(ind);
    img(:,:,1) = reshape(Rimg,[size(X,1),size(X,2)]);
    img(:,:,2) = reshape(Gimg,[size(X,1),size(X,2)]);
    img(:,:,3) = reshape(Bimg,[size(X,1),size(X,2)]);
    figure
    imshow(uint8(img))
end    


%% Wavelet decomposition

% % convert into
% img_gray = rgb2gray(patch_img);
% img_ycbcr = rgb2ycbcr(patch_img);
% img_gray= img_ycbcr(:,:,1);
% 
% % 
% %%% MulitWavelet plot
% level = 2;
% [C,S]=wavedec2(img_gray,level,'haar');
% 
% dec=plotwavelet2(C,S,level,'haar',255,'square');
% % 
% N = [1];
% NC = wthcoef2('h',C,S,N);
% dec=plotwavelet2(NC,S,level,'haar',255,'square');
% NC = wthcoef2('v',NC,S,N);
% dec=plotwavelet2(NC,S,level,'haar',255,'square');
% NC = wthcoef2('d',NC,S,N);
% dec=plotwavelet2(NC,S,level,'haar',255,'square');
% 
% low_res = waverec2(NC,S,'haar');
% low_res = reshape(low_res,[patch_size,patch_size]);
% 
% imtool(low_res,[]);
% %imshow(low_res)
% 
% figure
% N = [2];
% NC = wthcoef2('h',NC,S,N);
% dec=plotwavelet2(NC,S,level,'haar',255,'square');
% NC = wthcoef2('v',NC,S,N);
% dec=plotwavelet2(NC,S,level,'haar',255,'square');
% NC = wthcoef2('d',NC,S,N);
% dec=plotwavelet2(NC,S,level,'haar',255,'square');
% low_res2 = waverec2(NC,S,'haar');
% low_res2 = reshape(low_res2,[patch_size,patch_size]);
% 
% imtool(low_res2,[]);


% for i=1:eva.OptimalK
%     img = data_close_center(:,i);
%     img = reshape(img,[patch_size, patch_size]);
%     for j = 1:
%         
%     end
% end
