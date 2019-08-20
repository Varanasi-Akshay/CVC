clc;close all;clear all;imtool close all;

addpath('/home/akshay/Desktop/CVC/Multi_target/Code/Multi-level_Sampling');
addpath('/home/akshay/Desktop/Data/Biopsy');
addpath('/home/akshay/Desktop/Data');
X  = imread('26.jpg');

% if rgb
X = double(X);
X = X(1:960,1:960,:);

figure
imshow(uint8(X))
% X=X(3001:5000,3001:6000,:);
%X = X(1:300,01:300,:);
% figure
% imshow(uint8(X))

%% Parameters
% patch_size = 100;
%patch_size = 50;

patch_size = 64;
basis = 'haar';
max_level= wmaxlev([patch_size,patch_size],basis);


% no.of patches
m = size(X)/patch_size; 
% disp(m)


% patch_data = zeros([patch_size*patch_size,m(1)*m(2)]);
patch_data = zeros([patch_size*patch_size*3,m(1)*m(2)]);

% wavelet pixel data
% pixel_data = zeros(size(X,1),size(X,2),3*(1+3*max_level));

count = 0;
Img_reconst = zeros(size(X));
actual_zoom = zeros(m(1)*m(2),1);

for i=1:m(1)
    for j=1:m(2)
        patch = X(1+(i-1)*patch_size:(i)*patch_size,1+(j-1)*patch_size:(j)*patch_size,:);
        count = count + 1;
        patch_data(:,count) = reshape(patch,[patch_size*patch_size*3,1]);
        
        %%% 2d wavelet
        [C,S,actual_zoom(count)] = zoom_level(patch,basis);
%         [Ct,~] = wavelet_rgb(patch,basis);
        patch_wv_data(:,count) = C;%Ct

        
        %%% 3d wavelet
%         [Ct,actual_zoom(count)] = zoom_level_3d(patch,basis);
%         patch_wv_data(:,count) = Ct;
        
        %%% wavelet_pixel
%         for k=1:3
%             pixel_data(1+(i-1)*patch_size:(i)*patch_size,1+(j-1)*patch_size:(j)*patch_size,1+(k-1)*(1+3*max_level):k*(1+3*max_level)) = wavelet_pixel(patch(:,:,k),basis);        
%         end    
        
    end    
end


%% Clustering

rng('default');  % For reproducibility 'CalinskiHarabasz' 'gap','silhouette'
% num_trials = 10;
num_trials = 10;
temp = inf;

%% Finding how many components or clusters 



% till which level should I consider, 
t=1;
level_consider=1;
for i=1:level_consider
    t = t + 3*(2^(2*(i-1)));
end

% Original data
data = patch_wv_data(1,:)';

numComponents = 2;


%% Doing K-means


for i=1:num_trials
    state = rng(i);

%% Clustering using Kmeans
    [classes_temp, centers_temp, sumd_temp, D_temp ]=kmeans(data,numComponents);

    
    loss = sum(sumd_temp);
    if loss < temp
        temp = loss;
        state_save = state;
        classes = classes_temp; 
        centers = centers_temp;
        sumd = sumd_temp;
        D = D_temp;
    end  

% %% Clustering using GMM
%     eva = evalclusters(data,'gmdistribution','gap','KList',klist);
%     gmm = fitgmdist(data,eva.OptimalK,'SharedCovariance',true);
%     [classes_temp, nlogL_temp, P_temp,logpdf_temp, D_temp ]=cluster(gmm,data);   
%     state_save = state;
%     classes = classes_temp; 
%     D = D_temp;
%     eva_save = eva;
%     loss = nlogL_temp;
%     if loss < temp
%         temp = loss;
%         state_save = state;
%         classes = classes_temp; 
%         D = D_temp;
%         eva_save = eva;
%     end 

 
end


%% To visualize the patches with different clusters
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
colormap(jet(numComponents))
colorbar


% take the patch closest to class centers
[M,I] = min(D);
data_close_center = zeros(patch_size*patch_size*3,numComponents);
level_zoom = zeros(numComponents,1);
for i=1:numComponents
    patch_img = patch_data(:,I(i));
    data_close_center(:,i) = patch_img;
    patch_img = reshape(patch_img,[patch_size, patch_size,3]); 
    basis = 'haar';
    % 2D
    [C,S,level_zoom(i)] = zoom_level(patch_img,basis);
    
    %3D
%     [Ct,level_zoom(i)]=zoom_level_3d(patch_img,basis);
%     figure
%     imshow(uint8(patch_img))
end    
% 
% transfer the zoom level to other clusters
patch_zoom_level = level_zoom(classes);


for i=1:numComponents
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
    img_gray = rgb2gray(img);
%     BW = roicolor(img_gray,1,255);
%     figure
%     imshow(BW)
%     hold on
%     props = regionprops(BW, 'BoundingBox');
%     for k = 1 : length(props)
%      eva.OptimalK   thisBB = props(k).BoundingBox;
%         hold on;
%         rectangle('Position', thisBB,'EdgeColor','r','LineWidth',2);
%     end
end    

%% Checking the accuracy of the clustering by comparing with actual zoom
elem = find(actual_zoom == patch_zoom_level);
accuracy = nnz(elem)*100/(m(1)*m(2));



%% Hierarchical clustering with in each zoom level



%% To visualize the patches with zoom level
% count=0;
% vis_patch = zeros(size(X,1),size(X,2));
% for i=1:m(1)
%     for j=1:m(2)
%         count=count+1;
%         vis_patch(1+(i-1)*patch_size:(i)*patch_size,1+(j-1)*patch_size:(j)*patch_size)=patch_zoom_level(count);
%     end
% end    
% figure
% imagesc(vis_patch)
% title('Zoom level')
% colormap(jet(numComponents))
% colorbar
% 
% % Then we use pixels in each zoom level for hierarchical each pixel with in each zoom level 
% % 
% diff_zoom_levels = unique(vis_patch);






%%% Block wise
% 
% YP = cell(size(diff_zoom_levels,1),1);
% ZP=YP; ind_block=YP;
% for i=1:size(diff_zoom_levels,1)
%     ind_block{i}=find(patch_zoom_level==diff_zoom_levels(i)); %[row,col]
%     data_block = patch_wv_data(:,ind_block{i})';
%     YP{i} = pdist(data_block);
%     ZP{i} = linkage(YP{i});
%     figure
%     dendrogram(ZP{i})
%     title(['Hierarchy for Zoom level ' num2str(diff_zoom_levels(i))]);
% end












