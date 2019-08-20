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
% s = 30;
% r = 30;
% k= 10;

% no.of patches
m = size(X)/patch_size; 
% disp(m)


% patch_data = zeros([patch_size*patch_size,m(1)*m(2)]);
patch_data = zeros([patch_size*patch_size*3,m(1)*m(2)]);

% wavelet pixel data

pixel_data = zeros(size(X,1),size(X,2),3*(1+3*max_level));

%% Reconstructing each patch using CUR
count = 0;
Img_reconst = zeros(size(X));
actual_zoom = zeros(m(1)*m(2),1);

for i=1:m(1)
    for j=1:m(2)
        patch = X(1+(i-1)*patch_size:(i)*patch_size,1+(j-1)*patch_size:(j)*patch_size,:);
                
%         for l=1:3
%             patch_cur(:,:,l) = CUR_v2(patch(:,:,l),k,s,r);
%         end    
%         patch_cur = CUR_v2(patch,k,s,r);
        patch_cur = patch;
        XX(1+(i-1)*patch_size:(i)*patch_size,1+(j-1)*patch_size:(j)*patch_size,:) = patch_cur;
        count = count + 1;
%         disp(count)
        patch_data(:,count) = reshape(patch_cur,[patch_size*patch_size*3,1]);
        
        %%% 2d wavelet
%         [C,S,actual_zoom(count)] = zoom_level(patch,basis);
%         [Ct,~] = wavelet_rgb(patch,basis);
%         patch_wv_data(:,count) = Ct;

        
        %%% 3d wavelet
        [Ct,actual_zoom(count)] = zoom_level_3d(patch,basis);
        patch_wv_data(:,count) = Ct;
        
        %%% wavelet_pixel
        for k=1:3
            pixel_data(1+(i-1)*patch_size:(i)*patch_size,1+(j-1)*patch_size:(j)*patch_size,1+(k-1)*(1+3*max_level):k*(1+3*max_level)) = wavelet_pixel(patch(:,:,k),basis);        
        end    
        
    end    
end

% figure
% imshow(uint8(XX))
% figure
% v = scree_plot_new(patch(:,:,1));

%% SVD
%[U, S, V] = svd(patch_data);

figure
% [U,Sig,V,v] = scree_plot_new(patch_data);

%% Wavelet coeff
[U,Sig,V,v] = scree_plot_new(patch_wv_data);
% 
rank = 50;
Coeff = Sig(1:rank,1:rank)*V(:,1:rank)';
% 
% patch_data_svd = U(:,1:rank)*Coeff;



%% Clustering
n = 10; % till how many you want to fit
klist= 2;
% klist=2:n;%the number of clusters you want to try
% myfunc = @(X,K)(kmeans(X, K));
% eva = evalclusters(net.IW{1},myfunc,'CalinskiHarabasz','klist',klist)
% classes=kmeans(net.IW{1},eva.OptimalK);

%rng('default');  % For reproducibility 'CalinskiHarabasz' 'gap','silhouette'
% num_trials = 10;
num_trials = 10;
temp = inf;

%% Finding how many components or clusters 
% Coefficient
data = Coeff';

% Original data
% data = patch_wv_data;



AIC = zeros(num_trials,n);
GMModels = cell(num_trials,n);
for i = 1:n
    for j=1:num_trials
        state = rng(j);
        GMModels{j}{i} = fitgmdist(data,i,'Regularize', 1e-5);
        AIC(j,i)= GMModels{j}{i}.AIC;
    end
end

minAIC = min(AIC(:));
[trial,numComponents] = find(AIC==minAIC);
BestModel = GMModels{trial}{numComponents};



%% Hierarchical clustering
% 
% numComponents = 6;
% 
% Y = pdist(data);
% Z = linkage(Y);
% dendrogram(Z)
% classes_hier=cluster(Z,'maxclust',numComponents);
% 
% %% Checking the accuracy of the clustering by comparing with actual zoom
% 
% level_zoom = zeros(numComponents,1);
% % transfer the zoom level to other clusters
% patch_zoom_level = level_zoom(classes_hier);
% elem = find(actual_zoom == patch_zoom_level);
% accuracy = nnz(elem)*100/(m(1)*m(2));




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
%     [C,S,level_zoom(i)] = zoom_level(patch_img,basis);
    
    %3D
    [Ct,level_zoom(i)]=zoom_level_3d(patch_img,basis);
    figure
    imshow(uint8(patch_img))
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
    BW = roicolor(img_gray,1,255);
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


%% Fitting gaussian in each cluster
GMM_component = cell(1,numComponents);
norm_weights = cell(1,numComponents);
for i =1:numComponents
    ind=find(classes==i);
    if size(ind,1)< rank
        GMM_component{i} = 0;
        norm_weights{i} = 0;
    else
        GMM_component{i} = fitgmdist(data(ind,:),1); 
        Var = GMM_component{i}.Sigma;
        mean_data = mean(data(ind,:));
        norm_weights{i} = vecnorm((data(ind,:)-mean_data)/Var,2,2);
    end    
end
% gm = GMM_component{i}; 
% X1 = data(ind,:);
% scatter(X1(:,1),X1(:,2),10,'.') % Scatter plot with points of size 10
% hold on
% gmPDF = @(x,y)reshape(pdf(gm,[x(:) y(:)]),size(x));
% fcontour(gmPDF,[-8 6])


%% Hierarchical clustering with in each zoom level



%% To visualize the patches with zoom level
count=0;
vis_patch = zeros(size(X,1),size(X,2));
for i=1:m(1)
    for j=1:m(2)
        count=count+1;
        vis_patch(1+(i-1)*patch_size:(i)*patch_size,1+(j-1)*patch_size:(j)*patch_size)=patch_zoom_level(count);
    end
end    
figure
imagesc(vis_patch)
title('Zoom level')
colormap(jet(numComponents))
colorbar

% Then we use pixels in each zoom level for hierarchical each pixel with in each zoom level 
% 
diff_zoom_levels = unique(vis_patch);



%%% Pixel wise

% YP = cell(size(diff_zoom_levels,1),1);
% ZP=YP; ind_pixel=YP;
% for i=1:size(diff_zoom_levels,1)
%     ind_pixel{i}=find(vis_patch==diff_zoom_levels(i)); %[row,col]
%     data_pixel = pixel_data(ind_pixel);
%     YP{i} = pdist(data_pixel);
%     ZP{i} = linkage(YP{i});
%     dendrogram(ZP{i})
%     title(['Hierarchy for Zoom level ' num2str(i)]);
% end

% pixel_data = reshape(pixel_data,size(pixel_data,1)*size(pixel_data,2),size(pixel_data,3));
% cgo = clustergram(pixel_data(1:1000,:),'Standardize','Row')



%%% Block wise

YP = cell(size(diff_zoom_levels,1),1);
ZP=YP; ind_block=YP;
for i=1:size(diff_zoom_levels,1)
    ind_block{i}=find(patch_zoom_level==diff_zoom_levels(i)); %[row,col]
    data_block = patch_wv_data(:,ind_block{i})';
    YP{i} = pdist(data_block);
    ZP{i} = linkage(YP{i});
    figure
    dendrogram(ZP{i})
    title(['Hierarchy for Zoom level ' num2str(diff_zoom_levels(i))]);
end









% %% Movie Test.
%  
% %% Set up some function. 
% % Sine between -2*pi and 2*pi.
% x = (10*-pi:0.1:10*pi)'; % Note the transpose.
% y = sin(x);
% fid = figure;
% hold on
% % The final plot.
% plot(x,y, '*');
%  
% %% Set up the movie.
% writerObj = VideoWriter('out.avi'); % Name it.
% writerObj.FrameRate = 60; % How many frames per second.
% open(writerObj); 
%  
% for i=1:size(y)      
%     % We just use pause but pretend you have some really complicated thing here...
%     pause(0.1);
%     figure(fId); % Makes sure you use your desired frame.
%     plot(x(i),y(i),'or');
%  
%     %if mod(i,4)==0, % Uncomment to take 1 out of every 4 frames.
%         frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
%         writeVideo(writerObj, frame);
%     %end
%  
% end
% hold off
% close(writerObj); % Saves the movie



