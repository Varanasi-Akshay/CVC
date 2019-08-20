% start = [1 1];
% count = [5 3];
% data = h5read('example.h5','/g4/world',start,count)

%% Reading the data one row of chunks

clc;clear; close all;
addpath('/home/akshay/Desktop/Data/Biopsy');
filename = 'biopsy_converted.hdf5';
hinfo = h5info(filename); %,'V71Dimensions',true
datasetname = '/';

patch_size = 256;
wavelet = 'haar';
max_level = wmaxlev(patch_size,wavelet);
min_level = 2;

%%% we go till level 5 due to memory issue
% till which level should I consider, 

level_consider=5;
t=num_coeff(level_consider);


data_size = hinfo.Datasets(1).Dataspace.Size;
m = floor(data_size(1)/patch_size);
n = floor(data_size(2)/patch_size);

 
% m = 10; n = 10;

start = [1 1]; % starting point
chunk_size = [256 256]; % how many in each dimension (Chunk size)
%h5disp(filename)


%%% Number of patches
patch_data = zeros(m*n,3*t);
count=0;
patch_index=zeros(m*n,2);
patch_zoom = zeros(m*n,1);
T = zeros(patch_size,patch_size,3);


for i=1:m
   for j=1:n
       
       count=count+1;
       disp(count)
       start = [1+(i-1)*patch_size,1+(j-1)*patch_size ];
       %%% Red
       data = h5read(filename,strcat(datasetname,'r'),start,chunk_size);
       T(:,:,1)=data;
%        [C,~]= wavedec2(data,max_level,wavelet);
%        patch_data(count,1:t) = C(1:t);
       
       %%% Green
       data = h5read(filename,strcat(datasetname,'g'),start,chunk_size);
       T(:,:,2)=data;
%        [C,~]= wavedec2(data,max_level,wavelet);
%        patch_data(count,t+1:2*t) = C(1:t);
       
       %%% Blue
       data = h5read(filename,strcat(datasetname,'b'),start,chunk_size);
       T(:,:,3)=data;
%        [C,S]= wavedec2(data,max_level,wavelet);
%        patch_data(count,2*t+1:3*t) = C(1:t);
       
       [~,patch_zoom(count)] = zoom_level_3d(T,wavelet);
       patch_index(count,:) = [i,j];
   end
end 

save('biopsy_zoom.mat')
% 
% %% Building the tree
% 
% %%% Initial clustering is based on 
% start_level = 2;
% 
% %%% final clustering is based on
% end_level = 5;
% numComponents = 2;
% % cluster_level = cell(end_level-start_level+1);
% 
% i=start_level;
% while i<=end_level
%     tl=num_coeff(i);
%     ind = [1:tl,1*t+1:1*t+tl,2*t+1:2*t+tl];
%     data= patch_data(:,ind);
%     [classes_temp, centers_temp, sumd_temp, D_temp ]=kmeans(data,numComponents);
%     cluster_level{i-1} = classes_temp;
%     i=i+1;
%     for j=1:numComponents
%         tl=num_coeff(i);
%         ind = [1:tl,1*t+1:1*t+tl,2*t+1:2*t+tl];
%         ind_patch=find(classes_temp==j);
%         if (size(ind_patch,1) < size(ind,1))
%             continue
%         end    
%         data = patch_data(ind_patch,ind);
%         [classes_temp, centers_temp, sumd_temp, D_temp ]=kmeans(data,numComponents);
%         cluster_level{i-1,j} = classes_temp;
%     end    
%     break
% end

function t=num_coeff(level_consider)
t=1;
for i=1:level_consider
    t = t + 3*(2^(2*(i-1)));
end    
end

