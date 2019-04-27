function [rank_sing,small_dim]= rankcalc(in_path, type, out_path, color_dim)
clc

% Set default parameters
if ~exist('in_path','var')
    in_path = 'stains_775x522';
end
if ~exist('type','var')
    type = '*';
end
if ~exist('out_path','var')
    out_path = in_path;
end
if ~exist('color_dim','var')
    color_dim = 1024;
end

% File i/o
in_path = strcat(in_path,'/');
files = dir(strcat(in_path,type));

for i=1:length(files)
    
    filename = files(i).name;
    img=imread(strcat(in_path,filename));
    
    % Convert multispectral image to single value per pixel
    img_ind = rgb2ind(img, color_dim);
    M1 = double(img_ind); 
    mw = size(M1,1);
    mh = size(M1,2);
    
    block_sizes = round(linspace(1,min(mh,mw)/4,200));
    
    % Pre-allocate
    rank_sing = zeros(1,size(block_sizes,2));
    small_dim = ones(1,size(block_sizes,2));
    
    for i=2:size(block_sizes,2)
        % Compute rank of blocked matrix
        [rank_sing(i),small_dim(i)] = rankfind(M1,block_sizes(i));
    end
    
    % Plot and save
    f1 = figure('visible','off');
    plot(block_sizes,rank_sing./small_dim)
    xlabel('block size')
    ylabel('Rank/Full Rank')
    
    saveas(f1, strcat(out_path,'/',filename(1:end-4),'_out.png'));
end

