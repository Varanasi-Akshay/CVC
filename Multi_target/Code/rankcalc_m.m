clc
clear all

path = 'stains_775x522/';
color_dim = 1024*4;
files = dir(strcat(path,'*.bmp'));
for i=1:length(files)
    filename = files(i).name;
    img=imread(strcat(path,filename));
    img_ind = rgb2ind(img, color_dim);
    M1 = double(img_ind); 
    mw = size(M1,1);
    mh = size(M1,2);
    
    block_size = round(linspace(1,min(mh,mw)/2,100));
    
    rank_sing = zeros(1,size(block_size,2));
    small_dim = ones(1,size(block_size,2));
    for i=2:size(block_size,2)-1
        [rank_sing(i),small_dim(i)] = rankfind(M1,block_size(i));
    end
    f1 = figure('visible','off');
    plot(block_size,rank_sing./small_dim)
    xlabel('block size')
    ylabel('Rank/Full Rank')
    
    saveas(f1, strcat(filename(1:end-4),'.png'));

    %basis=orth(blocks_)
end

