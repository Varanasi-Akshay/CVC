clc
clear all

%path = 'stains_775x522';
color_dim = 1024;
%files = dir (strcat(path,'\*.bmp'));
%for ioh=1:length(files)
    img=imread('10-12813-01-2.bmp');
    img_ind = rgb2ind(img, color_dim);
    M1 = double(img_ind); 
    mw = size(M1,1);
    mh = size(M1,2);

    % at 1000 this takes a bit to run
    block_size = round(linspace(1,min(mh,mw),1000));
    
    rank_sing = zeros(1,size(block_size,2));
    small_dim = ones(1,size(block_size,2));
    for i=2:size(block_size,2)-1
        [rank_sing(i),small_dim(i)] = rankfind(M1,block_size(i));
    end
    plot(block_size,rank_sing./small_dim)
    xlabel('block size')
    ylabel('Rank/Full Rank')

    %basis=orth(blocks_)
%end

