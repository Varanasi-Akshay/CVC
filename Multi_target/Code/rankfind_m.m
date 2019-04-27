function [rank_sing,small_dim]= rankfind(M1,block_size)

dim1 = size(M1,1);
dim2 = size(M1,2);
cut1 = floor(dim1/block_size)*block_size;
cut2 = floor(dim2/block_size)*block_size;
M = M1(1:cut1,1:cut2);

cut_param1 = block_size*ones(1,cut1/block_size);
cut_param2 = block_size*ones(1,cut2/block_size);

blocks=mat2cell(M,cut_param1, cut_param2);
num_blocks = size(blocks,1)*size(blocks,2);
blocks_reshaped = zeros(num_blocks,block_size*block_size);
for i=1:num_blocks
    column = blocks{i}(:);
    blocks_reshaped(i,:) = column;
    
end    
[~, S, ~] = svd(blocks_reshaped,'econ');% columnbasis = U(:,logical(S))
t = 0.01;
singular_ratio = zeros(size(S,1)-1,1);
rank_sing = 0;
for i=1:size(S,1)-1
    singular_ratio(i)=S(i,i)/S(1,1);
    if singular_ratio(i)> t
        rank_sing = i;
    end    
end

small_dim=min(size(blocks_reshaped,1), size(blocks_reshaped,2));