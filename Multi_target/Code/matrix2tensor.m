function tensorblock = matrix2tensor(M1,block_size)

dim1 = size(M1,1);
dim2 = size(M1,2);
cut1 = floor(dim1/block_size)*block_size;
cut2 = floor(dim2/block_size)*block_size;
M = M1(1:cut1,1:cut2);


cut_param1 = block_size*ones(1,cut1/block_size);
cut_param2 = block_size*ones(1,cut2/block_size);

blocks=mat2cell(M,cut_param1, cut_param2);
num_blocks = size(blocks,1)*size(blocks,2);
tensorblock = zeros(block_size,block_size,num_blocks);
for i=1:num_blocks
    tensorblock(:,:,i) = blocks{i}(:,:);
end
