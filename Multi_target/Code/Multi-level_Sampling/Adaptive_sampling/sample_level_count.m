function [sample_level]=sample_level_count(I,IND,level)

% showing the sampling pattern
sample=zeros(size(I));
sample(IND)=1;
%imshow(sample)

% Finding no.of samples at each level
max_level=log2(size(I,1));
sample_level=zeros(level,1);
for i=1:level
    count=nnz(sample(1:2^(max_level-level+i),1:2^(max_level-level+i)));
    sample_count = count-sum(sample_level);
    sample_level(i) = sample_count;
end