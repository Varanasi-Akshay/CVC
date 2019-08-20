function TS2 = sparse_approx_3D(T,basis,per)
[n1,n2,n3] = size(T);

N = n1*n2;
T = reshape(T,N,n3);

level = floor(log2(min(n1,n2)));
l = floor(log2(n3));

for j = 1:n3
   
    I = T(:,j);
    I = reshape(I,n1,n2);
   
    [C,S] = wavedec2(I,level,basis);
   
    T1(:,j) = C(:);
   
end

T1 = T1.';

for j = 1:size(T1,2)
   
    [T2(:,j),L] = wavedec(T1(:,j),l,basis);
   
end

T2 = T2.';
T2 = T2(:);



% m = floor(per*n1*n2*n3);

m = floor(per*size(T2));

[~,ind] = sort(abs(T2),'descend');
TS = zeros(size(T2));
TS(ind(1:m)) = T2(ind(1:m));

TS = reshape(TS,size(T1,2),[]);

for j = 1:size(TS,2)
   
    C = TS(:,j);

    C = reshape(C,1,[]);
    I = waverec2(C,S,basis);
   
    TS1(:,j) = I(:);
   
end

TS1 = TS1.';

for j = 1:N
   
    TS2(:,j) = waverec(TS1(:,j),L,basis);
   
end

TS2 = TS2.';
TS2 = reshape(TS2,n1,n2,n3);