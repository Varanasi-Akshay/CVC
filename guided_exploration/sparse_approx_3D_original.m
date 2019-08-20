function TS = sparse_approx_3D_original(T,per)
% per: percentage of 3D coefficients to keep

[n1,n2,n3] = size(T);

N = n1*n2;
T = reshape(T,N,n3);

level = floor(log2(min(n1,n2)));
l = floor(log2(n3));

for j = 1:n3
   
    I = T(:,j);
    I = reshape(I,n1,n2);
   
    [C,S] = wavedec2(I,level,'haar');
   
    T(:,j) = C(:);
   
end

T = T.';

for j = 1:N
   
    [T(:,j),L] = wavedec(T(:,j),l,'haar');
   
end

T = T.';
T = T(:);

m = floor(per*n1*n2*n3);

[~,ind] = sort(abs(T),'descend');
TS = zeros(size(T));
TS(ind(1:m)) = T(ind(1:m));

TS = reshape(TS,N,[]);

for j = 1:size(TS,2)
   
    C = TS(:,j);

    C = reshape(C,1,[]);
    I = waverec2(C,S,'haar');
   
    TS(:,j) = I(:);
   
end

TS = TS.';

for j = 1:N
   
    TS(:,j) = waverec(TS(:,j),L,'haar');
   
end

TS = TS.';
TS = reshape(TS,n1,n2,n3);