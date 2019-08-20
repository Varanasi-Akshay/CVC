function [T2,level_zoom] = zoom_level_3d(T,basis)
% TS = sparse_approx_3D(T,per)
% per: percentage of 3D coefficients to keep
details = ['h','v','d'];
threshold = 0.99;
level_zoom = 0; 
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

I = eye(3);
w = zeros(4,3);
for i=1:3
    w(:,i)=wavedec(I(:,i),l,basis);
end

for j = 1:size(T1,2)
   
%     [T2(:,j),L] = wavedec(T1(:,j),l,basis); % can be replaced by matrix multiplication
    T2(:,j)=w*T1(:,j);
   
end

T2 = T2.';
T2 = T2(:);

T2_original= T2;
L1 = zeros(size(T1,2),size(T1,1));

% L2 = zeros(size(T1,2),size(T1,1));

for i = 1:level
    
    for j = 1:n3
   
        I = T(:,j);
        I = reshape(I,n1,n2);
        thres_levels = 1:i;
        [C,S] = wavedec2(I,level,basis);
        NC = wthcoef2(details(1),C,S,thres_levels);
        NC = wthcoef2(details(2),NC,S,thres_levels);
        NC = wthcoef2(details(3),NC,S,thres_levels);

   
        L1(:,j) = NC(:);
   
    end

    LL1 = L1.';

    for j = 1:size(LL1,2)
   
%         [L2(:,j),L] = wavedec(LL1(:,j),l,basis); % can be replaced by matrix multiplication
        L2(:,j)=w*LL1(:,j);
   
    end

    LL2 = L2.';
    L3 = LL2(:);
    

    value = sumsqr(L3)/sumsqr(T2_original);

  
    if value < threshold
        level_zoom = level-i;
        break
    end   
end
