%% give me a RGB named 'img' in double format
% make sure tensor toolbox is working

addpath(genpath('~/tensor_toolbox-master/'))
img=imread('Semal_(Bombax_ceiba)_flowers_in_Kolkata_W_IMG_4132.jpg');
img=double(img);
for i = 1:3
    
    img = double(img);
    
    rank = 20*i;
    
    R = img(:,:,1);
    
    [U,S,V] = svd(R,'econ');
    U = U(:,1:rank);
    S = S(1:rank,1:rank);
    V = V(:,1:rank);
    R = U*S*V.';
    
    G = img(:,:,2);
    
    [U,S,V] = svd(G,'econ');
    U = U(:,1:rank);
    S = S(1:rank,1:rank);
    V = V(:,1:rank);
    G = U*S*V.';
    
    B = img(:,:,3);
    
    [U,S,V] = svd(B,'econ');
    U = U(:,1:rank);
    S = S(1:rank,1:rank);
    V = V(:,1:rank);
    B = U*S*V.';
    
    I1 = zeros(size(img));
    I1(:,:,1) = R;
    I1(:,:,2) = G;
    I1(:,:,3) = B;
    
    % cp
    M = cp_als(tensor(img),rank,'maxiters',999);
    A = M.U{1};
    B = M.U{2};
    C = M.U{3};
    v = M.lambda;
    I2 = ktensor(v,{A,B,C},[1 2 3]);
    I2 = double(I2);
    
    % tucker
    M = tucker_als(tensor(img),[rank rank 3],'maxiters',999);
    G = M.core;
    A = M.U{1};
    B = M.U{2};
    C = M.U{3};    
    I3 = ttm(G,{A,B,C},[1 2 3]);
    I3 = double(I3);
    
    Collection{i,1} = I1;
    Collection{i,2} = I2;
    Collection{i,3} = I3;
    
end


%I1 = Collection{5,1};
%I2 = Collection{5,2};
%I3 = Collection{5,3};
%imshow([I1;I2;I3])
