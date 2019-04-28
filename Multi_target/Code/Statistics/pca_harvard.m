% load('img1.mat')
% I=ref;
% [a,b,c] = size(I);
% patch_size = 8;
% p=(a/patch_size);
% q=(b/patch_size);
% data = [];
% for i=1:p
%     for j=1:q
%        patch = I(1+(i-1)*patch_size:(i)*patch_size,1+(j-1)*patch_size:(j)*patch_size,:); 
%        X = reshape(patch,size(patch,1)*size(patch,2)*size(patch,3),1);
%        data = [data X];
%     end
% end   
load('img1_patch_data.mat')
patch_size = 8;
[coeff,score,latent] = pca(data,'NumComponents',100);
Ipc=reshape(score(:,1),patch_size,patch_size,31); %% First principal component of the patch
save('img1_patch_data.mat')