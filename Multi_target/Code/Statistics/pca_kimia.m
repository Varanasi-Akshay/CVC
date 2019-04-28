data = [];
D = '/home/akshay/Downloads/Multitarget/Codes/Statistics/KIMIAPath24RGB_zip/KIMIAPath24 RGB/Training/s';
for i=0:2
    folder = strcat(D,num2str(i));
    S = dir(fullfile(folder,'*.jpg')); % pattern to match filenames.

    for k = 1:numel(S)
        F = fullfile(folder,S(k).name);
        I = imread(F);
        I = imresize(I, [256, 256]);
        %imshow(I)
        I = double(I);
        X = reshape(I,size(I,1)*size(I,2)*size(I,3),1);
        data = [data X];
        %S(k).data = I; % optional, save data.
    end
end    



%% Patchwise 
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


[coeff,score,latent] = pca(data,'NumComponents',23);
Ipc=reshape(score(:,1),256,256,3); %% First principal component of the patch
imshow(Ipc)
save('kimia_data.mat')