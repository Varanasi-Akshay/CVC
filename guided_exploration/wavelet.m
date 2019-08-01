clc,clear all; close all

X  = imread('26.jpg');
X = double(X);
figure
imshow(uint8(X))
% colormap(map)
% title('Original')
wv = 'haar';
wt = wavedec3(X,3,wv);
[img] = waverec3(wt);
figure
imshow(uint8(img))
n=3;
img_rec = zeros(size(X));
for i=1:3
   [C,S] = wavedec2(X(:,:,i),n,wv); 
   img_rec(:,:,i) = waverec2(C,S,wv);    
end    
figure
imshow(uint8(img_rec))

[Ct,S] = wavelet_rgb(X,wv);
