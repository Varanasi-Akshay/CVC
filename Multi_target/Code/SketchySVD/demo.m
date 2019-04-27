close all;
clear all;
%image = imread('10-12813-01-2.bmp');
image = imread('09-1339-01-2.bmp');
img = rgb2gray(image);
img = double(img)/255;
figure(1)
imshow(img)
k=150;
s=300;
r=50;
[imgr]=sketch(img,k,s,r);

figure(2)
imshow(imgr)
err1 = norm(img-imgr,'inf')/norm(img,'inf')*100;
fprintf('Inf norm error sketchysvd %d %% \n', err1)

err1 = norm(img-imgr,2)/norm(img,2)*100;
fprintf('2 norm error sketchysvd %d %% \n', err1)

err1 = norm(img-imgr,'fro')/norm(img,'fro')*100;
fprintf('fro norm error sketchysvd %d %% \n', err1)


[u,sing,v]=svd(img);
imgsvd = u(:,1:r)*sing(1:r,1:r)*v(:,1:r)';

figure(3)
imshow(imgsvd)

err2 = norm(img-imgsvd,'inf')/norm(img,'inf')*100;
fprintf('Inf norm error svd %d %% \n', err2)


err2 = norm(img-imgsvd,2)/norm(img,2)*100;
fprintf('2 norm error svd %d %% \n', err2)


err2 = norm(img-imgsvd,'fro')/norm(img,'fro')*100;
fprintf('fro norm error svd %d %% \n', err2)
