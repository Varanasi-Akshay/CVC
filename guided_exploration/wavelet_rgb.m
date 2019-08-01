function [Ct,S] = wavelet_rgb(img,basis)
s = size(img(:,:,1));
level = wmaxlev(s,basis);
for i=1:3
   [C(:,i),S] = wavedec2(img(:,:,i),level,basis); 
end
Ct = reshape(C,[1,size(C,1)*size(C,2)]);
end