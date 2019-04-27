function [I,S] = sparse_approx_shearlet(I,scale,per)

I = double(I);
[n1,n2] = size(I);

useGPU = 0;

sigma = 30;
thresholdingFactors = [0 3 3 4 4];

% create shearlet system 
system = SLgetShearletSystem2D(useGPU,n1,n2,scale);

% co efficients
coeff = SLsheardec2D(data,system);

m = floor(per*n1*n2);

% [~,ind] = sort(abs(C),'descend');
% CC = C;
% C = zeros(size(CC));
% C(ind(1:m)) = CC(ind(1:m));



% thresholding
for j =1:system.nShearlets
   shearletIdx = system.shearletIdxs(j,:);
   disp(shearletIdx);
   coeff(:,:,j) = coeff(:,:,j).*(abs(coeff(:,:,j)) > ...
   thresholdingFactors(shearletIdx(2)+1)*system.RMS(j)*sigma);
end
% Reconstruction
I = SLshearrec2D(coeff,system);

