function PSNR = SLcomputePSNR(X,Xnoisy)
%SLcomputePSNR Compute peak signal to noise ratio (PSNR).
%
%Usage:
%
% PSNR = SLcomputePSNR(X, Xnoisy) 
%
%Input:
%
%      X: 2D or 3D signal.
% Xnoisy: 2D or 3D noisy signal.
%
%Output:
%
% PSNR: The peak signal to noise ratio (in dB).

    PSNR = 20*log10(255/(sqrt(sum(sum(sum(abs(X-Xnoisy).^2))))/sqrt(size(X,1)*size(X,2)*size(X,3))));
end

%
%  Copyright (c) 2014. Rafael Reisenhofer
%
%  Part of ShearLab3D v1.1
%  Built Mon, 10/11/2014
%  This is Copyrighted Material
%
%  If you use or mention this code in a publication please cite the website www.shearlab.org and the following paper:
%  G. Kutyniok, W.-Q. Lim, R. Reisenhofer
%  ShearLab 3D: Faithful Digital SHearlet Transforms Based on Compactly Supported Shearlets.
%  ACM Trans. Math. Software 42 (2016), Article No.: 5.
