% Set number of iterations and wavelet name. 
iter = 10;
wav = 'haar';

waveinfo('haar')

% Compute approximations of the wavelet function using the
% cascade algorithm. 
for i = iter 
    [phi,psi,xval] = wavefun(wav,i); 
    plot(xval,psi); 
    hold on 
end
title(['Approximations of the wavelet ',wav, ... 
       ' for 1 to ',num2str(iter),' iterations']); 
hold off