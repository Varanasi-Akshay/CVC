function y = partial_Fourier(x,S,IND,n1,n2,level,mode)

if mode == 1
    x = reshape(x,1,[]);
    z = waverec2(x,S,'db2');
    y = fftshift(fft2(z))/sqrt(n1*n2);
    y = y(IND);
else
    x = reshape(x,[],1);
    z = zeros(n1,n2);
    z(IND) = x;
    z = ifft2(ifftshift(z))*sqrt(n1*n2);
    [y,~] = wavedec2(z,level,'db2');
    y = reshape(y,[],1);
end