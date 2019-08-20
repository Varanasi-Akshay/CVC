% clc;close all;clear;imtool close all;
% load woman
% addpath('/home/akshay/Desktop/CVC/Multi_target/Code/Multi-level_Sampling');
% level =2;
% wavelet = 'haar';
% [C,S]=wavedec2(X,level,wavelet);
% [dec]=plotwavelet2(C,S,level,wavelet,255,'square')


function [Coeff_pixel]=wavelet_pixel(X,wavelet)
% function is for single band


level= wmaxlev(size(X),wavelet);
[C,S]=wavedec2(X,level,wavelet);

%% tracking the pixel

% one approx and 3 details in each level.
D = cell(1,level); H = D; V = D; 
n = 1+3*level;
Coeff_pixel = zeros(size(X,1),size(X,2),n);
% temp_x = zeros(1,level); 
% temp_y = zeros(1,level);

for i=1:size(X,1)
    for j=1:size(X,2)
        temp_x = i;
        temp_y = j;

        for k=1:level
            x = round(temp_x/(2)); 
            y = round(temp_y/(2));
            m=rem(temp_x,2);
            n=rem(temp_y,2);
            if (m==0 && n==0) %d
                h=-1;
                v=-1;
                d = 1;
            elseif (m==0 && n==1)  %c
                h = 1;
                v = -1;
                d = -1;
            elseif (m==1 && n==0) % b
                h = -1;
                v = 1;
                d = -1;
            else % a
                h = 1;
                v = 1;
                d = 1;
            end    
                
                
            
            [H{k}, V{k}, D{k}] = detcoef2('a',C,S,k); % details
            l=3*(k-1)+1;
            Coeff_pixel(i,j,l)=H{k}(x,y)*h;
            Coeff_pixel(i,j,l+1)=V{k}(x,y)*v;
            Coeff_pixel(i,j,l+2)=D{k}(x,y)*d;
            
            if k==level
                l=3*k+1;
                A = appcoef2(C,S,wavelet,k);
                Coeff_pixel(i,j,l)=A(x,y);
            end
            
            temp_x = x;
            temp_y = y;
        end
    end   
end    

