%% wavelet as matrix product

clc;clear;close all;
I = eye(3);
w = zeros(4,3);
for i=1:3
    w(:,i)=wavedec(I(:,i),2,'haar');
end    

%% checking for random matrix
A3 = rand(3,6);
wA = zeros(4,6);
wA_dec = zeros(4,6);
for i=1:6
    wA(:,i)=w*A3(:,i);
    wA_dec(:,i)= wavedec(A3(:,i),2,'haar');
end    

%% Works!!


%% 2D matrix
A = eye(4);
B=zeros(size(A)); C=B;
for i=1:4
   B(i,:)=wavedec(A(i,:),1,'haar');   
end

for i=1:4
   C(:,i)=wavedec(B(:,i),1,'haar');   
end

A1 = rand(4,4);
B1=zeros(size(A1)); C1=B1;
for i=1:4
   B1(i,:)=wavedec(A1(i,:),1,'haar');   
end

for i=1:4
   C1(:,i)=wavedec(B1(:,i),1,'haar');   
end

