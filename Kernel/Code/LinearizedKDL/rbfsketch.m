function [C] = rbfsketch(X,sigma,s)
n = size(X,1);
idx = sort(randsample(n,s));
C = rbf(X,X(idx,:),sigma);