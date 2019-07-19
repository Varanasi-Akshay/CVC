function [Ltrain,Ltest] = LKDL(Xtrain,Xtest,sigma)


n = size(Xtrain,2);
s = ceil(0.2*n); % 20% columns
k = ceil(0.8*s); % can be tuned

S = sort(randsample(n,s)); % uniform sampling
X_R = Xtrain(:,S);
gamma = 1/(2*sigma*sigma);
C_train = exp(-gamma .* pdist2(Xtrain',X_R','euclidean').^2);
W = C_train(S,:);
[UW, SW, VW] = svd(W);
SW = diag(SW); 
SW_inverse_root = 1./sqrt(SW(1:k));
Ltrain = diag(SW_inverse_root)*VW(:,1:k)'*C_train'; % K is approzimated by L*L'
C_test = exp(-gamma .* pdist2(Xtest',X_R','euclidean').^2);

Ltest = diag(SW_inverse_root)*VW(:,1:k)'*C_test';