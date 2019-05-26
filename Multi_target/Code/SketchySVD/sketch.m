%% SVD using the Sketchy SVD
function [Ar]=sketch(A,k,s,r)

[m,n] = size(A);

% Random Projections or test matrices
Gamma = rand(k,m)/sqrt(k);
Omega = rand(k,n)/sqrt(k);
Phi = rand(s,m)/sqrt(m);
Psi = rand(s,n)/sqrt(n);

X = Gamma*A;
Y = A*Omega';
Z = Phi*A*Psi';
[P,R1] = qr(X');
[Q,R2] = qr(Y);
C = pinv(Phi*Q)*Z*(pinv(Psi*P))';
[u,s,v] = svd(C);
Cr = u(:,1:r)*s(1:r,1:r)*v(:,1:r)';
Ar=Q*Cr*P';
end
