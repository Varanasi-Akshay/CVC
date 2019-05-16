function [AA,AAA] = SketchySVD(A,r)
% implementation of SketchySVD in [Streaming low-rank matrix approximation 
% with an application to scientific simulation]

[m,n] = size(A);

alpha = 1;

% parameters (5.3)
k = 4*r+alpha; 
s = 2*k+alpha; 

type = 2; % 1 - Gaussian; 2 - SSRFT; 

switch type 
    case 1
        % left sketch
        T1 = randn(k,m);
        X = T1*A; % k by n
        [P,~] = qr(X',0); % n by k
        % right sktech
        T2 = randn(k,n);
        Y = A*T2'; % m by k
        [Q,~] = qr(Y,0); % m by k
        % core sketch
        T3 = randn(s,m);
        T4 = randn(s,n);
        Z = T3*A*T4'; % s by s
        L = pinv(T3*Q); % k by s
        R = pinv(T4*P); % k by s
    case 2
        % left sketch
        d = rand(m,1);
        d = d > 0.5;
        d = 2*d-1;
        TEMP1 = bsxfun(@times,d,A);
        clear d
        TEMP2 = dct(TEMP1);
        clear TEMP1
        ind = randsample(m,k);
        ind = sort(ind);
        X = TEMP2(ind,:);
        [P,~] = qr(X',0);
        clear TEMP2
        %disp('left sktech')
        % right sketch
        d = rand(n,1);
        d = d > 0.5;
        d = 2*d-1;
        TEMP1 = bsxfun(@times,d,A');
        clear d
        TEMP2 = dct(TEMP1);
        clear TEMP1
        ind = randsample(n,k);
        ind = sort(ind);
        Y = (TEMP2(ind,:))';
        [Q,~] = qr(Y,0);
        clear TEMP2
        %disp('right sketch')
        % core sketch
        d = rand(m,1);
        d = d > 0.5;
        d = 2*d-1;
        TEMP1 = bsxfun(@times,d,A);
        TEMP3 = bsxfun(@times,d,Q);
        clear d
        TEMP2 = dct(TEMP1);
        TEMP4 = dct(TEMP3);
        clear TEMP1 TEMP3
        ind = randsample(m,s);
        ind = sort(ind);
        TEMP = TEMP2(ind,:);
        L = TEMP4(ind,:);
        L = pinv(L);
        clear TEMP2 TEMP4
        d = rand(n,1);
        d = d > 0.5;
        d = 2*d-1;
        TEMP1 = bsxfun(@times,d,TEMP');
        TEMP3 = bsxfun(@times,d,P);
        clear d TEMP
        TEMP2 = dct(TEMP1);
        TEMP4 = dct(TEMP3);
        clear TEMP1 TEMP3
        ind = randsample(n,s);
        ind = sort(ind);
        Z = (TEMP2(ind,:))';
        R = TEMP4(ind,:);
        R = pinv(R);
        clear TEMP2 TEMP4
        %disp('core sketch')
    otherwise
        error('Test Matrices of Type Not Supported!')
end

% initial approx
C = L*Z*R'; % k by k
AA = Q*C*P';
%disp('initial approx')

% truncated approx
[U,S,V] = svds(C,r);
Q = Q*U(:,1:r); % m by r
P = P*V(:,1:r); % n by r
s = diag(S);
AAA = Q*diag(s)*P';