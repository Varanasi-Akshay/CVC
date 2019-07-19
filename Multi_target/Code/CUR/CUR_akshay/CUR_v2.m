function XX = CUR_v2(X,k,s,r)
% CUR decomposition based on column/row lengths with various core matrix
% construction schemes
%
% ``CUR decompositions, approximations, and perturbations''

%% sampling based on column/row lengths
prob_col = sum(X.^2,1);
prob_row = sum(X.^2,2);
prob_col = prob_col/sum(prob_col);
prob_row = prob_row/sum(prob_row);
prob_col_cum = cumsum(prob_col);
prob_row_cum = cumsum(prob_row);

id_col = zeros(1,s);
id_row = zeros(r,1);

for j = 1:s
    p = rand(1);
    id = find(prob_col_cum>p);
    id_col(j) = id(1);
end
id_col = sort(id_col);

C = X(:,id_col);

for i = 1:r
    p = rand(1);
    id = find(prob_row_cum>p);
    id_row(i) = id(1);
end
id_row = sort(id_row);

R = X(id_row,:);

%% U matrix construction
method = 2; 
switch method
    case 1 % oracle
        U = pinv(C)*X*pinv(R);
    case 2 % intersection + low-rank approximation
        A = X(id_row,id_col);
        [Q,S,P] = svds(A,k);
        A = Q*S*P';
        U = pinv(A);
    %case 3 % constructed from C
    %    A = C'*C;
    %    [Q,S,P] = svds(A,k);
    %    A = Q*S*P';
    %    U = pinv(A);
    %    U = U*(C(id_row,:)');
end

%% approximation
XX = C*U*R;