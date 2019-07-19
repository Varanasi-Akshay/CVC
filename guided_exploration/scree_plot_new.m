function [U,S,V,v] = scree_plot_new(A)

[U,S,V] = svd(A,'econ');
s = diag(S);

n = length(s);

v = zeros(n-1,1);

for r = 1:n-1
    v(r) = sum(s(r+1:end).^2)/s(1)^2;
end

%close all;
plot(1:n-1,v)
end
% ------------------------------------------------
% v = scree_plot_new(X);
% plot(1:100,v(1:100)) % close inspection