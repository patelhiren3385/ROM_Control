function [xi,vv]=col(n)
% Gauus- Hermite nodes and weights

ii   = 1:n-1;
a   = sqrt(ii/2);
CM  = (sqrt(2))*(diag(a,1) + diag(a,-1));

% Determining the abscissas (x) and weights (w)
    % - since det(xI-CM)=L_n(x), the abscissas are the roots of the
    %   characteristic polynomial, i.e. the eigenvalues of CM;
    % - the weights can be derived from the corresponding eigenvectors.
[V L]   = eig(CM);
[xi ind] = sort(diag(L));
V       = V(:,ind)';

vv= V(:,1).^2;