function nrm = norm2(x,d)
% Computes the element wise l2-norm with respect to dimension d
if nargin < 2
    d = ndims(x);
end
nrm = sqrt(sum(abs(x).^2,d));