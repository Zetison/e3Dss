function nrm = norm2(x,conjugation)
% x is a MxN matrix. Computes the norms nrm(i) = norm(x(i,:))

M = size(x,1);

nrm = zeros(M,1);
if isa(x,'sym')
    nrm = vpa(nrm);
end
if nargin > 1
    nrm = sqrt(sum(x.*conj(x),2));
else
    nrm = sqrt(sum(x.*x,2));
end