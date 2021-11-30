function indices = indices_(nu,z,nu_a)
% Find indices at which the asymptotic expansion for large orders for the
% Besselfunctions is supposed to be used
if nu_a == -1
    indices = false(size(z));
else
    indices = and(real(-exponent_(1,nu,z,nu_a)) < log(sqrt(realmin(class(z)))), nu > nu_a);
end