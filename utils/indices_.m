function indices = indices_(nu,z,nu_a)
% Find indices at which the asymptotic expansion for large orders for the
% Besselfunctions is supposed to be used
if nu_a == -1
    if numel(nu) > numel(z)
        indices = false(size(nu));
    else
        indices = false(size(z));
    end
else
    indices = and(real(-exponent_(1,nu,z,nu_a)) < log(sqrt(realmin('double'))), nu > nu_a);
end