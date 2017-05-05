function d2Z = d2bessel_s(n,z,i,Z)
% Returns the second derivative of the n'th spherical bessel function of kind i
% evaluated at every element in z

d2Z = (n*(n-1)./z.^2-1).*Z{i,1} + 2./z.*Z{i,2};
if i == 1
    if n == 0
        d2Z(logical(abs(z) < eps)) = -1/3;
    elseif n == 2
        d2Z(logical(abs(z) < eps)) = 2/15;
    else
        d2Z(logical(abs(z) < eps)) = 0;
    end
end
if ~isa(z,'sym')
    if any(any(abs(d2Z) > 1e290))
        error('A Bessel function evaluation was too large')
    end
end