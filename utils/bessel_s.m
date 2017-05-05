function Z = bessel_s(n,z,i)
%Returns the n'th spherical bessel function of kind "type" evaluated at
%every element in z

if isa(z,'sym')
    Z = sqrt(vpa(pi)/2)./sqrt(z).*bessel_c(n+0.5,z,i);
else
    Z = sqrt(pi/2)./sqrt(z).*bessel_c(n+0.5,z,i);
end

if i == 1
    if n == 0
        Z(logical(abs(z) < eps)) = 1;
    else
        Z(logical(abs(z) < eps)) = 0;
    end
end
if ~isa(z,'sym')
    if any(any(abs(Z) > 1e290))
        error('A Bessel function evaluation was too large')
    end
end