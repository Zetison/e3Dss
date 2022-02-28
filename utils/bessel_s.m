function Z = bessel_s(n,z,i,nu_a,U_pol,u_k,v_k,Eps)
%Returns the n'th spherical bessel function of kind "i" evaluated at
%every element in z

if isa(z,'sym')
    tiny = realmin('double');
elseif isa(z,'mp')
    tiny = realmin('double');
else
    tiny = eps;
end
PI = getC(class(z),'pi');
if i == 1
    indices = logical(abs(z) < tiny);
    z(indices) = ones(1,class(z)); % Avoid symbolic division by zero. These values will be replaced anyways
end

Z = sqrt(PI/2)./sqrt(z).*bessel_c(n+0.5,z,i,nu_a,U_pol,u_k,v_k,Eps);

if i == 1
    if n == 0
        Z(indices) = ones(1,class(z));
    else
        Z(indices) = zeros(1,class(z));
    end
end
if ~(isa(z,'sym') || isa(z,'mp'))
    if any(isinf(Z(:)))
        warning('e3Dss:infBessel','A Bessel function evaluation was too large')
    end
end