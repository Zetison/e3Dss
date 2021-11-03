function Z = bessel_s(n,z,type,nu_a,U_pol,u_k,v_k)
%Returns the n'th spherical bessel function of kind "type" evaluated at
%every element in z

if isa(z,'sym')
    PI = vpa('pi');
    tiny = realmin('double');
elseif isa(z,'mp')
    PI = mp('pi');
    tiny = realmin('double');
else
    tiny = eps;
    PI = pi;
end
if type == 1
    indices = logical(abs(z) < tiny);
    z(indices) = 1; % Avoid symbolic division by zero. These values will be replaced anyways
end

Z = sqrt(PI/2)./sqrt(z).*bessel_c(n+0.5,z,type,nu_a,U_pol,u_k,v_k);

if type == 1
    if n == 0
        Z(indices) = 1;
    else
        Z(indices) = 0;
    end
end
if ~(isa(z,'sym') || isa(z,'mp'))
    if any(isinf(Z(:)))
        warning('e3Dss:infBessel','A Bessel function evaluation was too large')
    end
end