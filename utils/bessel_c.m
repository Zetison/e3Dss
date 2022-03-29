function f = bessel_c(nu,z,i,nu_a,U_pol,u_k,v_k,Eps)

indices = indices_(nu,z,nu_a);
f = zeros(size(z),class(z));
if numel(nu) > 1
    f(~indices) = bessel_std(nu(~indices),z,i,nu_a);
    if any(indices(:))
        f(indices) = bessel_asy(nu(indices),z,i,U_pol,u_k,v_k,Eps);
    end
else
    f(~indices) = bessel_std(nu,z(~indices),i,nu_a);
    if any(indices(:))
        f(indices) = bessel_asy(nu,z(indices),i,U_pol,u_k,v_k,Eps);
    end
end
end

function f = bessel_std(nu,z,i,nu_a)
if i == 1 % besselj
    f = besselj(nu,z,1).*exp(abs(imag(z)) + exponent_(i,nu,z,nu_a));
elseif i == 2 % bessely
    f = bessely(nu,z,1).*exp(abs(imag(z)) + exponent_(i,nu,z,nu_a));
else
    f = besselh(nu,1,z,1).*exp(1i*z + exponent_(i,nu,z,nu_a)); % Hankel function of first kind
end
end

function f = bessel_asy(nu,z,i,U_pol,u_k,v_k,Eps)
% Reference: https://dlmf.nist.gov/10.20
prec = class(z);
PI = getC(prec,'pi');
oneThirds = getC(prec,'1/3');
twoThirds = getC(prec,'2/3');
fiveThirds = getC(prec,'5/3');
K = numel(U_pol)/2-1;
sum1 = zeros(size(z),class(z));
sum2 = zeros(size(z),class(z));
useScaling = 1;
% useScaling = 0;
% sum1_k_arr = zeros(K+1,max(length(z),length(nu)),class(z));
% sum2_k_arr = zeros(K+1,max(length(z),length(nu)),class(z));
zdnu = z./nu;
zeta23 = zeta23_(zdnu);
zeta = zeta_(zdnu,zeta23);
arg1mz2 = (1-zdnu.^2).^(-1/2);
for k = 0:K
    nu2k = nu.^(2*k);
    sum1_k = A_k(k,U_pol,v_k,     zeta23,arg1mz2)./nu2k;
    sum2_k = B_k(k,U_pol,u_k,zeta,zeta23,arg1mz2)./nu2k;
%     sum1_k_arr(k+1,:) = sum1_k;
%     sum2_k_arr(k+1,:) = sum2_k;
    sum1 = sum1 + sum1_k;
    sum2 = sum2 + sum2_k;
    if all(abs(sum1_k)./abs(sum1) < Eps) && all(abs(sum2_k)./abs(sum2) < Eps) 
        break
    end
end
% semilogy(0:K,abs(sum1_k_arr),0:K,abs(sum2_k_arr))
if k == K
    error('e3Dss:divergeBessel','Series did not converge')
end
nu23zeta = nu.^twoThirds.*zeta;
if i == 3
    nu23zeta = nu23zeta*exp(2*PI*1i/3);
    airytype = 0;
else
    airytype = 2*(i-1);
end
sum1 = airy(airytype,  nu23zeta,useScaling)./nu.^oneThirds.*sum1;
sum2 = airy(airytype+1,nu23zeta,useScaling)./nu.^fiveThirds.*sum2;
if i == 3
    sum2 = sum2*exp(2*PI*1i/3);
end

f = (-1)^(i+1)*(4*zeta./(1-zdnu.^2)).^(1/4).*(sum1+sum2);
if i == 3
    f = f.*2*exp(-PI*1i/3);
end
end

function A = A_k(k,U_pol,v_k,zeta23,arg1mz2)

A = zeros(size(zeta23),class(zeta23));
for j = 0:2*k
    A = A + v_k(j+1)*zeta23.^(-j).*U_k(2*k-j,arg1mz2,U_pol);
end
end

function B = B_k(k,U_pol,u_k,zeta,zeta23,arg1mz2)

B = zeros(size(zeta23),class(zeta23));
for j = 0:(2*k+1)
    B = B + u_k(j+1)*zeta23.^(-j).*U_k(2*k-j+1,arg1mz2,U_pol);
end
B = -1i*(-zeta).^(-1/2).*B;
end

function U = U_k(k,z,U_pol)
U = polyval(U_pol{k+1},z);
end
