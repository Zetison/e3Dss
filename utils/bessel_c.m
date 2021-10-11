function f = bessel_c(nu,z,type)

if type == 1 % besselj
    f = besselj(nu,z);
else % bessely
    f = bessely(nu,z);
end
% 
% indices = nu > z;
% f = zeros(size(z));
% if type == 1 % besselj
%     f(~indices) = besselj(nu,z(~indices));
%     if any(indices)
%         f(indices) = bessel_J_asy(nu,z(indices));
%     end
% else % bessely
%     f(~indices) = bessely(nu,z(~indices));
%     if any(indices)
%         f(indices) = bessel_Y_asy(nu,z(indices));
%     end
% end

end
function J = bessel_J_asy(nu,z)
load('../miscellaneous/U_pol.mat')
Eps = eps;
nu = nu + 1/2;
K = numel(U_pol)/2-1;
sum1 = 0;
sum2 = 0;
sum1_k_arr = zeros(K+1,1);
sum2_k_arr = zeros(K+1,1);
for k = 0:K
    sum1_k = A_k(k,z./nu,U_pol)./nu.^(2*k);
    sum2_k = B_k(k,z./nu,U_pol)./nu.^(2*k);
    sum1_k_arr(k+1) = sum1_k;
    sum2_k_arr(k+1) = sum2_k;
    sum1 = sum1 + sum1_k;
    sum2 = sum2 + sum2_k;
    if all(abs(sum1_k)./abs(sum1) < Eps) && all(abs(sum2_k)./abs(sum2) < Eps) 
        break
    end
end
semilogy(0:K,abs(sum1_k_arr),0:K,abs(sum2_k_arr))
if k == K
    error('Series did not converge')
end
sum1 = airy(0,nu.^(2/3).*zeta(z./nu),1)./nu.^(1/3).*sum1;
sum2 = airy(1,nu.^(2/3).*zeta(z./nu),1)./nu.^(5/3).*sum2;

J = (4*zeta(z./nu)./(1-(z./nu).^2)).^(1/4).*(sum1+sum2);
end

function Y = bessel_Y_asy(nu,z)
load('../miscellaneous/U_pol.mat')
Eps = eps;
nu = nu + 1/2;
K = numel(U_pol)/2-1;
sum1 = 0;
sum2 = 0;
for k = 0:K
    sum1_k = A_k(k,zeta(z./nu),U_pol)./nu.^(2*k);
    sum2_k = B_k(k,zeta(z./nu),U_pol)./nu.^(2*k);
    sum1 = sum1 + sum1_k;
    sum2 = sum2 + sum2_k;
    if all(abs(sum1_k)./abs(sum1) < Eps) && all(abs(sum2_k)./abs(sum2) < Eps) 
        break
    end
end
if k == K
    error('Series did not converge')
end
sum1 = airy(2,nu.^(2/3).*zeta(z./nu),1)./nu.^(1/3).*sum1;
sum2 = airy(3,nu.^(2/3).*zeta(z./nu),1)./nu.^(5/3).*sum2;

Y = -(4*zeta(z./nu)./(1-(z./nu).^2)).^(1/4).*(sum1+sum2);
end

function A = A_k(k,z,U_pol)

A = 0;
for j = 0:2*k
    A = A + (3/2)^j*v_k(j)*zeta(z).^(-3*j/2).*U_k(2*k-j,(1-z.^2).^(-1/2),U_pol);
end
end

function B = B_k(k,z,U_pol)

B = 0;
for j = 0:(2*k+1)
    B = B + (3/2)^j*u_k(j)*zeta(z).^(-3*j/2).*U_k(2*k-j+1,(1-z.^2).^(-1/2),U_pol);
end
B = -1i*(-zeta(z)).^(-1/2).*B;
end

function U = U_k(k,z,U_pol)
U = polyval(U_pol{k+1},z);
end

function u = u_k(k)
u = prod((2*k+1):2:(6*k-1))/216^k/factorial(k);
end

function v = v_k(k)
v = (6*k+1)/(1-6*k)*u_k(k);
end

function zeta = zeta(z)
zeta = (3/2*(log((1+sqrt(1-z.^2))./z)-sqrt(1-z.^2))).^(2/3);
end