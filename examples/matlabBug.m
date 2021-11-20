close all
a = 4.3;
z = a + 1i*a*sqrt(3)
% z = a + 1i*a
% z = a
zeta = 2/3*z^(3/2)
Ai_MATLAB = airy(0,z,0)
Ai_dlmf = airy_asympt(0,z)
Ai_as = airy_ascs(0,z)
% Ai_check = 1/pi*sqrt(z/3)*besselk(1/3,zeta)
Bi_MATLAB = airy(2,z,0)
Bi_dlmf = airy_asympt(2,z)
Bi_as = airy_ascs(2,z)
% Ai_check = sqrt(z/3)*(besseli(-1/3,zeta)+besseli(1/3,zeta))
expArg = 2/3*z^(3/2)
expFact = exp(expArg)

check_MATLAB = (exp(pi*1i/6)*airy(0,z*exp(2*pi*1i/3)) + exp(-pi*1i/6)*airy(0,z*exp(-2*pi*1i/3)) - airy(2,z))/airy(2,z)
check_dlmf = (exp(pi*1i/6)*airy_asympt(0,z*exp(2*pi*1i/3)) + exp(-pi*1i/6)*airy_asympt(0,z*exp(-2*pi*1i/3)) - airy_asympt(2,z))/airy_asympt(2,z)
check_as = (exp(pi*1i/6)*airy_ascs(0,z*exp(2*pi*1i/3)) + exp(-pi*1i/6)*airy_ascs(0,z*exp(-2*pi*1i/3)) - airy_ascs(2,z))/airy_ascs(2,z)


a = -linspace(0,40,1000);
z = a + 1i*a*sqrt(3);
semilogy(abs(z),abs(airy_asympt(0,z)),abs(z),abs(airy(0,z)),abs(z),abs(airy_ascs(0,z)),abs(z),abs(airy_asympt(2,z)),abs(z),abs(airy(2,z)),abs(z),abs(airy_ascs(2,z)))
legend('Ai Asymptotic expansion','Ai Asymptotic expansion (MATLAB)','Ai Ascending series', 'Bi Asymptotic expansion','Bi Asymptotic expansion (MATLAB)','Bi Ascending series')

type = 0;
b = 40;
x = linspace(-b*2,b*2,2000);
y = linspace(-b,b,1000);
[X,Y] = ndgrid(x,y);
Z = X+1i*Y;
if type == 0
    % imagesc(x,y,abs(exp(2/3.*Z.^(3/2))).')
    imagesc(x,y,log10(abs(exp(-2/3.*Z.^(3/2)))).')
else
    % imagesc(x,y,exp(-abs(2/3.*real(Z.^(3/2)))).')
    imagesc(x,y,log10(exp(abs(2/3.*real(Z.^(3/2))))).')
end
colorbar
axis equal
set(gca,'YDir','normal')
xlabel x
ylabel y
figure
imagesc(x,y,log10(abs(airy(type,Z))).')
colorbar
axis equal
set(gca,'YDir','normal')
xlabel x
ylabel y

function A = airy_asympt(type,z)

A = zeros(size(z));
N = 100;
zeta = 2/3*z.^(3/2);
indices = true(size(z));
if type == 0
    for k = 0:N
        A_k = (-1)^k*u_k(k)./zeta(indices).^k;
        A(indices) = A(indices) + A_k;
        indices(indices) = abs(A_k./A(indices)) > eps;
        if ~any(indices)
            break
        end
    end
    A = A.*exp(-2/3*z.^(3/2))./(2*sqrt(pi)*z.^(1/4));
elseif type == 2
    for k = 0:N
        A_k = u_k(k)./zeta(indices).^k;
        A(indices) = A(indices) + A_k;
        indices(indices) = abs(A_k./A(indices)) > eps;
        if ~any(indices)
            break
        end
    end
    A = A.*exp(2/3*z.^(3/2))./(sqrt(pi)*z.^(1/4));
end
if k == N
    warning('Series did not converge for %d out of %d elements.',sum(indices),numel(z))
end

end

function A = airy_ascs(type,z)

% Implementation of 10.4.2 and 10.4.3 in Abramowitz and Stegun
f = zeros(size(z));
g = zeros(size(z));
indices_f = true(size(z));
indices_g = true(size(z));
N = 100;
for k = 0:N
    if any(indices_f)
        alpha = 0;
        f_k = prod((3*alpha+1):3:(3*alpha+3*k-2))*z(indices_f).^(3*k)/factorial(3*k);
        f(indices_f) = f(indices_f) + f_k;
        indices_f(indices_f) = abs(f_k./f(indices_f)) > eps;
    end
    if any(indices_g)
        alpha = 1/3;
        g_k = prod((3*alpha+1):3:(3*alpha+3*k-2))*z(indices_g).^(3*k+1)/factorial(3*k+1);
        g(indices_g) = g(indices_g) + g_k;
        indices_g(indices_g) = abs(g_k./g(indices_g)) > eps;
    end
    if ~(any(indices_f) || any(indices_g))
        break
    end
end
if k == N
    warning('Series did not converge for %d out of %d elements.',sum(or(indices_f,indices_g)),numel(z))
end
c1 = 3^(-2/3)/gamma(2/3);
c2 = 3^(-1/3)/gamma(1/3);
if type == 0
    A = c1*f - c2*g;
elseif type == 2
    A = sqrt(3)*(c1*f + c2*g);
end

end

function u = u_k(k)
u = prod((2*k+1):2:(6*k-1))./216.^k./factorial(k);
end