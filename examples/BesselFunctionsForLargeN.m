close all
clear all %#ok

startup
resultsFolder = [folderName '/BesselFunctionsForLargN'];
if ~exist(resultsFolder, 'dir')
    mkdir(resultsFolder);
end

%% Create plot of bessel functinos for large n
xi = 500;
n = 900;
if false
    i_max = 10000;
    U_pol = cell(i_max,1);
    for i = 1:i_max
        U_pol{i} = U_p(i-1,U_pol);
    end
    save('../miscellaneous/U_pol.mat','U_pol')
else
    load('../miscellaneous/U_pol.mat')
end
% U_p(1,U_pol)-[-5,0,3,0]/24
% U_p(2,U_pol)-[385,0,-462,0,81,0,0]/1152
% U_p(3,U_pol)-[-425425,0,765765,0,-369603,0,30375,0,0,0]/414720
% bessel_j_asy(n,xi)
% bessel_s(n,xi,1)

type = 2;
switch type
    case 1
        N = 550;
    case 2
        N = 1050;
end
% N = 150;
npts = N;
B = zeros(npts,1);
N_arr = linspace(1,N,npts).';
B1 = bessel_s(N_arr,xi,1);
B2 = bessel_s(N_arr,xi,2);
% semilogy(N_arr,1./abs(B1),N_arr,abs(B2),N_arr,j_asm)
% semilogy(N_arr,abs(B1),N_arr,abs(B2),N_arr,abs(asc_j(N_arr,xi)),N_arr,abs(asc_y(N_arr,xi)))
% semilogy(N_arr,abs(B1),'blue',N_arr,abs(bessel_j_asy(N_arr,xi,U_pol)),'green')
% semilogy(N_arr,abs(B1),'blue',N_arr,abs((4*zeta(xi./(N_arr+0.5))./(1-(xi./(N_arr+0.5)).^2)).^(1/4)./(N_arr+0.5).^(1/3).*airy((N_arr+0.5).^(2/3).*zeta(xi./(N_arr+0.5)))),'green')
% semilogy(N_arr,abs(B1),'blue',N_arr,abs(airy((N_arr+0.5).^(2/3).*zeta(xi./(N_arr+0.5)))./(N_arr+0.5).^(1/3)),'green')
% semilogy(N_arr,abs((N_arr+0.5).^(2/3).*zeta(xi./(N_arr+0.5))),'green')
% plot(N_arr,abs((N_arr+0.5).^(2/3).*zeta(xi./(N_arr+0.5))),'green')
semilogy(N_arr,airy(0,N_arr,1),N_arr,Ai(N_arr),'green')
% semilogy(N_arr,abs(B1),'blue',N_-arr,abs(exp(-2/3*((N_arr+0.5).^(2/3).*zeta(xi./(N_arr+0.5))).^(3/2))),'green')
xlim([0,N])
xlabel('$n$','interpreter','latex')
ylabel('Magnitude of Bessel functions')
switch type
    case 1
        ylim([2e-10,5e3])
        savefig([resultsFolder '/Figure4a'])
    case 2
%         ylim([1e-223,1e223])
        savefig([resultsFolder '/Figure4b'])
end
% legend({'$$|\mathrm{j}_n(500)|$$',...
%         '$$|\mathrm{y}_n(500)|$$',...
%         'Asymptotic $$|\mathrm{j}_n(500)|$$',...
%         'Asymptotic $$|\mathrm{y}_n(500)|$$'},'interpreter','latex','location','northwest')
legend({'$$|\mathrm{j}_n(500)|$$',...
        'Asymptotic $$|\mathrm{j}_n(500)|$$'},'interpreter','latex','location','northwest')
% 
% switch type
%     case 1
%         ylim([1e-10 1e4])
%         printResultsToFile([resultsFolder '/besseljForLargeN1'], {'x', N_arr, 'y', abs(B1)})
%         printResultsToFile([resultsFolder '/besselyForLargeN1'], {'x', N_arr, 'y', abs(B2)})
%     case 2
%         ylim([1e-220 1e220])
%         printResultsToFile([resultsFolder '/besseljForLargeN2'], {'x', N_arr, 'y', abs(B1)})
%         printResultsToFile([resultsFolder '/besselyForLargeN2'], {'x', N_arr, 'y', abs(B2)})
% end

function j = bessel_j_asy(nu,z,U_pol)
nu = nu + 1/2;
K = 0;
sum1 = 0;
sum2 = 0;
for k = 0:K
    sum1 = sum1 + A_k(k,zeta(z./nu),U_pol)./nu.^(2*k);
    sum2 = sum2 + B_k(k,zeta(z./nu),U_pol)./nu.^(2*k);
end
sum1 = airy(nu.^(2/3).*zeta(z./nu))./nu.^(1/3).*sum1;
sum2 = airy(1,nu.^(2/3).*zeta(z./nu))./nu.^(5/3).*sum2;
% sum1 = 1./nu.^(1/3).*sum1;
% sum2 = 1./nu.^(5/3).*sum2;
j = (4*zeta(z./nu)./(1-(z./nu).^2)).^(1/4).*(sum1+sum2);
j = sqrt(pi./(2*z)).*j;
end

function A = A_k(k,z,U_pol)

A = 0;
for j = 0:2*k
    A = A + (3/2)^j*v_k(j)*zeta(z).^(-3*j/2).*U_k(2*k-j,(1-z.^2).^(-1/2),U_pol);
%     A = A + (3/2)^j*v_k(j)*1i^(3*j).*(-zeta(z)).^(-3*j/2).*U_k(2*k-j,(z.^2-1).^(-1/2),U_pol);
end
end

function B = B_k(k,z,U_pol)

B = 0;
for j = 0:(2*k+1)
    B = B + (3/2)^j*u_k(j)*zeta(z).^(-3*j/2).*U_k(2*k-j+1,(1-z.^2).^(-1/2),U_pol);
%     B = B + (3/2)^j*u_k(j)*1i^(3*j).*(-zeta(z)).^(-3*j/2).*U_k(2*k-j+1,(z.^2-1).^(-1/2),U_pol);
end
B = -1i*(-zeta(z)).^(-1/2).*B;
end

function U = U_k(k,z,U_pol)
U = polyval(U_pol{k+1},z);
end

function U = U_p(k,U_pol)
if k == 0
    U = 1;
else
    U = zeros(1,3*k+1);
    a = conv([-0.5,0,0.5,0,0], polyder(U_pol{k}));
    a = a(1:end-1);
    b = 1/8*polyint(conv([-5,0,1],U_pol{k}));
    U(1:length(a)) = a;
    U(1:length(b)) = U(1:length(b)) + b;
end
end

function j = asc_j(n,z)

j = sqrt(pi./(2*z)).*1./sqrt(2*pi*(n+0.5)).*(exp(1)*z./(2*(n+0.5))).^(n+0.5);

end

function y = asc_y(n,z)

y = -sqrt(pi./(2*z)).*sqrt(2./(pi*(n+0.5))).*(exp(1)*z./(2*(n+0.5))).^(-n-0.5);

end

function u = u_k(K)

u = 1;
for k = 1:K
    u = (6*k-5)*(6*k-3)*(6*k-1)/(2*k-1)/216/k*u;
end

end

function v = v_k(K)

if K == 0
    v = 1;
else
    v = (6*K+1)/(1-6*K)*u_k(K);
end

end

function zeta = zeta(z)
zeta = zeros(size(z));
z_i = z(z < 1);
zeta(z < 1) = (3/2*(log((1+sqrt(1-z_i.^2))./z_i)-sqrt(1-z_i.^2))).^(2/3);
z_i = z(z >= 1);
zeta(z >= 1) = -(3/2*(sqrt(z_i.^2-1)-asec(z_i))).^(2/3);
end

function A = Ai(z)

A = 0;
for n = 0:20
    A = A + (-1)^n*gamma(n+5/6)*gamma(n+1/6)*(3/4)^n./(2*pi*factorial(n)*z.^(3*n/2));
end
% A = A.*exp(-2/3*z.^(3/2))./(2*sqrt(pi)*z.^(1/4));
A = A./(2*sqrt(pi)*z.^(1/4));

end