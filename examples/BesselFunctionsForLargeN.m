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
x = linspace(0,1,100000);
% v_k(2)
% plot(x, zeta(x), x, 1-2^(-1/3)*x+3/10*2^(-1/3)*x.^2+1/700*x.^3)
% plot(x, x-zeta(1-2^(-1/3)*x+3/10*2^(-1/3)*x.^2+1/700*x.^3))
% k_arr = 0:100;
% semilogy(k_arr,u_k(k_arr),k_arr,abs(v_k(k_arr)))
% return
% for k = 0:1 %numel(U_pol)-1
% %     plot(x, polyval(U_pol{k+1},x))
% %     plot(x, polyval(U_pol{k+1},(1-x.^2).^(-1/2)))
% %     plot(x, polyval(U_pol{k+1},(1-x.^2).^(-1/2)))
% %     plot(x, A_k(k,x,U_pol))
%     A_k(k,0,U_pol)
%     
%     ylim([-10,10])
%     hold on
% end
% return
% U_p(1,U_pol)-[-5,0,3,0]/24
% U_p(2,U_pol)-[385,0,-462,0,81,0,0]/1152
% U_p(3,U_pol)-[-425425,0,765765,0,-369603,0,30375,0,0,0]/414720
% bessel_j_asy(n,xi)
% bessel_s(n,xi,1)

% g = @(n) exp((n+0.5).*(log(n+0.5+sqrt((n+0.5).^2-xi.^2)./xi)-sqrt(1-(xi./(n+0.5)).^2)) ...
%              - (n+1.5).*(log(n+0.5+sqrt((n+1.5).^2-xi.^2)./xi)-sqrt(1-(xi./(n+1.5)).^2)));

type = 2;
switch type
    case 1
        N = 550;
    case 2
        N = 1050;
end
npts = N;
B = zeros(npts,1);
N_arr = linspace(1,N,npts).';
% semilogy(N_arr,abs(g(N_arr)));
% semilogy(k_arr,u_k(k_arr),k_arr,abs(v_k(k_arr)))
% semilogy(N_arr,abs(bessel_c(N_arr,xi,1)),N_arr,abs(besselj(N_arr,xi)))
% semilogy(N_arr,abs(bessel_c(N_arr,xi,2)),N_arr,abs(bessely(N_arr,xi)))
% semilogy(N_arr,abs(bessel_c(N_arr,xi,2)),N_arr,abs(bessely(N_arr,xi)))
% figure
semilogy(N_arr,abs(bessel_c(N_arr,xi,2)-bessely(N_arr,xi))./abs(bessely(N_arr,xi)))
hold on
semilogy(N_arr,abs(bessel_c(N_arr,xi,1)-besselj(N_arr,xi))./abs(besselj(N_arr,xi)))
if false
    xi_arr = linspace(1,1000,1000);
    C = zeros(numel(N_arr),numel(xi_arr));
    for n = N_arr.'
        n
        for j = 1:numel(xi_arr)
            xi = xi_arr(j);
            if n > (xi+32)/0.97 %xi < n
                C(n,j) = log10(abs(bessel_c(n,xi,2)-bessely(n,xi))./abs(bessely(n,xi)));
            end
        end
    end
    imagesc(xi_arr,N_arr,C)
    hold on
    plot(10:1000,((10:1000)+32)/0.97,'magenta') % n > (xi+30)/0.97
    xlabel('xi')
    ylabel('n')
    colorbar
    
    npts = 1000;
    z = linspace(1,1000,npts);
    N_arr = linspace(1,1000,npts).';
    nu = N_arr+0.5;
    C = zeros(numel(N_arr),numel(z));
    indices = z < nu;
    C(indices) = 1;
    indices = abs(nu.*zeta23_(z./nu)) > log(sqrt(realmax(class(z))));
    C(indices) = 2;
    indices = and(indices, z < nu);
    C(indices) = 3;
    imagesc(z,N_arr,C)
    hold on
%     plot(10:1000,((10:1000)+32)/0.97,'magenta') % n > (xi+30)/0.97
    xlabel('z')
    ylabel('n')
    colorbar
end
return
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
    sum1 = sum1 + A_k(k,z./nu,U_pol)./nu.^(2*k);
    sum2 = sum2 + B_k(k,z./nu,U_pol)./nu.^(2*k);
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
    U = conv([-0.5,0,0.5,0,0], polyder(U_pol{k}));
    U = U(1:end-1);
    U = U + 1/8*polyint(conv([-5,0,1],U_pol{k}));
end
end

function j = asc_j(n,z)

j = sqrt(pi./(2*z)).*1./sqrt(2*pi*(n+0.5)).*(exp(1)*z./(2*(n+0.5))).^(n+0.5);

end

function y = asc_y(n,z)

y = -sqrt(pi./(2*z)).*sqrt(2./(pi*(n+0.5))).*(exp(1)*z./(2*(n+0.5))).^(-n-0.5);

end

function u = u_k(k)
u = prod((2*k+1):2:(6*k-1))./216.^k./factorial(k);
end

function v = v_k(k)
v = (6*k+1)./(1-6*k).*u_k(k);
end

% function zeta = zeta(z)
% zeta = zeros(size(z));
% z_i = z(z < 1);
% zeta(z < 1) = (3/2*(log((1+sqrt(1-z_i.^2))./z_i)-sqrt(1-z_i.^2))).^(2/3);
% z_i = z(z >= 1);
% zeta(z >= 1) = -(3/2*(sqrt(z_i.^2-1)-asec(z_i))).^(2/3);
% end
function zeta = zeta(z)
zeta = (3/2*(log((1+sqrt(1-z.^2))./z)-sqrt(1-z.^2))).^(2/3);
end

function A = Ai(z)

A = 0;
for n = 0:20
    A = A + (-1)^n*gamma(n+5/6)*gamma(n+1/6)*(3/4)^n./(2*pi*factorial(n)*z.^(3*n/2));
end
% A = A.*exp(-2/3*z.^(3/2))./(2*sqrt(pi)*z.^(1/4));
A = A./(2*sqrt(pi)*z.^(1/4));

end