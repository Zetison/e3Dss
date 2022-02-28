close all
clear all %#ok

startup
homeDir = expanduser('~');
folderName = [homeDir '/Dropbox/Apps/Overleaf/createFigures/data/e3Dss_article2'];
resultsFolder = [folderName '/BesselFunctionsForLargN'];
if ~exist(resultsFolder, 'dir')
    mkdir(resultsFolder);
end

%% Create plot of bessel functinos for large n
load('miscellaneous/U_pol_double.mat','U_pol','u_k','v_k')
xi = 500;
n = 900;
x = linspace(0,1,100000);
prec = class(xi);
negOneThirds = getC(prec,'-1/3');
oneThirds = getC(prec,'1/3');
twoThirds = getC(prec,'2/3');
% v_k(2)
% plot(x, zeta(x), x, 1-2^negOneThirds*x+3/10*2^negOneThirds*x.^2+1/700*x.^3)
% plot(x, x-zeta(1-2^negOneThirds*x+3/10*2^negOneThirds*x.^2+1/700*x.^3))
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
if 0
    x = 1000;
    nu_a = 100;
    nu = 100+1/2;
    
    j_n = bessel_c(nu,x,1,nu_a,U_pol,u_k,v_k);
    y_n = bessel_c(nu,x,2,nu_a,U_pol,u_k,v_k);
    h_n = bessel_c(nu,x,3,nu_a,U_pol,u_k,v_k);
    j_n/exp(exponent_(1,nu,x,nu_a)) + 1i*y_n/exp(exponent_(2,nu,x,nu_a))
    h_n/exp(exponent_(3,nu,x,nu_a))
    
    besselj(nu,x) + 1i*bessely(nu,x)
    besselh(nu,x)
    return
end
type = 2;
switch type
    case 1
        N = 550;
    case 2
        N = 1050;
end
npts = N;
nu_a = 100;
B = zeros(npts,1);
N_arr = linspace(1,N,npts).';
if 0
    N_arr = 10.^linspace(0,9,1000).';
%     N_arr = 10.^linspace(log10(xi),10,1000).';
    figure
    for x = 10.^(0:4)
        figure(45)
        g1 = @(n) exp( (n+0.5).*zeta23_(x./(n+0.5)) - (n+1.5).*zeta23_(x./(n+1.5)));
        scaledg1 = abs(g1(N_arr))./(x./(2*N_arr));
        semilogx(N_arr,scaledg1,'DisplayName',['$$g_n^{(1)}(x)/(x/(2n))$$ with $$x=10^' num2str(log10(x)) '$$']); 
        printResultsToFile([resultsFolder '/scaledg_' num2str(log10(x))], {'x', N_arr, 'y', scaledg1, 'xlabel','n', 'ylabel','scaledg'})
%         loglog(N_arr,abs(1-abs(g1(N_arr))./(x./(2*N_arr))),'DisplayName',['$$|1-g_n^{(1)}(x)/(x/(2n))|$$ with $$x=10^' num2str(log10(x)) '$$']); 
        hold on
        figure(46)
        g2 = @(n) exp( -abs(real((n+0.5).*zeta23_(x./(n+0.5)))) + abs(real((n+1.5).*zeta23_(x./(n+1.5)))));
        scaledg2 = abs(g2(N_arr))./(2*N_arr./x);
        semilogx(N_arr,scaledg2,'DisplayName',['$$g_n^{(2)}(x)/(2n/x)$$ with $$x=10^' num2str(log10(x)) '$$']); 
        printResultsToFile([resultsFolder '/scaledg_' num2str(log10(x))], {'x', N_arr, 'y', scaledg2, 'xlabel','n', 'ylabel','scaledg'})
        ylim([0,5])
%         loglog(N_arr,abs(1-abs(g1(N_arr))./(x./(2*N_arr))),'DisplayName',['$$|1-g_n^{(1)}(x)/(x/(2n))|$$ with $$x=10^' num2str(log10(x)) '$$']); 
        hold on
        figure(47)
        g3 = @(n) exp( -(n+0.5).*zeta23_(x./(n+0.5)) + (n+1.5).*zeta23_(x./(n+1.5)));
        scaledg3 = abs(g3(N_arr))./(2*N_arr./x);
        ylim([0,5])
        semilogx(N_arr,scaledg3,'DisplayName',['$$g_n^{(2)}(x)/(2n/x)$$ with $$x=10^' num2str(log10(x)) '$$']); 
        printResultsToFile([resultsFolder '/scaledg_' num2str(log10(x))], {'x', N_arr, 'y', scaledg3, 'xlabel','n', 'ylabel','scaledg'})
%         loglog(N_arr,abs(1-abs(g1(N_arr))./(x./(2*N_arr))),'DisplayName',['$$|1-g_n^{(1)}(x)/(x/(2n))|$$ with $$x=10^' num2str(log10(x)) '$$']); 
        hold on
    end
    for i = [45,46,47]
        figure(i)
        xlabel('$$n$$','interpreter','latex')
        legend('interpreter','latex')
    end
end
if 0
    n = 900;
    g1 = @(x) exp( (n+0.5).*zeta23_(x./(n+0.5)) - (n+1.5).*zeta23_(x./(n+1.5)));
    x = linspace(-1,1,100);
    y = x;
    [X,Y] = ndgrid(x,y);
    Z = X+1i*Y;
    imagesc(x,y,abs(g1(Z).'))
    set(gca,'YDir','normal')
    
    return
end
% semilogy(k_arr,u_k(k_arr),k_arr,abs(v_k(k_arr)))
% semilogy(N_arr,abs(bessel_c(N_arr,xi,1)),N_arr,abs(besselj(N_arr,xi)))
% semilogy(N_arr,abs(bessel_c(N_arr,xi,2)),N_arr,abs(bessely(N_arr,xi)))
% semilogy(N_arr,abs(bessel_c(N_arr,xi,2)),N_arr,abs(bessely(N_arr,xi)))
% figure
% semilogy(N_arr,abs(bessel_c(N_arr,xi,2,nu_a,U_pol,u_k,v_k)-bessely(N_arr,xi))./abs(bessely(N_arr,xi)))
% hold on
% semilogy(N_arr,abs(bessel_c(N_arr,xi,1,nu_a,U_pol,u_k,v_k)-besselj(N_arr,xi))./abs(besselj(N_arr,xi)))
if 0
    xi_arr = linspace(1,1000,1000);
    C = zeros(numel(N_arr),numel(xi_arr));
    for n = N_arr.'
        n
        for j = 1:numel(xi_arr)
            xi = xi_arr(j);
            if n > (xi+32)/0.97 %xi < n
                C(n,j) = log10(abs(bessel_c(n,xi,2,nu_a,U_pol,u_k,v_k)-bessely(n,xi))./abs(bessely(n,xi)));
            end
        end
    end
    imagesc(xi_arr,N_arr,C)
    hold on
    plot(10:1000,((10:1000)+32)/0.97,'magenta') % n > (xi+30)/0.97
    xlabel('xi')
    ylabel('n')
    colorbar
end
if 0
    npts = 1000;
    z = linspace(1,1000,npts);
    N_arr = linspace(1,1000,npts).';
    nu = N_arr+0.5;
    C = zeros(numel(N_arr),numel(z));
    indices = z < nu;
    C(indices) = 1;
    indices = abs(nu.*real(zeta23_(z./nu))) > log(sqrt(realmax(class(z))));
    C(indices) = 2;
    indices = indices_(nu,z,nu_a);
    C(indices) = 3;
    imagesc(z,N_arr,C)
    hold on
    xlabel('z')
    ylabel('n')
    colorbar
    return
end

if true
    % semilogy(N_arr,1./abs(B1),N_arr,abs(B2),N_arr,j_asm)
    % semilogy(N_arr,abs(B1),N_arr,abs(B2),N_arr,abs(asc_j(N_arr,xi)),N_arr,abs(asc_y(N_arr,xi)))
    % semilogy(N_arr,abs(B1),'blue',N_arr,abs(bessel_j_asy(N_arr,xi,U_pol)),'green')
    % semilogy(N_arr,abs(B1),'blue',N_arr,abs((4*zeta(xi./(N_arr+0.5))./(1-(xi./(N_arr+0.5)).^2)).^(1/4)./(N_arr+0.5).^oneThirds.*airy((N_arr+0.5).^twoThirds.*zeta(xi./(N_arr+0.5)))),'green')
    % semilogy(N_arr,abs(B1),'blue',N_arr,abs(airy((N_arr+0.5).^twoThirds.*zeta(xi./(N_arr+0.5)))./(N_arr+0.5).^oneThirds),'green')
    % semilogy(N_arr,abs((N_arr+0.5).^twoThirds.*zeta(xi./(N_arr+0.5))),'green')
    % plot(N_arr,abs((N_arr+0.5).^twoThirds.*zeta(xi./(N_arr+0.5))),'green')
    % semilogy(N_arr,airy(0,N_arr,1),N_arr,Airy(N_arr),'green')
    % semilogy(N_arr,abs(B1),'blue',N_-arr,abs(exp(-twoThirds*((N_arr+0.5).^twoThirds.*zeta(xi./(N_arr+0.5))).^(3/2))),'green')
    
    %% Create plot of bessel functinos for large n
    xi = 500;
    type = 1;
    nu_a = -1;
    switch type
        case 1
            N = 550;
        case 2
            N = 1050;
    end
    npts = N;
    N_arr = linspace(1,N,npts).';
    B1 = bessel_s(N_arr,xi,1,nu_a,U_pol,u_k,v_k,eps);
    B2 = bessel_s(N_arr,xi,2,nu_a,U_pol,u_k,v_k,eps);
    semilogy(N_arr,abs(B1),N_arr,abs(B2))
    xlim([0,N])
    xlabel('$n$','interpreter','latex')
    ylabel('Magnitude of Bessel functions')
    switch type
        case 1
            ylim([2e-10,5e3])
%             savefig([resultsFolder '/Figure4a'])
        case 2
            ylim([1e-223,1e223])
%             savefig([resultsFolder '/Figure4b'])
    end
    legend({'$$|\mathrm{j}_n(500)|$$','$$|\mathrm{y}_n(500)|$$'},'interpreter','latex','location','northwest')
%     switch type
%         case 1
%             ylim([1e-10 1e4])
%             printResultsToFile([resultsFolder '/besseljForLargeN1'], {'x', N_arr, 'y', abs(B1)})
%             printResultsToFile([resultsFolder '/besselyForLargeN1'], {'x', N_arr, 'y', abs(B2)})
%         case 2
%             ylim([1e-220 1e220])
%             printResultsToFile([resultsFolder '/besseljForLargeN2'], {'x', N_arr, 'y', abs(B1)})
%             printResultsToFile([resultsFolder '/besselyForLargeN2'], {'x', N_arr, 'y', abs(B2)})
%     end
end

function A = airy_dlmf(type,z,u_k)

A = 0;
N = 100;
twoThirds = getC(prec,'2/3');
zeta = twoThirds*z.^(3/2);
if type == 0
    for k = 0:N
        A_n = (-1)^k*u_k(k+1)./zeta.^k;
        A = A + A_n;
        if abs(A_n./A) < eps
            break
        end
    end
    if k == N
        error('Series did not converge')
    end
    A = A.*exp(-twoThirds*z.^(3/2))./(2*sqrt(pi)*z.^(1/4));
elseif type == 2
    for k = 0:N
        A_n = u_k(k+1)./zeta.^k;
        A = A + A_n;
        if abs(A_n./A) < eps
            break
        end
    end
    if k == N
        error('Series did not converge')
    end
    A = A.*exp(twoThirds*z.^(3/2))./(sqrt(pi)*z.^(1/4));
end

end

