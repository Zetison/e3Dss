%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This study correspond to Figure 7a, Figure 7b and Figure 8 in Venas2019e3s
% Venas2019e3s is available at https://doi.org/10.1016/j.cma.2018.02.015 (open access version at http://hdl.handle.net/11250/2493754)
% It is based on the example in Chang1994soa Figure 16 and Figure 17
% Chang1994soa is available at https://www.oden.utexas.edu/media/reports/1994/9412.pdf

close all
clear all %#ok

startup
addpath ../ASIGA/integration % Get the ASIGA toolbox here: https://github.com/Zetison/ASIGA

%% Chang and Demkowiz (1994) example
P_inc = 1; % Amplitude of incident wave

%%%%%%%%%
BC = 'SSBC';
% BC = 'SHBC';
if strcmp(BC,'SSBC')
    noDomains = 2;
else
    noDomains = 1;
end
layer = setChangParameters(noDomains);
% layer{1}.c_f = 1460; % 1460
% layer{1}.rho = 1000; % 1000
% 
% layer{2}.E = 2.0e10; % 2.0e11
% layer{2}.nu = 0.1;   % 0.3
% layer{2}.rho = 780; % 7800

R_1 = layer{1}.R;
k = [10/R_1, 20/R_1];
omega = k*layer{1}.c_f;

d_vec = [0,0,1].';
p_inc = @(v) P_inc*exp(1i*dot3(v,d_vec)*k);
theta = linspace(0,pi,2000).';
[Q,W] = gaussLegendreQuad(64);
theta = (Q+1)*pi/2;
% theta = 0;
phi = 0;
X = layer{1}.R*[sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)];
options = struct('d_vec', d_vec, ...
                 'BC', BC, ...
                 'omega', omega, ...
                 'saveRelTermMax', false, ... 
                 'Display', 'none', ... 
                 'nu_a', -1, ...
                 'P_inc', P_inc);
             
layer{1}.X     	= X;       % Evaluation points
layer{1}.calc_p = true;

warning('off', 'e3Dss:N_max_reached')


if 0
%     k = 10/R_1;
    k = 20/R_1;
    omega = k*layer{1}.c_f;
    options.omega = omega;
    theta = 0;
    Eps = eps;
    omega = options.omega;
    nu_a = -1;
    load('miscellaneous/U_pol_double.mat','U_pol','u_k','v_k')
    layer{1}.theta     	= 0.1; 
    layer{1}.X = 1.1*layer{1}.R*[sin(layer{1}.theta), 0, cos(layer{1}.theta)];
    N_arr = 0:190;
    d1 = zeros(size(N_arr));
    d2 = zeros(size(N_arr));
    frac = zeros(size(N_arr));
    frac2 = zeros(size(N_arr));
    
    [Z_zeta,Z_zeta_1,Z_xi_1,Z_eta_1,Z_xi_2,Z_eta_2] = allocateArrays(layer);
    counter = 1;
    for n = N_arr
        [Z_zeta,Z_zeta_1,Z_xi_1,Z_eta_1,Z_xi_2,Z_eta_2] = iterate_Zs(layer,options,n,Z_zeta,Z_zeta_1,Z_xi_1,Z_eta_1,Z_xi_2,Z_eta_2,nu_a,U_pol,u_k,v_k,Eps);
        [A1n3,d1(counter),d2(counter),frac(counter),frac2(counter)] = A1n3_(n,layer,options,Z_zeta_1,Z_xi_1,Z_eta_1,Z_xi_2,Z_eta_2);
        counter = counter + 1;
    end
    % semilogy(N_arr,abs(d1),N_arr,abs(d2),N_arr,abs(frac2))
    semilogy(N_arr,abs(frac2))
    hold on
    semilogy(ones(2,1).*ceil(2*k*layer{1}.R+eps),ylim.'.*ones(1,2),'--','color','black')
    semilogy(xlim.'.*ones(1,2),ones(2,1).*1/sqrt(3),'--','color','black')
    return
    figure
    semilogy(N_arr,abs(d1),N_arr,abs(d2))
    semilogy(N_arr,abs(d2./d1))
    % p = p_es(layer,options)
    % layer = e3Dss(layer, options);
    % p = layer{1}.p
    
    return
end
noRuns = 60;
% noRuns = round(k(1)*a+2);
layer = e3Dss(layer, options);
count = 1;
Error = zeros(noRuns,2);
N_arr = 0:noRuns-1;
p = layer{1}.p;
for N = N_arr
    options.N_max = N;
    layer = e3Dss(layer, options);
    p_N = layer{1}.p;
    Error(count,:) = sqrt(2*pi*sum(abs(p - p_N).^2.*sin(theta).*W,1));
    fprintf('Completed %d out of %d runs\n',count,noRuns)
    count = count + 1;
end

% x = k*layer{1}.R;
digits(100)
x = k*vpa(layer{1}.R);
xr = k.*vpa(norm(X(1,:)));
PI = getC(class(x),'pi');
third = getC(class(x),'1/3');
sbesselsquared = zeros(noRuns-1,2,class(x));
for n = N_arr(1:end-1) % sum from n=0 to N-1
    j_spherical = sqrt(PI./(2*x)).*besselj(n+1/2,x);
    j_spherical2 = sqrt(PI./(2*x)).*besselj(n+1/2+1,x);
    dj_spherical = n./x.*j_spherical - j_spherical2;
    sbesselsquared_n = (2*n+1)*dj_spherical.^2;
    if n == 0
        sbesselsquared(n+1,:) = sbesselsquared_n;
    else
        sbesselsquared(n+1,:) = sbesselsquared(n,:) + sbesselsquared_n;
    end
end
errorBound = double(sqrt(4/3*P_inc^2*PI*abs(third - sbesselsquared)));
close all
if 0
    figure
    err_diff2 = zeros(noRuns-1,size(k,2),class(PI));
    err_diff3 = zeros(noRuns-1,size(k,2),class(PI));
    for N = N_arr(2:end)
        for n = N:50 % sum from n=N to inf
            j_spherical = sqrt(PI./(2*x)).*besselj(n+1/2,x);
            j_spherical2 = sqrt(PI./(2*x)).*besselj(n+1/2+1,x);
            y_spherical = sqrt(PI./(2*x)).*bessely(n+1/2,x);
            y_spherical2 = sqrt(PI./(2*x)).*bessely(n+1/2+1,x);
            dj_spherical = n./x.*j_spherical - j_spherical2;
            dy_spherical = n./x.*y_spherical - y_spherical2;
            dhn = dj_spherical + 1i*dy_spherical;
            j_spherical_r = sqrt(PI./(2*xr)).*besselj(n+1/2,xr);
            y_spherical_r = sqrt(PI./(2*xr)).*bessely(n+1/2,xr);
            hn = j_spherical_r + 1i*y_spherical_r;
            err_diff2(N,:) = err_diff2(N,:) + (2*n+1)*abs(dj_spherical.*hn./dhn).^2;
            err_diff3(N,:) = err_diff3(N,:) + (2*n+1)*abs(dj_spherical).^2;
        end
        fprintf('Completed %d out of %d runs\n',N,noRuns-1)
    end
    errorBound2 = double(sqrt(4*P_inc^2*PI*err_diff2));
    errorBound3 = double(sqrt(P_inc^2*PI*err_diff2));
    semilogy(N_arr(2:end), errorBound2(:,1))
    semilogy(N_arr(2:end), errorBound3(:,1))
    return
end
figure
% semilogy(N_arr, Error(:,1),N_arr, Error(:,2))
semilogy(N_arr, Error(:,1),'color','blue')
hold on
semilogy(N_arr, Error(:,2),'color','red')
% semilogy(N_arr(2:end), errorBound(:,1),N_arr(2:end), errorBound(:,2))
semilogy(N_arr(2:end), errorBound(:,1),'--','color','blue')
semilogy(N_arr(2:end), errorBound(:,2),'--','color','red')
% semilogy(N_arr(2:end), errorBound2(:,1),N_arr(2:end), errorBound2(:,2))
% semilogy(N_arr(2:end), errorBound(:,1))
% semilogy(ones(2,1).*ceil(x*2+eps),ylim.'.*ones(1,2),'--','color','black')
validN = [ceil(x(1)*2+eps),53];
semilogy(ones(2,1).*validN(1),ylim.'.*ones(1,1),':','color','blue')
semilogy(ones(2,1).*validN(2),ylim.'.*ones(1,1),':','color','red')
x = double(x);

legend({['$$k_1R_1 = ' num2str(x(1)) '$$'], ['$$k_1R_1 = ' num2str(x(2)) '$$'], ...
        ['Bound for $$k_1R_1 = ' num2str(x(1)) '$$'], ['Bound for $$k_1R_1 = ' num2str(x(2)) '$$'], ...
        '$$N = 2\cdot 10$$', '$$N = 53$$'},'interpreter','latex')
xlabel('$N$','interpreter','latex')
ylabel('$\| p-p^{(N)}\|_{L^2(\Gamma_{\tilde{r}})}$','interpreter','latex')
set(gcf,'color','w');
% export_fig('/home/zetison/Dropbox/Apps/Overleaf/e3Dss_article3/graphics/Chang1994soa_2', '-pdf')

% legend({['$$k_1 = ' num2str(k(1)) '\mathrm{m}^{-1}$$'], '$$k_1 = 20\mathrm{m}^{-1}$$', ...
%         ['Bound for $$k_1 = ' num2str(k(1)) '\mathrm{m}^{-1}$$'], 'Bound for $$k_1 = 20\mathrm{m}^{-1}$$', ...
%         ['Bound2 for $$k_1 = ' num2str(k(1)) '\mathrm{m}^{-1}$$'], 'Bound2 for $$k_1 = 20\mathrm{m}^{-1}$$', ...
%         ['Bound3 for $$k_1 = ' num2str(k(1)) '\mathrm{m}^{-1}$$'], 'Bound3 for $$k_1 = 20\mathrm{m}^{-1}$$'},'interpreter','latex')
% legend({['$$k_1 = ' num2str(k(1)) '\mathrm{m}^{-1}$$'], ...
%         ['Bound for $$k_1 = ' num2str(k(1)) '\mathrm{m}^{-1}$$'], ...
%         ['Bound2 for $$k_1 = ' num2str(k(1)) '\mathrm{m}^{-1}$$'], ...
%         ['Bound3 for $$k_1 = ' num2str(k(1)) '\mathrm{m}^{-1}$$']},'interpreter','latex')

function p = p_es(layer,options)
nu_a = -1;
load('miscellaneous/U_pol_double.mat','U_pol','u_k','v_k')
Eps = eps;
omega = options.omega;
c_f = layer{1}.c_f;
k_1 = omega/c_f;


xi_1 = k_1*norm2(layer{1}.X);

[Z_zeta,Z_zeta_1,Z_xi_1,Z_eta_1,Z_xi_2,Z_eta_2,P,dP,d2P,p] = allocateArrays(layer);

N_max = 60;
for n = 0:N_max
    [Z_zeta,Z_zeta_1,Z_xi_1,Z_eta_1,Z_xi_2,Z_eta_2] = iterate_Zs(layer,options,n,Z_zeta,Z_zeta_1,Z_xi_1,Z_eta_1,Z_xi_2,Z_eta_2,nu_a,U_pol,u_k,v_k,Eps);

    A1n3 = A1n3_(n,layer,options,Z_zeta_1,Z_xi_1,Z_eta_1,Z_xi_2,Z_eta_2);
    [P, dP, d2P] = legendreDerivs(n, cos(layer{1}.theta), P, dP, d2P);
    hn = bessel_s(n,xi_1,3,nu_a,U_pol,u_k,v_k,Eps);
    p = p + A1n3*P(2,:).*hn;
end
end

function [Z_zeta,Z_zeta_1,Z_xi_1,Z_eta_1,Z_xi_2,Z_eta_2] = iterate_Zs(layer,options,n,Z_zeta,Z_zeta_1,Z_xi_1,Z_eta_1,Z_xi_2,Z_eta_2,nu_a,U_pol,u_k,v_k,Eps)
omega = options.omega;
c_f = layer{1}.c_f;
k_1 = omega/c_f;
xi_1 = k_1*norm2(layer{1}.X);
E_2 = layer{2}.E;
nu_2 = layer{2}.nu;
R_1 = layer{1}.R;
R_2 = layer{2}.R;
K_2 = E_2/(3*(1-2*nu_2));
G_2 = E_2/(2*(1+nu_2));
rho_2 = layer{2}.rho;
c_s12 = sqrt((3*K_2+4*G_2)/(3*rho_2));
c_s22 = sqrt(G_2/rho_2);
a_2 = omega/c_s12;
b_2 = omega/c_s22;
xi_11 = k_1*R_1;
xi_21 = a_2*R_1;
xi_22 = a_2*R_2;
eta_21 = b_2*R_1;
eta_22 = b_2*R_2;

besselIndices = [true,false,true];
Z_zeta = iterate_Z(n,xi_1,Z_zeta,1,besselIndices,nu_a,U_pol,u_k,v_k,Eps);
Z_zeta_1 = iterate_Z(n,xi_11,Z_zeta_1,1,besselIndices,nu_a,U_pol,u_k,v_k,Eps);
besselIndices = [true,true,false];
Z_xi_1 = iterate_Z(n,xi_21,Z_xi_1,1,besselIndices,nu_a,U_pol,u_k,v_k,Eps);
Z_eta_1 = iterate_Z(n,eta_21,Z_eta_1,1,besselIndices,nu_a,U_pol,u_k,v_k,Eps);
Z_xi_2 = iterate_Z(n,xi_22,Z_xi_2,1,besselIndices,nu_a,U_pol,u_k,v_k,Eps);
Z_eta_2 = iterate_Z(n,eta_22,Z_eta_2,1,besselIndices,nu_a,U_pol,u_k,v_k,Eps);
end


function [Z_zeta,Z_zeta_1,Z_xi_1,Z_eta_1,Z_xi_2,Z_eta_2,P,dP,d2P,p] = allocateArrays(layer)
Z_xi_1 = cell(3,2);
Z_eta_1 = cell(3,2);
Z_xi_2 = cell(3,2);
Z_eta_2 = cell(3,2);
Z_zeta = cell(3,2);
Z_zeta_1 = cell(3,2);
for m = 1:2
    switch layer{m}.media
        case 'solid'
            for i = 1:2
                for j = 1:2
                    Z_xi_1{i,j} = 0;
                    Z_eta_1{i,j} = 0;
                    Z_xi_2{i,j} = 0;
                    Z_eta_2{i,j} = 0;
                end
            end
        case 'fluid'
            i = 3;
            n_X = size(layer{m}.X,1);
            P = zeros(2,n_X); 
            dP = zeros(2,n_X); 
            d2P = zeros(2,n_X);
            p = zeros(1,n_X);
            for j = 1:2
                Z_zeta{i,j} = 0;
                Z_zeta_1{i,j} = 0;
            end
    end
end
end

function [A1n3,d1,d2,frac,frac2] = A1n3_(n,layer,options,Z_zeta_1,Z_xi_1,Z_eta_1,Z_xi_2,Z_eta_2)
gxi = {1,1};
geta = {1,1};
gzeta = {1,1,1};
omega = options.omega;
P_inc = options.P_inc;
c_f = layer{1}.c_f;
a_1 = omega/c_f;
E_2 = layer{2}.E;
nu_2 = layer{2}.nu;
R_1 = layer{1}.R;
R_2 = layer{2}.R;
K_2 = E_2/(3*(1-2*nu_2));
G_2 = E_2/(2*(1+nu_2));
rho_1 = layer{1}.rho;
rho_2 = layer{2}.rho;
c_s12 = sqrt((3*K_2+4*G_2)/(3*rho_2));
c_s22 = sqrt(G_2/rho_2);
a_2 = omega/c_s12;
b_2 = omega/c_s22;

xi_11 = a_1*R_1;
xi_21 = a_2*R_1;
xi_22 = a_2*R_2;
eta_21 = b_2*R_1;
eta_22 = b_2*R_2;

if n == 0
    H = [S_(1, 1, n, xi_21, eta_21, Z_xi_1,gxi), S_(1, 2, n, xi_21, eta_21, Z_xi_1,gxi);
         S_(5, 1, n, xi_21, eta_21, Z_xi_1,gxi), S_(5, 2, n, xi_21, eta_21, Z_xi_1,gxi);
         S_(5, 1, n, xi_22, eta_22, Z_xi_2,gxi), S_(5, 2, n, xi_22, eta_22, Z_xi_2,gxi)];
else
    H = [S_(1, 1, n, xi_21, eta_21, Z_xi_1,gxi), S_(1, 2, n, xi_21, eta_21, Z_xi_1,gxi), T_(1, 1, n, eta_21, Z_eta_1,geta), T_(1, 2, n, eta_21, Z_eta_1,geta);
         S_(5, 1, n, xi_21, eta_21, Z_xi_1,gxi), S_(5, 2, n, xi_21, eta_21, Z_xi_1,gxi), T_(5, 1, n, eta_21, Z_eta_1,geta), T_(5, 2, n, eta_21, Z_eta_1,geta);
         S_(7, 1, n, xi_21, eta_21, Z_xi_1,gxi), S_(7, 2, n, xi_21, eta_21, Z_xi_1,gxi), T_(7, 1, n, eta_21, Z_eta_1,geta), T_(7, 2, n, eta_21, Z_eta_1,geta);
         S_(7, 1, n, xi_22, eta_22, Z_xi_2,gxi), S_(7, 2, n, xi_22, eta_22, Z_xi_2,gxi), T_(7, 1, n, eta_22, Z_eta_2,geta), T_(7, 2, n, eta_22, Z_eta_2,geta);
         S_(5, 1, n, xi_22, eta_22, Z_xi_2,gxi), S_(5, 2, n, xi_22, eta_22, Z_xi_2,gxi), T_(5, 1, n, eta_22, Z_eta_2,geta), T_(5, 2, n, eta_22, Z_eta_2,geta)];
end
d1 = det(H(2:end,:));
d2 = det(H([1,3:end],:));
sdjn = dbessel_s(n,xi_11,1,Z_zeta_1,true,gzeta);
jn = Z_zeta_1{1,1};
sdhn = dbessel_s(n,xi_11,3,Z_zeta_1,true,gzeta);
hn = Z_zeta_1{3,1};

frac =   (2*G_2*sdjn*d1 + rho_1*omega^2*R_1^2*jn*d2)...
       ./(2*G_2*sdhn*d1 + rho_1*omega^2*R_1^2*hn*d2);
frac2 =   (2*G_2*xi_11*hn*d1 + rho_1*omega^2*R_1^2*hn*jn/sdjn*d2)...
       ./(2*G_2*sdhn*d1 + rho_1*omega^2*R_1^2*hn*d2);
A1n3 = -P_inc*1i^n*(2*n+1)*frac;
end



