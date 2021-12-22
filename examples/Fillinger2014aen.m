%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This study correspond to Figure 4 in Fillinger2014aen
% Fillinger2014aen is available at https://repository.tudelft.nl/view/tno/uuid%3A5b611e92-39cf-4c3e-b889-61c7f8d2f1d4

close all
clear all %#ok

startup
resultsFolder = [folderName '/Fillinger2014aen'];
if ~exist(resultsFolder, 'dir')
    mkdir(resultsFolder);
end

P_inc = 1; % Amplitude of incident wave
c_f = 1500;
k10start = -1;
k10end = 2;
k = 10.^linspace(k10start,k10end,3000);
R = 1; % Outer radius of shell

alpha_s = 0;
beta_s = 0;

beta_f = beta_s;
alpha_f = alpha_s;
X = R*[cos(beta_f)*cos(alpha_f); cos(beta_f)*sin(alpha_f); sin(beta_f)*ones(size(alpha_f))]';

d_vec = -[cos(beta_s)*cos(alpha_s);
          cos(beta_s)*sin(alpha_s);
          sin(beta_s)];

omega = c_f*k; % Angular frequency
layer{1} = struct('media', 'fluid', ...
                  'X', X, ...
                  'R', R, ...
                  'c_f', c_f, ...
                  'rho', 1500);
layer{1}.calc_p_0 = true; % Calculate the far field pattern

options = struct('BC', 'SHBC', ...
                 'd_vec', d_vec,... 
                 'omega', omega, ...
                 'P_inc', P_inc);
             
             
layer = e3Dss(layer, options);
TS_inf = 20*log10(R/2);
TS = 20*log10(abs(layer{1}.p_0));

figure(1)
semilogx(k, TS-TS_inf,'DisplayName','Analytic')
ylim([-30,10])
hold on
k = 10.^linspace(k10start,-0.35,2);
TS = 20*log10(5/6*k.^2*R^3);
semilogx(k, TS-TS_inf,'green','DisplayName','Asymptotic as $k \to 0$')

k = 10.^linspace(0.2,k10end,2);
TS = 10*log10(R^2/4)*ones(size(k));
semilogx(k, TS-TS_inf,'red','DisplayName','Asymptotic as $k \to \infty$')


xlabel('$$k\cdot a$$','interpreter','latex')
ylabel('$$\mathrm{TS}-\mathrm{TS}_{\infty}$$','interpreter','latex')

legend('off');
legend('show','interpreter','latex');

savefig([resultsFolder '/Figure1'])
