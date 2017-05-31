close all
clear all %#ok

addpath ..
addpath ../utils
addpath ../models

startMatlabPool
% In order to reproduce the plots in the paper: Comment away line 19-21 in
% bessel_s.m
setFenderParameters
R_i = R_o - t;
nFreqs = 2000;
k_max = 450;
k = [linspace(1e-300,k_max/nFreqs,200), linspace((k_max+1)/nFreqs,k_max,nFreqs)];
omega = k*c_f(1);   % Wave number for outer fluid domain

x = k*R_o(1);

K = E./(3*(1-2*nu));
G = E./(2*(1+nu));
c_s_1 = sqrt((3*K+4*G)./(3*rho_s));
c_s_2 = sqrt(G./rho_s);

Upsilon = min([R_i./c_s_1, R_i./c_s_2, R_o./c_f(1:end-1)]);

NN = 0:600;
N = zeros(size(x));
for i = 1:length(x)
    temp = abs(bessel_s(NN,omega(i)*Upsilon,2)) - 10^290;
    indices = find(temp > 0);
    N(i) = NN(indices(1));
end
hold on
plot(x,N)

xlabel('x')
ylabel('N')
legend({'Plot of $$N$$ as a function of $$x$$ for when $$y_N(x) = 10^{290}$$'}, 'interpreter','latex')
% printResultsToFile('../results/Fender_Bessely_Bound', x.', N.', [], 0, 1)
