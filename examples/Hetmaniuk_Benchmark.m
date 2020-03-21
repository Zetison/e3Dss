%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This study correspond to Figure 8 and 17 in Hetmaniuk2012raa
% Hetmaniuk2012raa is available at https://doi.org/10.1002/nme.4271

close all
clear all %#ok

pathToResults = '../../../results/e3Dss/';
% pathToResults = '../results';

layer = setHetmaniukParameters();
BC = 'SHBC';
layer = layer(1:end-1);

k = (9:0.05:36)';
omega = k*layer{1}.c_f;
options.omega = omega;
options.BC = BC;
options.P_inc = -1; % Amplitude of incident wave

layer{1}.X = layer{1}.R_i*[0,0,1];
layer{1}.calc_p = true;
layer{1}.calc_p_0 = false;
layer = e3Dss(layer, options);

figure(1)
real_p_Hetmaniuk = importdata('../models/Hetmaniuk2012raa/Figure8.csv');
plot(k, real(layer{1}.p),'DisplayName','Exact')
title('Figure 17 in Hetmaniuk2012raa')
hold on
plot(real_p_Hetmaniuk(:,1), real_p_Hetmaniuk(:,2),'DisplayName','Hetmaniuk2012raa')
xlabel('Wavenumber')
xlim([5, 40])
ylim([-2, 2])
ylabel('Real part of pressure')  
legend('show');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

layer = setHetmaniukParameters();
BC = 'SSBC';

if false
    f = (1430:12:4290)';
    omega = 2*pi*f;   % Wave number for outer fluid domain
    k = omega/layer{1}.c_f;
else
    k = (6:0.05:18)';
    omega = k*layer{1}.c_f;
    f = omega/(2*pi);
end
d_vec = [0,0,1].';
options = struct('applyLoad', 'planeWave', ...
                 'd_vec', d_vec, ...
                 'BC', BC, ...
                 'omega', omega, ...
                 'P_inc', 1);
             

layer{1}.X = layer{1}.R_i*[0,0,-1];
layer{1}.calc_p = true;
layer{1}.calc_p_0 = false;
layer = e3Dss(layer, options);

figure(3)
real_p_Hetmaniuk = importdata('../models/Hetmaniuk2012raa/Figure17.csv');
plot(f, real(layer{1}.p),'DisplayName','Exact')
title('Figure 17 in Hetmaniuk2012raa')
hold on
plot(real_p_Hetmaniuk(:,1), real_p_Hetmaniuk(:,2),'DisplayName','Hetmaniuk2012raa')
xlabel('Frequency [Hz]')
xlim([f(1), f(end)])
ylim([-15, 15])
ylabel('Real part of pressure')  
legend('show');
