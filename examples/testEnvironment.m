close all
clear all %#ok

startup

layer = setHetmaniukParameters();
BC = 'SHBC';
layer = layer(1:end-1);

k = (9:0.05:36)';
omega = k*layer{1}.c_f;
d_vec = [0,0,1].';

theta_arr = linspace(0,pi,2000)'; % Set of angles used for plotting
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
omega = 1e2;
options = struct('applyLoad', 'surfExcitation', ...
                 'd_vec', d_vec, ...
                 'r_s', layer{1}.R_i, ...
                 'theta_s', pi/2+[0,10]*pi/180, ...
                 'BC', BC, ...
                 'omega', omega, ...
                 'P_inc', -1);
             

layer{2}.X = layer{1}.R_i*[sin(theta_arr), zeros(size(theta_arr)), cos(theta_arr)]; % Evaluate physical location of plotting points;
layer{2}.calc_sigma_s = [1,0,0,0,0,0];
layer = e3Dss(layer, options);

figure(12)
plot(theta_arr, imag(layer{2}.sigma_rr),'DisplayName','Exact')
title('Figure 12 in Hetmaniuk2012raa')
hold on
xlabel('Frequency [Hz]')
xlim([theta_arr(1), theta_arr(end)])
% ylim([-200, 200])
ylabel('Real part of pressure')  
legend('show');