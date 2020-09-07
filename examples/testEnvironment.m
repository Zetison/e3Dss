%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This study correspond to Figure 4 in Fillinger2014aen
% Fillinger2014aen is available at https://repository.tudelft.nl/view/tno/uuid%3A5b611e92-39cf-4c3e-b889-61c7f8d2f1d4

close all
clear all %#ok

startup

P_inc = 1; % Amplitude of incident wave
c_f = 1500;
k10start = -1;
k10end = 4;
k = 10.^linspace(k10start,k10end,3000);
k = 1;
R = 1; % Outer radius of shell

alpha_s = 0;
beta_s = 0;

beta_f = beta_s;
alpha_f = linspace(0,pi,1000);
X = 1.3*R*[cos(beta_f)*cos(alpha_f); cos(beta_f)*sin(alpha_f); sin(beta_f)*ones(size(alpha_f))]';

d_vec = -[cos(beta_s)*cos(alpha_s);
          cos(beta_s)*sin(alpha_s);
          sin(beta_s)];

omega = c_f*k; % Angular frequency
layer{1} = struct('media', 'fluid', ...
                  'X', X, ...
                  'R_i', R, ...
                  'c_f', c_f, ...
                  'rho', 1500);
layer{1}.calc_p_0 = false; % Calculate the far field pattern
layer{1}.calc_p = true;
layer{2} = struct('media', 'solid', ...
                  'R_i', 0.5, ...
                  'E', 207e9, ...
                  'nu', 0.3, ...
                  'rho', 7669);

options = struct('BC', 'SHBC', ...
                 'd_vec', d_vec,... 
                 'r_s', R, ...
                 'theta_s', [pi/4,3*pi/4], ...
                 'applyLoad','mechExcitation', ...
                 'omega', omega, ...
                 'P_inc', P_inc);
             
             
layer = e3Dss(layer, options);

figure(1)
plot(alpha_f, abs(layer{1}.p),'DisplayName','Analytic')


