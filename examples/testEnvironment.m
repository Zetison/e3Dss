%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This study correspond to Figure 4 in Fillinger2014aen
% Fillinger2014aen is available at https://repository.tudelft.nl/view/tno/uuid%3A5b611e92-39cf-4c3e-b889-61c7f8d2f1d4

close all
clear all %#ok

startup

P_inc = 1; % Amplitude of incident wave
c_f = 1500;
k10start = -1;
k10end = 3;
k = 10.^linspace(k10start,k10end,1000);
% k = 100;
% k = 10;
R = 1; % Outer radius of shell

alpha_s = 240*pi/180;
beta_s = 30*pi/180;

beta_f = beta_s;
alpha_f = linspace(0,pi,1000);
alpha_f = 0.443;
r = 1.3*R;
% r = R;
X = r*[cos(beta_f)*cos(alpha_f); cos(beta_f)*sin(alpha_f); sin(beta_f)*ones(size(alpha_f))]';

d_vec = -[cos(beta_s)*cos(alpha_s);
          cos(beta_s)*sin(alpha_s);
          sin(beta_s)];

omega = c_f*k; % Angular frequency
layer{1} = struct('media', 'fluid', ...
                  'X', X, ...
                  'R_i', R, ...
                  'c_f', c_f, ...
                  'rho', 1500);
layer{1}.calc_p_0 = 1; % Calculate the far field pattern
layer{1}.calc_p = 0;
layer{2} = struct('media', 'fluid', ...
                  'R_i', 0,...
                  'c_f', c_f, ...
                  'rho', 100);
layer = setS135Parameters('double');
layer{1}.X = X;
layer{1}.calc_p_0 = 1; % Calculate the far field pattern
layer{1}.calc_p = 0;

options = struct('BC', 'NNBC', ...
                 'd_vec', d_vec,... 
                 'omega', omega, ...
                 'P_inc', P_inc);
             
             
[layer,N_eps,flag,relTermMaxArr] = e3Dss(layer, options);

figure(1)
options.nu_a = -1;
layer2 = e3Dss(layer, options);

if 1
    TS = 20*log10(abs(layer{1}.p_0));
    TS2 = 20*log10(abs(layer2{1}.p_0));
%     plot(alpha_f*180/pi, TS, alpha_f*180/pi, TS2)
    plot(k, TS, k, TS2)
    figure(2)
%     semilogy(alpha_f*180/pi, abs(layer2{1}.p_0-layer{1}.p_0)./abs(layer2{1}.p_0),'DisplayName','Analytic')
    semilogy(k, abs(layer2{1}.p_0-layer{1}.p_0)./abs(layer2{1}.p_0),'DisplayName','Analytic')
else
    layer{1}.p
    layer2{1}.p
%     err = abs(layer2{1}.p-layer{1}.p)./abs(layer2{1}.p)
end


