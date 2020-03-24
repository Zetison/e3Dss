%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This study correspond to Figure 10.4 and Figure 10.5 in Skelton1997tao
% Skelton1997tao is available at https://doi.org/10.1142/9781848160750_0001

close all
clear all %#ok

%% Define parameters
R = 1;
t_steel = 0.02;
t_coating = 0.02;
theta = 180*pi/180;
theta_s = theta;

layer{1}.media = 'fluid';
layer{1}.R_i = R+t_coating;
layer{1}.rho = 1000;
layer{1}.c_f = 1500;
layer{1}.calc_p_0 = true; % Calculate the far field pattern

rho_c = 800;    % density of coating
eta = 0.1;  	% loss factor
c_l = 123*sqrt(1-1i*eta);
c_s = 33*sqrt(1-1i*eta);
E = c_s^2*rho_c*(3*c_l^2-4*c_s^2)/(c_l^2-c_s^2);
nu = (c_l^2-2*c_s^2)/(2*(c_l^2-c_s^2));
layer{2}.media = 'solid';
layer{2}.R_i = R;
layer{2}.rho = 800;
layer{2}.E = E; % 0.260e7
layer{2}.nu = nu; % 0.460

rho_s = 7700;   % density of steel
eta = 0.01;  	% loss factor
layer{3}.media = 'solid';
layer{3}.R_i = R-t_steel;
layer{3}.rho = rho_s;
layer{3}.E = 0.195e12*sqrt(1-1i*eta);
layer{3}.nu = 0.290*sqrt(1-1i*eta);

%% Calculate dependent parameters
npts = 2000;
f_max = 1000;
f = linspace(f_max/npts,1000,1000);

d_vec = [sin(theta_s),0,cos(theta_s)].';

options = struct('BC', 'SSBC', ...
                 'd_vec', d_vec, ...
                 'omega', 2*pi*f);

layer{1}.X = -layer{1}.R_i*[sin(theta),0,cos(theta)]; % Compute backscattered pressure

layerSSBC = layer([1,3]);
layerSSBC{1}.R_i = R;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run simulation
layerCoating = e3Dss(layer, options);
TS = 20*log10(abs(layerCoating{1}.p_0));
plot(f,TS,'DisplayName','Coating')
hold on

layerSSBC = e3Dss(layerSSBC, options);
TS = 20*log10(abs(layerSSBC{1}.p_0));
plot(f,TS,'DisplayName','SSBC')

options.BC = 'SHBC';
layerSHBC = e3Dss(layerSSBC(1), options);
TS = 20*log10(abs(layerSHBC{1}.p_0));

plot(f,TS,'DisplayName','Hard')
set(0,'defaulttextinterpreter','latex')
xlabel('Frequency (Hz)')
ylabel('Target strength')
xlim([0 f_max])
ylim([-30,20])

TS_SHBC = importdata('../models/Skelton1997tao/Figure10.4_SHBC.csv');
plot(TS_SHBC(:,1),TS_SHBC(:,2),'DisplayName','Ref SHBC')
TS_SSBC = importdata('../models/Skelton1997tao/Figure10.4_SSBC.csv');
plot(TS_SSBC(:,1),TS_SSBC(:,2),'DisplayName','Ref SSBC')
TS_Coating = importdata('../models/Skelton1997tao/Figure10.5.csv');
plot(TS_Coating(:,1),TS_Coating(:,2),'DisplayName','Ref Coating')



