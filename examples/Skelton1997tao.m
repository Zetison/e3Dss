%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This study correspond to Figure 10.4, Figure 10.5, Figure 10.6 and Figure 10.7 in Skelton1997tao
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
f = linspace(f_max/npts,f_max,npts);

d_vec = -[sin(theta_s),0,cos(theta_s)].';

options = struct('BC', 'SSBC', ...
                 'd_vec', d_vec, ...
                 'omega', 2*pi*f);

layer{1}.X = -options.d_vec.'; % Compute backscattered pressure

layerSSBC = layer([1,3]);
layerSSBC{1}.R_i = R;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run simulation
figure(1)
layerCoating = e3Dss(layer, options);
TS = 20*log10(abs(layerCoating{1}.p_0));
plot(f,TS,'DisplayName','e3Dss Coating')
hold on

layerSSBC = e3Dss(layerSSBC, options);
TS = 20*log10(abs(layerSSBC{1}.p_0));
plot(f,TS,'DisplayName','e3Dss SSBC')

options.BC = 'SHBC';
layerSHBC = e3Dss(layerSSBC(1), options);
TS = 20*log10(abs(layerSHBC{1}.p_0));

plot(f,TS,'DisplayName','e3Dss Hard')
set(0,'defaulttextinterpreter','latex')
xlabel('Frequency (Hz)')
ylabel('Target strength')
xlim([0 f_max])
ylim([-30,20])
title('Figure 10.4 and Figure 10.5')

TS_SHBC = importdata('../models/Skelton1997tao/Figure10.4_SHBC.csv');
plot(TS_SHBC(:,1),TS_SHBC(:,2),'DisplayName','Ref SHBC')
TS_SSBC = importdata('../models/Skelton1997tao/Figure10.4_SSBC.csv');
plot(TS_SSBC(:,1),TS_SSBC(:,2),'DisplayName','Ref SSBC')
TS_Coating = importdata('../models/Skelton1997tao/Figure10.5.csv');
plot(TS_Coating(:,1),TS_Coating(:,2),'DisplayName','Ref Coating')
legend show


figure(2)
npts = 2000;
f_max = 25e3;
f = linspace(f_max/npts,f_max,npts);
options = struct('BC', 'SSBC', ...
                 'd_vec', d_vec, ...
                 'omega', 2*pi*f);

layerSSBC = e3Dss(layerSSBC, options);
TS = 20*log10(abs(layerSSBC{1}.p_0));
plot(f,TS,'DisplayName','e3Dss')
TS_SSBC1 = importdata('../models/Skelton1997tao/Figure10.4_SSBC.csv');
TS_SSBC2 = importdata('../models/Skelton1997tao/Figure10.6.csv');
TS_SSBC = [TS_SSBC1(TS_SSBC1(:,1) < 1000,1), TS_SSBC1(TS_SSBC1(:,1) < 1000,2);
            TS_SSBC2(TS_SSBC2(:,1) > 1000,1), TS_SSBC2(TS_SSBC2(:,1) > 1000,2)];
hold on
plot(TS_SSBC(:,1),TS_SSBC(:,2),'DisplayName','Ref')
ylim([-30,20])
title('Figure 10.6')
xlabel('Frequency (Hz)')
ylabel('Target strength')
legend show

% this case is not perfectly reproducable due to lacking parameters
figure(3)
f_c = 15e3; % pulse center freq.
N = 2^12;       % **NOTE** it is not clear what Skelton uses for this parameter
T = 10*5e-3;    % **NOTE** it is not clear what Skelton uses for this parameter
B = N/T; % bandwidth
f_L = -B/2;
f_R = B/2;
df = 1/T;
f = linspace(0,f_R-df,N/2);
omega = 2*pi*f;
omega_c = 2*pi*f_c;
P_inc = 1;
layerSSBC{1}.X = -d_vec.'; % Compute backscattered pressure
plotP_inc = 1;

options = struct('BC', 'SSBC', ...
                 'd_vec', d_vec, ...
                 'omega', omega(2:end));
options.P_inc = @(omega) P_inc_(omega, omega_c,P_inc,3);

layerSSBC = e3Dss(layerSSBC, options);

npts = 1;
p_0 = zeros(npts,N/2);
for n = 0:N-1
    f_n = f_L + (f_R-f_L)/N*n;
    omega_n = 2*pi*f_n;
    if n >= N/2+1
        k = omega_n/layer{1}.c_f;
        k_vec = options.d_vec*k;
        p_0(:,n-N/2+1) = layerSSBC{1}.p_0(:,n-N/2);
    end
end
startIdx = round(0.9668*N);       % **NOTE** it is not clear what Skelton uses for this parameter
p_0_t = 2/T*fft(p_0,N,2);
temp = p_0_t;
p_0_t(:,1:N-startIdx+1) = temp(:,startIdx:end);
p_0_t(:,N-startIdx+2:end) = temp(:,1:startIdx-1);
dt = T/N;
t = dt*(0:N-1);
plot(1000*t,1000*real(p_0_t),'DisplayName','e3Dss')
p_ref = importdata('../models/Skelton1997tao/Figure10.7.csv');
hold on
plot(1000*p_ref(:,1),1000*p_ref(:,2),'DisplayName','Ref')
ylim([-600,600])
xlim([0,5])
title('Figure 10.7')
xlabel('Time (ms)')
ylabel('Pressure (mPa)')
legend show



