close all
clear all %#ok

%% Define parameters
nFreqs = 1000;
a = 1;

layer{1}.media = 'fluid';
layer{1}.R_i = a;
layer{1}.rho = 1025;
layer{1}.c_f = 1531;
layer{1}.calc_p_0 = true; % Calculate the far field pattern

%% Calculate dependent parameters

k = 10.^linspace(-2,2,nFreqs)';
% k = 1;
c_f = layer{1}.c_f;
omega = k*c_f;

d_vec = [0,0,1].';

%%%%%%%%%
%% Run simulation
options = struct('BC','SHBC',...
                 'd_vec', d_vec, ...
                 'omega', omega);

layer{1}.X = [0,0,-1]; % Compute backscattered pressure

options.omega = k*layer{1}.c_f;
layerSHBC = e3Dss(layer, options);
options.BC = 'SSBC';
layerSSBC = e3Dss(layer, options);

%% Plot results
figure(1)   
loglog(k*a, abs(layerSHBC{1}.p_0),'DisplayName','SHBC')
hold on
loglog(k*a, abs(layerSSBC{1}.p_0),'DisplayName','SSBC')
set(0,'defaulttextinterpreter','latex')
xlabel('$$k_1 a$$')
ylabel('$$\frac{\sigma}{\pi a^2}$$')
xlim(a*[k(1) k(end)])
% ylim([1e-4,1e4])
options.BC = 'IBC';
for i = 10.^linspace(-2,2,10)
    c_f = layer{1}.c_f;
    rho = layer{1}.rho;
    options.z = i*rho*c_f*(1+1i);
    layerIBC = e3Dss(layer, options);
    loglog(k*a, abs(layerIBC{1}.p_0),'DisplayName',['z = ' num2str(options.z)])
end
legend('show','location','best')