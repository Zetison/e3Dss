%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This study correspond to Fig. 1 and Fig. 4 in Sage1979mri
% Sage1979mri is available at https://doi.org/10.1121/1.382928

close all
clear all %#ok

startup
resultsFolder = [folderName '/Sage1979mri'];
if ~exist(resultsFolder, 'dir')
    mkdir(resultsFolder);
end

%% Define parameters
nFreqs = 5000;
a = 1;

layer{1}.media = 'fluid';
layer{1}.R_i = a;
layer{1}.rho = 1025;
layer{1}.c_f = 1531;
layer{1}.calc_p_0 = true; % Calculate the far field pattern

layer{2}.media = 'fluid';
layer{2}.R_i = 0;
layer{2}.rho = 1.293;
layer{2}.c_f = 346.2;

%% Calculate dependent parameters

k = [10.^linspace(-3,0,nFreqs)'; linspace(1+5/nFreqs,5,nFreqs)'];
% k = 1;
c_f = layer{1}.c_f;
omega = k*c_f;

d_vec = [0,0,1].';

%%%%%%%%%
%% Run simulation
options = struct('BC','NNBC',...
                 'd_vec', d_vec, ...
                 'omega', omega);

layer{1}.X = [0,0,-1]; % Compute backscattered pressure

if 0
    f = @(k)-objFunc(k,layer,options);
    specialValues = findExtremas(f, k(1), k(end), 100000)';
    delta = 1e-5*k(end);
    specialValues = sort([specialValues; (specialValues-delta); (specialValues+delta)]);
    save('../miscellaneous/Sage_extremas', 'specialValues')
else
    load('../miscellaneous/Sage_extremas')
end
k = unique(sort([k; specialValues]));
options.omega = k*layer{1}.c_f;
layer = e3Dss(layer, options);

%% Plot results
SPL_Sage1 = importdata('../models/Sage1979mri/Figure1.csv');
SPL_Sage2 = importdata('../models/Sage1979mri/Figure4.csv');
SPL_Sage = [SPL_Sage1(SPL_Sage1(:,1) < 0.4,1), SPL_Sage1(SPL_Sage1(:,1) < 0.4,2);
            SPL_Sage2(SPL_Sage2(:,1) > 0.4,1), SPL_Sage2(SPL_Sage2(:,1) > 0.4,2)];
     
figure(1)   
sigma_s = 4*pi*abs(layer{1}.p_0).^2;
loglog(k*a, sigma_s/(pi*a^2)/pi,'DisplayName','Present work')
hold on
loglog(SPL_Sage(:,1),SPL_Sage(:,2),'DisplayName','Reference Solution from Sage (1979)')
warning('y axis scaled with another pi which is not consistent with Sage1979mri...')
set(0,'defaulttextinterpreter','latex')
xlabel('$$k_1 a$$')
ylabel('$$\frac{\sigma}{\pi a^2}$$')
xlim(a*[k(1) k(end)])
ylim([1e-4,1e4])

figure(4)
sigma_s = 4*pi*abs(layer{1}.p_0).^2;
semilogy(k*a, sigma_s/(pi*a^2)/pi,'DisplayName','Present work')
hold on
semilogy(SPL_Sage(:,1),SPL_Sage(:,2),'DisplayName','Reference Solution from Sage (1979)')
set(0,'defaulttextinterpreter','latex')
legend({'Present work', 'Reference Solution from Sage (1979)'})
xlabel('$$k_1 a$$')
ylabel('$$\frac{\sigma}{\pi a^2}$$')
xlim(a*[k(1) k(end)])
ylim([1e-4,100])

layer = layer(1);
options.BC = 'SSBC';
layer = e3Dss(layer, options);

figure(1) 
sigma_s = 4*pi*abs(layer{1}.p_0).^2;
loglog(k*a, sigma_s/(pi*a^2)/pi,'--','color','black','DisplayName','SSBC')
legend show
savefig([resultsFolder '/figure1.fig'])

figure(4) 
sigma_s = 4*pi*abs(layer{1}.p_0).^2;
semilogy(k*a, sigma_s/(pi*a^2)/pi,'--','color','black','DisplayName','SSBC')
legend off
legend show
savefig([resultsFolder '/figure4.fig'])

function sigma_s = objFunc(k,layer,options)

options.omega = k*layer{1}.c_f;
layer = e3Dss(layer, options);
sigma_s = 4*pi*abs(layer{1}.p).^2/abs(options.P_inc)^2;

end

