%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This study correspond to Figure 7a, Figure 7b and Figure 8 in Venas2019e3s
% Venas2019e3s is available at https://doi.org/10.1016/j.cma.2018.02.015 (open access version at http://hdl.handle.net/11250/2493754)
% It is based on the example in Chang1994soa Figure 16 and Figure 17
% Chang1994soa is available at https://www.oden.utexas.edu/media/reports/1994/9412.pdf

close all
clear all %#ok

startup
resultsFolder = [folderName '/Chang1994soa'];
if ~exist(resultsFolder, 'dir')
    mkdir(resultsFolder);
end

%% Chang and Demkowiz (1994) example
P_inc = 1; % Amplitude of incident wave

%%%%%%%%%
layer = setChangParameters();
k = [15, 20];
a = layer{1}.R_i;
h = layer{1}.R_i-layer{2}.R_i;
omega = k*layer{1}.c_f;

d_vec = [0,0,1].';
p_inc = @(v) P_inc*exp(1i*dot3(v,d_vec)*k);
p_tot_Chang15 = importdata('../models/Chang1994voa2/Figure16.csv');
p_tot_Chang20 = importdata('../models/Chang1994voa2/Figure17.csv');
theta = linspace(0,pi,2000).';
phi = 0;
X = layer{1}.R_i*[sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)];
options = struct('d_vec', d_vec, ...
                 'BC', 'SSBC', ...
                 'omega', omega, ...
                 'P_inc', P_inc);
             
layer{1}.X     	= X;       % Evaluation points
layer{1}.calc_p = true;

[layer,N_eps] = e3Dss(layer,options); % Compute solution
p_tot = layer{1}.p + p_inc(X);

figure(16)
plot(theta*180/pi, real(p_tot(:,1)), p_tot_Chang15(:,1), p_tot_Chang15(:,2))
title(sprintf('$$h/a = %.2f, k = %.1f, N_{eps} = %d$$', h/a, k(1), N_eps(1)),'interpreter','latex')
xlabel('$\vartheta$','interpreter','latex')
ylabel('Real part of pressure [Pa]')
ylim([-2 2])
legend({'Present work', 'Reference Solution from Chang (1994)'},'Location','northwest')
xlim([0 180])
xtickformat('degrees')
savefig([resultsFolder '/Figure7a'])

figure(17)
plot(theta*180/pi, real(p_tot(:,2)), p_tot_Chang20(:,1), p_tot_Chang20(:,2))
title(sprintf('$$h/a = %.2f, k = %.1f, N_{eps} = %d$$', h/a, k(2), N_eps(2)),'interpreter','latex')
xlabel('$\vartheta$','interpreter','latex')
xtickformat('degrees')
ylabel('Real part of pressure [Pa]')
ylim([-2 2])
legend({'Present work', 'Reference Solution from Chang (1994)'})
xlim([0 180])
xtickformat('degrees')
savefig([resultsFolder '/Figure7b'])


folderName = '../results';
if ~exist(folderName, 'dir')
    mkdir(folderName);
end

aspect = '0';
elevation = 'm90';
waveNumber = num2str(15);
BC = 'SSBC';

varCol.alpha_s = 0;
varCol.beta_s = -90;
scatteringCase = 'BI';
model = 'Chang';
varCol.scatteringCase = scatteringCase;



saveName = [model '_' BC '_' scatteringCase '_A'  aspect '_E' elevation];
filename = [folderName '/' saveName];
filename = [];
figure(42)
createConvergencePlot('2D',layer,options,60,filename)
savefig([resultsFolder '/Figure8'])