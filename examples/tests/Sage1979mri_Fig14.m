function tasks = Sage1979mri_Fig14(plotResults)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This study correspond to Fig. 1 and Fig. 4 in Sage1979mri
% Sage1979mri is available at https://doi.org/10.1121/1.382928

if nargin < 1
    plotResults = false;
end

startup
folderName = [homeDir '/Dropbox/Apps/Overleaf/createFigures/data/e3Dss_article2'];
resultsFolder = [folderName '/Sage1979mri'];
if ~exist(resultsFolder, 'dir')
    mkdir(resultsFolder);
end

%% Define parameters
nFreqs = 5000;
layer = setSage1979mriParameters();
layer{1}.calc_p_0 = true; % Calculate the far field pattern
R = layer{1}.R;

%% Calculate dependent parameters

kR = [10.^linspace(-3,0,nFreqs)'; linspace(1+5/nFreqs,5,nFreqs)'];
% k = 1;
c_f = layer{1}.c_f;
omega = kR*c_f;

d_vec = [0,0,1].';

%%%%%%%%%
%% Run simulation
options = struct('BC','NNBC',...
                 'P_inc', 1, ...
                 'd_vec', d_vec, ...
                 'Display', 'none', ...
                 'omega', omega);

layer{1}.X = [0,0,-1]; % Compute backscattered pressure

load('miscellaneous/Sage_extremas')

kR = unique(sort([kR; specialValues]));
options.omega = kR/layer{1}.R*layer{1}.c_f;
layer = e3Dss(layer, options);

%% Plot results
SPL_Sage1 = importdata('models/Sage1979mri/Figure1.csv');
SPL_Sage2 = importdata('models/Sage1979mri/Figure4.csv');
SPL_Sage = [SPL_Sage1(SPL_Sage1(:,1) < 0.4,1), SPL_Sage1(SPL_Sage1(:,1) < 0.4,2);
            SPL_Sage2(SPL_Sage2(:,1) > 0.4,1), SPL_Sage2(SPL_Sage2(:,1) > 0.4,2)];

sigma_s = 4*pi*abs(layer{1}.p_0).^2/abs(options.P_inc)^2;
if plotResults
    figure(1)   
    loglog(kR, sigma_s/(pi*R^2)/pi,'DisplayName','Present work')
    hold on
    set(0,'defaulttextinterpreter','latex')
    xlabel('$$k_1 a$$')
    ylabel('$$\frac{\sigma}{\pi a^2}\frac{1}{\pi}$$')
    xlim([kR(1) kR(end)])
    ylim([1e-4,1e4])

    figure(4)
    semilogy(kR, sigma_s/(pi*R^2)/pi,'DisplayName','Present work')
    hold on
    set(0,'defaulttextinterpreter','latex')
    xlabel('$$k_1 a$$')
    ylabel('$$\frac{\sigma}{\pi a^2}\frac{1}{\pi}$$')
    xlim([kR(1) kR(end)])
    ylim([1e-4,100])
end
tasks(1).sigma_s_scaled = sigma_s/(pi*R^2)/pi;
tasks(1).SPL_Sage = SPL_Sage;

layer = layer(1);
options.BC = 'SSBC';
layer = e3Dss(layer, options);

sigma_s = 4*pi*abs(layer{1}.p_0).^2/abs(options.P_inc)^2;
if plotResults
    colors = get(gca,'colororder');
    figure(1) 
    loglog(kR, sigma_s/(pi*R^2)/pi,'--','color','black','DisplayName','SSBC')
    savefig([resultsFolder '/figure1.fig'])

    figure(4) 
    semilogy(kR, sigma_s/(pi*R^2)/pi,'--','color','black','DisplayName','SSBC')
    savefig([resultsFolder '/figure4.fig'])


    figure(1) 
    hold on
    loglog(SPL_Sage(:,1),SPL_Sage(:,2),'color',colors(2,:),'DisplayName','Reference Solution from Sage (1979)')
    legend off
    legend show

    figure(4) 
    hold on
    semilogy(SPL_Sage(:,1),SPL_Sage(:,2),'color',colors(2,:),'DisplayName','Reference Solution from Sage (1979)')
    legend off
    legend show
end
tasks(2).sigma_s_scaled = sigma_s/(pi*R^2)/pi;
