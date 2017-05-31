close all
clear all %#ok

addpath ..
addpath ../utils
addpath ../models

startMatlabPool

%% Fender (1972) example
setFenderParameters
R_i = R_o - t; % Inner radius of shell
nFreqs = 3000;
k = linspace(32/nFreqs,32,nFreqs)';
omega = k*c_f(1);

d_vec = -[0,0,1].';

%%%%%%%%%
SPL_Fender0 = importdata('../models/Fender1972sfa/Fig2.csv');
SPL_Fender180 = importdata('../models/Fender1972sfa/Fig3.csv');
theta = 0;
options = struct('d_vec', d_vec, ...
                 'omega', omega, ...
                 'R_i', R_i, ...
                 'R_o', R_o, ...
                 'E', E, ...
                 'nu', nu, ...
                 'rho_s', rho_s, ...
                 'rho_f', rho_f, ...
                 'c_f', c_f);
if 0
    v = [0,0,R_o*cos(0)];
    f = @(k)-objFunc(k,options,v,c_f(1),1);
    specialValues = findExtremas(f, 2/nFreqs, 32, 100000)';
    v = [0,0,R_o*cos(pi)];
    f = @(k)-objFunc(k,options,v,c_f(1),1);
    specialValues = [specialValues; findExtremas(f, 2/nFreqs, 32, 100000)'];
    save('../results/Fender_extremas', 'specialValues')
    delta = 1e-4;
    specialValues = sort([specialValues; (specialValues-delta); (specialValues+delta)]);
else
    load('../results/Fender_extremas')
end
k = unique(sort([k; specialValues]));
omega = k*c_f(1);   % Wave number for outer fluid domain
options.omega = omega;
v = [0,0,R_o*cos(0);
     0,0,R_o*cos(pi)];
tic
data = e3Dss(v, options);
toc
p_inc = @(v) exp(1i*dot3(v,d_vec)*k.');
p_tot = data(1).p + p_inc(v);
SPL = 20*log10(abs(p_tot)); % sound pressure level
figure(2)
plot(k*R_o, SPL(1,:), SPL_Fender0(:,1), SPL_Fender0(:,2))
set(0,'defaulttextinterpreter','latex')
title('Predicted total pressure as a function of $$ka$$ at the surface of the shell, $$\theta = 0^\circ$$')
xlabel('$$k_Oa$$')
ylabel('Surface sound pressure level')
ylim([-80 120])
legend({'Present work', 'Reference Solution from Fender (1972) - Figure 2'})
xlim(R_o*[k(1) k(end)])

%%%%%%%%
figure(3)
plot(k*R_o, SPL(2,:), SPL_Fender180(:,1), SPL_Fender180(:,2))
set(0,'defaulttextinterpreter','latex')
title('Predicted total pressure as a function of $$ka$$ at the surface of the shell, $$\theta = 180^\circ$$')
xlabel('$$k_Oa$$')
ylabel('Surface sound pressure level')
ylim([-60 120])
legend({'Present work', 'Reference Solution from Fender (1972) - Figure 2'})
xlim(R_o*[k(1) k(end)])

%         
%         folderName = '../results';
%         if ~exist(folderName, 'dir')
%             mkdir(folderName);
%         end
%         
%         aspect = '0';
%         elevation = '90';
%         frequency = 'S';
%         BC = 'NNBC';
%         
%         varCol.alpha_s = 0;
%         varCol.beta_s = pi/2;
%         scatteringCase = 'Sweep';
%         model = 'Fender';
%         varCol.scatteringCase = scatteringCase;
%         
%         varCol.f_arr = 2*pi*omega;
%         
%         saveName = [model '_' BC '_' scatteringCase '_A'  aspect '_E' elevation '_F' frequency '_Fender1'];
%         varCol.saveName = saveName;
%         filename = [folderName '/' saveName];
%         printResultsToFile(filename, SPL_Fender0(:,1), SPL_Fender0(:,2), varCol, 1, 0, 'Fender', 'Results using WebPlotDigitizer to scan the results from Fender (1972)')
%         
%         
%         saveName = [model '_' BC '_' scatteringCase '_A'  aspect '_E' elevation '_F' frequency '_Fender2'];
%         varCol.saveName = saveName;
%         filename = [folderName '/' saveName];
%         printResultsToFile(filename, SPL_Fender180(:,1), SPL_Fender180(:,2), varCol, 1, 0, 'Fender', 'Results using WebPlotDigitizer to scan the results from Fender (1972)')
%         
%         
%         saveName = [model '_' BC '_' scatteringCase '_A'  aspect '_E' elevation '_F' frequency '_1'];
%         varCol.saveName = saveName;
%         filename = [folderName '/' saveName];
%         printResultsToFile(filename, k*a, SPL(1,:).', varCol, 1, 0, 'NTNU_FFI', 'Analytic solution')
%         
%         saveName = [model '_' BC '_' scatteringCase '_A'  aspect '_E' elevation '_F' frequency '_2'];
%         varCol.saveName = saveName;
%         filename = [folderName '/' saveName];
%         printResultsToFile(filename, k*a, SPL(2,:).', varCol, 1, 0, 'NTNU_FFI', 'Analytic solution')

if 0
    figure(42)
    nFreqs = 500;
    k = linspace(32/nFreqs,32,nFreqs)'; % wave number
    k = unique(sort([k; specialValues]));
    omega = k*c_f(1);   % Wave number for outer fluid domain
    options.omega = omega;
    createConvergencePlot('3D',options,v,75, '../../../LaTeX/createFigures/contents/e3Dss/FenderError')
    savefig('results/FenderError.fig')
end

% figure(4)
% nFreqs = 2000;
% k_max = 500;
% k = linspace(k_max/nFreqs,k_max,nFreqs);
% omega = k*c_f(1);   % Wave number for outer fluid domain
% options.omega = omega;
% v = [0,0,R_o*cos(0);
%      0,0,R_o*cos(pi)];
% tic
% data = e3Dss(v, options);
% toc
% plot(k*R_o,data(1).N_eps)
% xlabel('$k_1R_0$')
% ylabel('$N$')
% printResultsToFile('../results/Fender_Bound', k(~isnan(data(1).N_eps)).'*R_o, data(1).N_eps(~isnan(data(1).N_eps)), [], 0, 1)
