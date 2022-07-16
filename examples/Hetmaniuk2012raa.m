%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This study correspond to Figure 8, 12 and 17 in Hetmaniuk2012raa
% Hetmaniuk2012raa is available at https://doi.org/10.1002/nme.4271

close all
clear all %#ok

startup
folderName = [homeDir '/Dropbox/Apps/Overleaf/createFigures/data/e3Dss_article2'];
resultsFolder = [folderName '/Hetmaniuk2012raa'];
if ~exist(resultsFolder, 'dir')
    mkdir(resultsFolder);
end

layer = setHetmaniukParameters();
BC = 'SHBC';
layer = layer(1:end-1);

k = (9:0.05:36)';
omega = k*layer{1}.c_f;
d_vec = [0,0,1].';
options = struct('applyLoad', 'planeWave', ...
                 'd_vec', d_vec, ...
                 'Display', 'final', ...
                 'BC', BC, ...
                 'omega', omega, ...
                 'P_inc', -1);

layer{1}.X = layer{1}.R*[0,0,1];
layer{1}.calc_p = true;
layer = e3Dss(layer, options);

figure(8)
real_p_Hetmaniuk = importdata('models/Hetmaniuk2012raa/Figure8.csv');
plot(k, real(layer{1}.p),'DisplayName','Exact')
title('Figure 8 in Hetmaniuk2012raa')
hold on
plot(real_p_Hetmaniuk(:,1), real_p_Hetmaniuk(:,2),'DisplayName','Hetmaniuk2012raa')
xlabel('Wavenumber')
xlim([5, 40])
ylim([-2, 2])
ylabel('Real part of pressure')  
legend('show');
savefig([resultsFolder '/Figure8'])

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
options = struct('applyLoad', 'mechExcitation', ...
                 'd_vec', d_vec, ...
                 'r_s', layer{2}.R, ...
                 'Display', 'final', ...
                 'BC', BC, ...
                 'omega', omega, ...
                 'nu_a', 100, ...
                 'saveRelTermMax', true, ... 
                 'P_inc', -1);

layer{1}.X = layer{1}.R*[0,0,1];
layer{1}.calc_p = true;
[layer,~,~,relTermMaxArr] = e3Dss(layer, options);

figure(12)
plot(f, real(layer{1}.p),'DisplayName','Exact')
title('Figure 12 in Hetmaniuk2012raa')
hold on
real_p_Hetmaniuk = importdata('models/Hetmaniuk2012raa/Figure12.csv');
plot(real_p_Hetmaniuk(:,1), real_p_Hetmaniuk(:,2),'DisplayName','Hetmaniuk2012raa')
xlabel('Frequency [Hz]')
xlim([f(1), f(end)])
ylim([-200, 200])
ylabel('Real part of pressure')  
legend('show');
savefig([resultsFolder '/Figure12'])
printResultsToFile([resultsFolder '/Figure12'], {'x', real_p_Hetmaniuk(:,1), 'y', real_p_Hetmaniuk(:,2), 'xlabel','f', 'ylabel','realp'})
printResultsToFile([resultsFolder '/Figure12_e3Dss'], {'x', f, 'y', real(layer{1}.p).', 'xlabel','f', 'ylabel','realp'})


figure(42)
skipTermsInPlot = 50;
for i = 1:floor(size(relTermMaxArr,1)/skipTermsInPlot)
    semilogy(0:(size(relTermMaxArr,2)-1),relTermMaxArr(skipTermsInPlot*i,:).','DisplayName',['\omega = ' num2str(omega(10*i))])
    hold on
end
xlabel('Number of terms in truncated series, N')
ylabel('Relative magnitude of term N w.r.t. the sum of the n first terms')
hold off
legend('show');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

layer = setHetmaniukParameters();
BC = 'SSBC';

options = struct('applyLoad', 'planeWave', ...
                 'd_vec', d_vec, ...
                 'Display', 'final', ...
                 'BC', BC, ...
                 'omega', omega, ...
                 'saveRelTermMax', true, ... 
                 'P_inc', 1);
             

layer{1}.X = layer{1}.R*[0,0,-1];
layer{1}.calc_p = true;
[layer,~,~,relTermMaxArr] = e3Dss(layer, options);

figure(17)
real_p_Hetmaniuk = importdata('models/Hetmaniuk2012raa/Figure17.csv');
plot(f, real(layer{1}.p),'DisplayName','Exact')
title('Figure 17 in Hetmaniuk2012raa')
hold on
plot(real_p_Hetmaniuk(:,1), real_p_Hetmaniuk(:,2),'DisplayName','Hetmaniuk2012raa')
xlabel('Frequency [Hz]')
xlim([f(1), f(end)])
ylim([-15, 15])
ylabel('Real part of pressure')  
legend('show');
savefig([resultsFolder '/Figure17'])
printResultsToFile([resultsFolder '/Figure17'], {'x', real_p_Hetmaniuk(:,1), 'y', real_p_Hetmaniuk(:,2), 'xlabel','f', 'ylabel','realp'})
printResultsToFile([resultsFolder '/Figure17_e3Dss'], {'x', f, 'y', real(layer{1}.p).', 'xlabel','f', 'ylabel','realp'})


figure(43)
for i = 1:floor(size(relTermMaxArr,1)/skipTermsInPlot)
    semilogy(relTermMaxArr(skipTermsInPlot*i,:).','DisplayName',['\omega = ' num2str(omega(10*i))])
    hold on
end
xlabel('Number of terms in truncated series, N')
ylabel('Relative magnitude of term N w.r.t. the sum of the n first terms')
hold off
legend('show');



