function tasks = Hetmaniuk2012raa_Fig81217(plotResults)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This study correspond to Figure 8, 12 and 17 in Hetmaniuk2012raa
% Hetmaniuk2012raa is available at https://doi.org/10.1002/nme.4271

if nargin < 1
    plotResults = false;
end

startup
folderName = [homeDir '/Dropbox/Apps/Overleaf/createFigures/data/e3Dss_article2'];
resultsFolder = [folderName '/Hetmaniuk2012raa'];
if ~exist(resultsFolder, 'dir')
    mkdir(resultsFolder);
end

layer = setHetmaniukParameters(1);
BC = 'SHBC';

k = (9:0.05:36)';
omega = k*layer{1}.c_f;
d_vec = [0,0,1].';
options = struct('applyLoad', 'planeWave', ...
                 'd_vec', d_vec, ...
                 'Display', 'none', ...
                 'BC', BC, ...
                 'omega', omega, ...
                 'P_inc', -1);

layer{1}.X = layer{1}.R*[0,0,1];
layer{1}.calc_p = true;
layer = e3Dss(layer, options);

real_p_Hetmaniuk = importdata('models/Hetmaniuk2012raa/Figure8.csv');
if plotResults
    figure(8)
    plot(k, real(layer{1}.p),'DisplayName','Exact')
    title('Figure 8 in Hetmaniuk2012raa')
    hold on
    plot(real_p_Hetmaniuk(:,1), real_p_Hetmaniuk(:,2),'DisplayName','Hetmaniuk2012raa')
    xlabel('Wavenumber')
    xlim([5, 40])
    ylim([-2, 2])
    ylabel('Real part of pressure')  
    legend('show');
end
tasks(1).re_p = real(layer{1}.p);
tasks(1).real_p_Hetmaniuk = real_p_Hetmaniuk;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

layer = setHetmaniukParameters(2);
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
                 'BC', BC, ...
                 'Display', 'none', ...
                 'omega', omega, ...
                 'nu_a', 100, ...
                 'saveRelTermMax', true, ... 
                 'P_inc', -1);

layer{1}.X = layer{1}.R*[0,0,1];
layer{1}.calc_p = true;
[layer,~,~,relTermMaxArr] = e3Dss(layer, options);

real_p_Hetmaniuk = importdata('models/Hetmaniuk2012raa/Figure12.csv');
if plotResults
    figure(12)
    plot(f, real(layer{1}.p),'DisplayName','Exact')
    title('Figure 12 in Hetmaniuk2012raa')
    hold on
    plot(real_p_Hetmaniuk(:,1), real_p_Hetmaniuk(:,2),'DisplayName','Hetmaniuk2012raa')
    xlabel('Frequency [Hz]')
    xlim([f(1), f(end)])
    ylim([-200, 200])
    ylabel('Real part of pressure')  
    legend('show');
    
    figure(42)
    for i = 1:24
        semilogy(0:(size(relTermMaxArr,2)-1),relTermMaxArr(10*i,:).','DisplayName',['\omega = ' num2str(omega(10*i))])
        hold on
    end
    xlabel('n')
    hold off
    legend('show');
end
tasks(2).re_p = real(layer{1}.p);
tasks(2).real_p_Hetmaniuk = real_p_Hetmaniuk;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

layer = setHetmaniukParameters(2);
BC = 'SSBC';

options = struct('applyLoad', 'planeWave', ...
                 'd_vec', d_vec, ...
                 'BC', BC, ...
                 'Display', 'none', ...
                 'omega', omega, ...
                 'saveRelTermMax', true, ... 
                 'P_inc', 1);
             

layer{1}.X = layer{1}.R*[0,0,-1];
layer{1}.calc_p = true;
[layer,~,~,relTermMaxArr] = e3Dss(layer, options);

real_p_Hetmaniuk = importdata('models/Hetmaniuk2012raa/Figure17.csv');
if plotResults
    figure(17)
    plot(f, real(layer{1}.p),'DisplayName','Exact')
    title('Figure 17 in Hetmaniuk2012raa')
    hold on
    plot(real_p_Hetmaniuk(:,1), real_p_Hetmaniuk(:,2),'DisplayName','Hetmaniuk2012raa')
    xlabel('Frequency [Hz]')
    xlim([f(1), f(end)])
    ylim([-15, 15])
    ylabel('Real part of pressure')  
    legend('show');
end
tasks(3).re_p = real(layer{1}.p);
tasks(3).real_p_Hetmaniuk = real_p_Hetmaniuk;

if plotResults
    figure(43)
    for i = 1:24
        semilogy(relTermMaxArr(10*i,:).','DisplayName',['\omega = ' num2str(omega(10*i))])
        hold on
    end
    xlabel('n')
    hold off
    legend('show');
end


