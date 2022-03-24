function tasks = Skelton1997tao_Fig10567(plotResults)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This study correspond to Figure 10.4, Figure 10.5, Figure 10.6 and Figure 10.7 in Skelton1997tao
% Skelton1997tao is available at https://doi.org/10.1142/9781848160750_0001

if nargin < 1
    plotResults = false;
end
startup
folderName = [homeDir '/Dropbox/Apps/Overleaf/createFigures/data/e3Dss_article2'];
resultsFolder = [folderName '/Skelton1997tao'];
if ~exist(resultsFolder, 'dir')
    mkdir(resultsFolder);
end


%% Define parameters
theta = 180*pi/180;
theta_s = theta;
layer = setSkelton1997taoParameters();
layer{1}.calc_p_0 = true; % Calculate the far field pattern

%% Calculate dependent parameters
npts = 2000;
f_max = 1000;
f = linspace(f_max/npts,f_max,npts);
omega = 2*pi*f;
d_vec = [0,0,1].';

withDebug = 0;
options = struct('BC', 'SSBC', ...
                 'd_vec', d_vec, ...
                 'debug', withDebug, ...
                 'Display', 'none', ...
                 'omega', omega);

layer{1}.X = -options.d_vec.'; % Compute backscattered pressure

layerSSBC = layer([1,3]);
layerSSBC{1}.R = layer{2}.R;

load('miscellaneous/Skelton_extremas')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run simulation
f = sort(unique([f, specialValues(specialValues <= f_max).']));
options.omega = 2*pi*f;
layerCoating = e3Dss(layer, options);
TS = 20*log10(abs(layerCoating{1}.p_0));
if plotResults
    figure(1)
    plot(f,TS,'DisplayName','e3Dss with coating')
    hold on
end
tasks(1).TS_coating_e3Dss = TS;
% return
layerSSBC = e3Dss(layerSSBC, options);
TS = 20*log10(abs(layerSSBC{1}.p_0));
if plotResults
    plot(f,TS,'DisplayName','e3Dss without coating')
end
tasks(1).TS_SSBC_e3Dss = TS;

options.BC = 'SHBC';
layerSHBC = e3Dss(layerSSBC(1), options);
TS = 20*log10(abs(layerSHBC{1}.p_0));
tasks(1).TS_SHBC_e3Dss = TS;

if plotResults
    plot(f,TS,'DisplayName','e3Dss SHBC')
    set(0,'defaulttextinterpreter','latex')
    xlabel('Frequency (Hz)')
    ylabel('Target strength')
    xlim([0 f_max])
    ylim([-30,20])
    title('Figure 10.4 and Figure 10.5')
end

TS_Coating = importdata('models/Skelton1997tao/Figure10.5.csv');
TS_SSBC = importdata('models/Skelton1997tao/Figure10.4_SSBC.csv');
TS_SHBC = importdata('models/Skelton1997tao/Figure10.4_SHBC.csv');
if plotResults
    plot(TS_Coating(:,1),TS_Coating(:,2),'DisplayName','Ref Coating')
    plot(TS_SSBC(:,1),TS_SSBC(:,2),'DisplayName','Ref SSBC')
    plot(TS_SHBC(:,1),TS_SHBC(:,2),'DisplayName','Ref SHBC')
    legend show
end
tasks(1).TS_Coating = TS_Coating;
tasks(1).TS_SSBC = TS_SSBC;
tasks(1).TS_SHBC = TS_SHBC;

% npts = 2;
npts = 2000;
% npts = 100;
f_max = 25e3;
f = linspace(f_max/npts,f_max,npts);
f = sort(unique([f, specialValues(specialValues <= f_max).']));
omega = 2*pi*f;
options = struct('BC', 'SSBC', ...
                 'd_vec', d_vec, ...
                 'debug', withDebug, ...
                 'Display', 'none', ....
                 'omega', omega);

layerSSBC = e3Dss(layerSSBC, options);
TS = 20*log10(abs(layerSSBC{1}.p_0));
TS_SSBC1 = importdata('models/Skelton1997tao/Figure10.4_SSBC.csv');
TS_SSBC2 = importdata('models/Skelton1997tao/Figure10.6.csv');
TS_SSBC = [TS_SSBC1(TS_SSBC1(:,1) < 1000,1), TS_SSBC1(TS_SSBC1(:,1) < 1000,2);
            TS_SSBC2(TS_SSBC2(:,1) > 1000,1), TS_SSBC2(TS_SSBC2(:,1) > 1000,2)];
if plotResults
    figure(2)
    plot(f,TS,'DisplayName','e3Dss')
    hold on
    plot(TS_SSBC(:,1),TS_SSBC(:,2),'DisplayName','Ref')
    ylim([-30,20])
    title('Figure 10.6')
    xlabel('Frequency (Hz)')
    ylabel('Target strength')
    legend show
end
tasks(2).TS = TS;
tasks(2).TS_ref = TS_SSBC;

% return

f_max = 25e4;
% f_max = 16e4;
% npts = 2*npts;
f = [f, linspace(f(end)+(f_max-f(end))/npts,f_max,npts)];
f = sort(unique([f, specialValues(specialValues <= f_max).']));
omega = 2*pi*f;
k = omega./layer{1}.c;
options = struct('BC', 'SSBC', ...
                 'd_vec', d_vec, ...
                 'debug', withDebug, ...
                 'Display', 'none', ....
                 'omega', omega);
layerSSBC1 = e3Dss(layerSSBC, options);
TS = 20*log10(abs(layerSSBC1{1}.p_0));
if plotResults
    figure(22)
    plot(f,TS,'DisplayName','With scaling')
    ylim([-30,20])
    xlabel('Frequency (Hz)')
    ylabel('Target strength')
end
options.nu_a = -1;
[layerSSBC2,~,flag] = e3Dss(layerSSBC, options);
TS2 = 20*log10(abs(layerSSBC2{1}.p_0));

tasks(3).TS = TS;
tasks(3).TS2 = TS2;
if plotResults
    hold on
    plot(TS_SSBC(:,1),TS_SSBC(:,2),'DisplayName','Ref')
    plot(f(~flag),TS2(~flag),'DisplayName','Without scaling')
    legend show
end
