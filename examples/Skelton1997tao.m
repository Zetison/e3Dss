%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This study correspond to Figure 10.4, Figure 10.5, Figure 10.6 and Figure 10.7 in Skelton1997tao
% Skelton1997tao is available at https://doi.org/10.1142/9781848160750_0001

close all
clear all %#ok

startup
folderName = [homeDir '/Dropbox/Apps/Overleaf/createFigures/data/e3Dss_article2'];
resultsFolder = [folderName '/Skelton1997tao'];
if ~exist(resultsFolder, 'dir')
    mkdir(resultsFolder);
end

withDebug = true;

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

options = struct('BC', 'SSBC', ...
                 'd_vec', d_vec, ...
                 'debug', withDebug, ...
                 'Display', 'final', ...
                 'omega', omega);

layer{1}.X = -options.d_vec.'; % Compute backscattered pressure

layerSSBC = layer([1,3]);
layerSSBC{1}.R_i = layer{2}.R_i;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run simulation
figure(1)
layerCoating = e3Dss(layer, options);
TS = 20*log10(abs(layerCoating{1}.p_0));
plot(f,TS,'DisplayName','e3Dss with coating')
hold on
printResultsToFile([resultsFolder '/Figure10.45_e3Dss_coating'], {'x', f.', 'y', TS.', 'xlabel','f', 'ylabel','TS'})
% return
layerSSBC = e3Dss(layerSSBC, options);
TS = 20*log10(abs(layerSSBC{1}.p_0));
plot(f,TS,'DisplayName','e3Dss without coating')
printResultsToFile([resultsFolder '/Figure10.45_e3Dss_SSBC'], {'x', f.', 'y', TS.', 'xlabel','f', 'ylabel','TS'})

options.BC = 'SHBC';
layerSHBC = e3Dss(layerSSBC(1), options);
TS = 20*log10(abs(layerSHBC{1}.p_0));
printResultsToFile([resultsFolder '/Figure10.45_e3Dss_SHBC'], {'x', f.', 'y', TS.', 'xlabel','f', 'ylabel','TS'})

plot(f,TS,'DisplayName','e3Dss SHBC')
set(0,'defaulttextinterpreter','latex')
xlabel('Frequency (Hz)')
ylabel('Target strength')
xlim([0 f_max])
ylim([-30,20])
title('Figure 10.4 and Figure 10.5')

TS_Coating = importdata('models/Skelton1997tao/Figure10.5.csv');
plot(TS_Coating(:,1),TS_Coating(:,2),'DisplayName','Ref Coating')
printResultsToFile([resultsFolder '/Figure10.45_coating'], {'x', TS_Coating(:,1), 'y', TS_Coating(:,2), 'xlabel','f', 'ylabel','TS'})
TS_SSBC = importdata('models/Skelton1997tao/Figure10.4_SSBC.csv');
plot(TS_SSBC(:,1),TS_SSBC(:,2),'DisplayName','Ref SSBC')
printResultsToFile([resultsFolder '/Figure10.45_SSBC'], {'x', TS_SSBC(:,1), 'y', TS_SSBC(:,2), 'xlabel','f', 'ylabel','TS'})
TS_SHBC = importdata('models/Skelton1997tao/Figure10.4_SHBC.csv');
plot(TS_SHBC(:,1),TS_SHBC(:,2),'DisplayName','Ref SHBC')
printResultsToFile([resultsFolder '/Figure10.45_SHBC'], {'x', TS_SHBC(:,1), 'y', TS_SHBC(:,2), 'xlabel','f', 'ylabel','TS'})
legend show
savefig([resultsFolder '/figure1.fig'])

figure(2)
% npts = 2;
npts = 2000;
% npts = 100;
f_max = 25e3;
f = linspace(f_max/npts,f_max,npts);
omega = 2*pi*f;
options = struct('BC', 'SSBC', ...
                 'd_vec', d_vec, ...
                 'debug', withDebug, ...
                 'omega', omega);

layerSSBC = e3Dss(layerSSBC, options);
TS = 20*log10(abs(layerSSBC{1}.p_0));
plot(f,TS,'DisplayName','e3Dss')
TS_SSBC1 = importdata('models/Skelton1997tao/Figure10.4_SSBC.csv');
TS_SSBC2 = importdata('models/Skelton1997tao/Figure10.6.csv');
TS_SSBC = [TS_SSBC1(TS_SSBC1(:,1) < 1000,1), TS_SSBC1(TS_SSBC1(:,1) < 1000,2);
            TS_SSBC2(TS_SSBC2(:,1) > 1000,1), TS_SSBC2(TS_SSBC2(:,1) > 1000,2)];
hold on
plot(TS_SSBC(:,1),TS_SSBC(:,2),'DisplayName','Ref')
ylim([-30,20])
title('Figure 10.6')
xlabel('Frequency (Hz)')
ylabel('Target strength')
legend show
savefig([resultsFolder '/figure2.fig'])
printResultsToFile([resultsFolder '/Figure10.6'], {'x', TS_SSBC(:,1), 'y', TS_SSBC(:,2), 'xlabel','f', 'ylabel','TS'})
printResultsToFile([resultsFolder '/Figure10.6_e3Dss_1'], {'x', f.', 'y', TS.', 'xlabel','f', 'ylabel','TS'})

% return

figure(22)
f_max = 25e4;
% f_max = 16e4;
% npts = 2*npts;
f = [f, linspace(f(end)+(f_max-f(end))/npts,f_max,npts)];
omega = 2*pi*f;
k = omega./layer{1}.c_f;
options = struct('BC', 'SSBC', ...
                 'd_vec', d_vec, ...
                 'debug', withDebug, ...
                 'omega', omega);
layerSSBC1 = e3Dss(layerSSBC, options);
TS = 20*log10(abs(layerSSBC1{1}.p_0));
plot(f,TS,'DisplayName','With scaling')
ylim([-30,20])
xlabel('Frequency (Hz)')
ylabel('Target strength')
printResultsToFile([resultsFolder '/Figure10.6_e3Dss_2'], {'x', f.', 'y', TS.', 'xlabel','f', 'ylabel','TS'})
options.nu_a = -1;
layerSSBC2 = e3Dss(layerSSBC, options);
TS = 20*log10(abs(layerSSBC2{1}.p_0));
hold on
plot(TS_SSBC(:,1),TS_SSBC(:,2),'DisplayName','Ref')
plot(f,TS,'DisplayName','Without scaling')
legend show
printResultsToFile([resultsFolder '/Figure10.6_e3Dss_3'], {'x', f.', 'y', TS.', 'xlabel','f', 'ylabel','TS'})
savefig([resultsFolder '/figure2b.fig'])

if false
    % this case is not perfectly reproducable due to lacking parameters
    figure(3)
    f_c = 15e3; % pulse center freq.
    N = 2^12;       % **NOTE** it is not clear what Skelton uses for this parameter ("The time scale, in milliseconds, has an arbitrary origin.")
    T = 10*5e-3;    % **NOTE** it is not clear what Skelton uses for this parameter ("The time scale, in milliseconds, has an arbitrary origin.")
    B = N/T; % bandwidth
    f_L = -B/2;
    f_R = B/2;
    df = 1/T;
    f = linspace(0,f_R-df,N/2);
    omega = 2*pi*f;
    omega_c = 2*pi*f_c;
    P_inc = 1;
    layerSSBC{1}.X = -d_vec.'; % Compute backscattered pressure
    plotP_inc = true;
    type = 3; % The time variation of the plane wave is a f_c=15kHz sine wave on for 1 cycle

    options = struct('BC', 'SSBC', ...
                     'd_vec', d_vec, ...
                     'debug', withDebug, ...
                     'omega', omega(2:end));
    options.P_inc = @(omega) P_inc_(omega, omega_c,P_inc,type);

    if plotP_inc
        k_c = omega_c/layer{1}.c_f;
        %             t = linspace(-1/f_c,6/f_c,1000);
        tt = linspace(0,1/f_c,1000);
        tt = [-1e-3,tt,3e-3];
        Pt_inc = Pt_inc_(tt,0,omega_c,k_c,P_inc,type);
        plot(tt,Pt_inc)
        ft = linspace(0,f_R,2000);
        omegat = 2*pi*ft;
        P_incArr = P_inc_(omegat,omega_c,P_inc,type);

        figure(5)
        plot(omegat,abs(P_incArr))
        ylabel('$|P_{\mathrm{inc}}(\omega)|$ [Pa]')
        xlabel('$\omega$ [$\mathrm{s}^{-1}$]')

        figure(6)
        plot(tt,Pt_inc)
        ylabel('$\breve{P}_{\mathrm{inc}}(t)$ [Pa]')
        xlabel('$t$ [s]')
    end

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
    p_ref = importdata('models/Skelton1997tao/Figure10.7.csv');
    hold on
    plot(1000*p_ref(:,1),1000*p_ref(:,2),'DisplayName','Ref')
    ylim([-600,600])
    xlim([0,5])
    title('Figure 10.7')
    xlabel('Time (ms)')
    ylabel('Pressure (mPa)')
    legend show
    savefig([resultsFolder '/figure3.fig'])

end

