close all
clear all %#ok

startup
folderName = [homeDir '/Dropbox/Apps/Overleaf/e3Dss_article2/data'];
resultsFolder = [folderName '/Benchmark_timeDomain_0D'];
if ~exist(resultsFolder, 'dir')
    mkdir(resultsFolder);
end

set(0,'defaultTextInterpreter','latex');
startMatlabPool
plotP_inc = 0;
playP_inc = 1;

% applyLoad = 'pointCharge';
applyLoad = 'planeWave';
% applyLoad = 'radialPulsation';

model = 'Skelton1997tao';
% model = 'S35';
% model = 'S5';
type = 8;
filename = [folderName '/' model '_Type' num2str(type)];
f_c = 1500; % (300) source pulse center freq.
% f_c = 2000*2*pi; % center frequency
ss = 2^5;
% npts = 10;
Fs = 44100;
% Fs = 200000;
dt = 1/Fs;
T = 16; %60/f_c
if type == 6
    f_c = 3300;
end
if type == 7 || type == 8
%     f_c = 337.5;
    f_c = 533.527;
%     f_c = 12720.3;
%     f_c = 14597.8;
%     f_c = 135050;
%     f_c = 137725;
end
% N = 2^11*ss;
N = T*Fs;
B = N/T; % bandwidth
f_L = -B/2;
f_R = B/2;
df = 1/T;
f = linspace(0,f_R-df,N/2);
omega = 2*pi*f;
f(end)
T
dt = T/N;
d_vec = -[0,0,1].';
omega_c = 2*pi*f_c;
c_f = 1500;
k_c = omega_c/c_f;
P_inc = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%
if plotP_inc
    figure(1); clf;
    %             t = linspace(-1/f_c,6/f_c,1000);
    tt = linspace(0,1/f_c,1000);
    tt = [-1e-3,tt,3e-3];
    Pt_inc = Pt_inc_(tt,0,omega_c,k_c,P_inc,type);
    subplot(311), plot(tt,Pt_inc)
    ft = linspace(0,f_R,2000);
    omegat = 2*pi*ft;
    P_incArr = P_inc_(omegat,omega_c,P_inc,type);
    if 0
        subplot(312), plot(ft,abs(P_incArr))
        subplot(313), plot(ft,atan2(imag(P_incArr),real(P_incArr)));
        drawnow
        %             printResultsToFile([pathToResults 'Pt_inc'], tt.', Pt_inc.', [], 0, 1)
        %             printResultsToFile([pathToResults 'P_inc'], omegat.', abs(P_inc).', [], 0, 1)
        figure(2)
        plot(omegat,abs(P_incArr))
        ylabel('$|P_{\mathrm{inc}}(\omega)|$ [Pa]')
        xlabel('$\omega$ [$\mathrm{s}^{-1}$]')
%         savefig([pathToResults 'Figure17b'])    
    end

    figure(3)
    plot(tt,Pt_inc)
    ylabel('$\breve{P}_{\mathrm{inc}}(t)$ [Pa]')
    xlabel('$t$ [s]')
%     savefig([pathToResults 'Figure17a'])
    return
end
tt = (0:dt:(T-dt)).';
if playP_inc
    Pt_inc = Pt_inc_(tt,0,omega_c,k_c,P_inc,type);
    plot(tt,Pt_inc)
    ylabel('$p_{\mathrm{inc}}$ [Pa]','interpreter','latex')
    xlabel('$t$ [s]','interpreter','latex')
    savefig([folderName '/p_inc_Type' num2str(type) '.fig'])
    p_inc = P_inc_(omega(2:end), omega_c,P_inc,type);
    figure
%     semilogy(f(2:end),abs(p_inc))
    semilogy(f(2:end),abs(p_inc))
    
    sound(Pt_inc,Fs);
%     return
end
%             return
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot dynamic time domain solution in 1D
BC = 'NNBC';
switch model
    case 'S15'
        layer = setS15Parameters();
        layer = defineBCstring(layer,BC);
    case 'S35'
        layer = setS35Parameters();
        layer = defineBCstring(layer,BC);
    case 'S5'
        layer = setS5Parameters();
        layer = defineBCstring(layer,BC);
    case 'Skelton1997tao'
        layer = setSkelton1997taoParameters();
        BC = 'SSBC';
end
useFarField = true;
layer{1}.calc_p = ~useFarField;
layer{1}.calc_p_0 = useFarField;
layer{1}.calc_p_inc = false;


layer{1}.X = [0, 0, 2*layer{1}.R];

if ~exist('options','var')
    options = struct('BC', BC,...
                     'd_vec', d_vec, ...
                     'Display', 'iter', ...
                     'applyLoad',applyLoad);
end
options.omega = omega(2:end);
% options.omega = options.omega(46787);
% omega = 2*pi*2924.1875;
% options.omega = omega;
options.N_max = Inf;
options.d_vec = d_vec;
options.P_inc = @(omega) P_inc_(omega, omega_c,P_inc,type);
% options.P_inc = 1;

if false
    layer{1}.X = [1,0,0];
    layer = e3Dss(layer, options);
    plot(omega(2:end),20*log10(abs(layer{1}.p)))
    plot(omega(2:end),abs(layer{1}.p))
    return
else
    [layer,N_eps,flag] = e3Dss(layer, options);
end
% return
startIdx = 2000;
startIdx = 2*round(1900/2);
if N < 100
    startIdx = 1;
end
if type == 4
    startIdx = 1;
end
if type == 5
    startIdx = 1;
end
if type == 6
    startIdx = N-2000;
end
totField1 = zeros(1,N/2);
PincField = zeros(1,N/2);

for n = 0:N-1
    if n >= N/2+1
        if layer{1}.calc_p_inc
            PincField(:,n-N/2+1) = layer{1}.p_inc(:,n-N/2);
        end
% 
%         totField1(:,n-N/2+1) = PincField(:,n-N/2+1) + layer{1}.p(:,n-N/2);
        if useFarField
            totField1(:,n-N/2+1) = layer{1}.p_0(:,n-N/2);
        else
            totField1(:,n-N/2+1) = layer{1}.p(:,n-N/2);
        end
    end
end
totFieldTime = 2/T*fft(totField1,N,2);
PincFieldTime = 2/T*fft(PincField,N,2);

temp = totFieldTime;
totFieldTime(:,1:N-startIdx+1) = temp(:,startIdx:end);
totFieldTime(:,N-startIdx+2:end) = temp(:,1:startIdx-1);

temp = PincFieldTime;
PincFieldTime(:,1:N-startIdx+1) = temp(:,startIdx:end);
PincFieldTime(:,N-startIdx+2:end) = temp(:,1:startIdx-1);

y = real(totFieldTime(end,:));
figure
plot(tt,y)
ys = y/max(abs(y));
if 0
    load handel.mat
    sound(y,Fs);
end
audiowrite([filename '.wav'],ys,Fs);
ylabel('$\mathrm{Re}(p_0)$ [Pa]','interpreter','latex')
xlabel('$t$ [s]','interpreter','latex')
savefig([filename '.fig'])
% [ys,Fs] = audioread(filename);
if playP_inc
    sound(Pt_inc,Fs);
    filename = [folderName '/p_inc_Type' num2str(type) '.wav'];
    audiowrite(filename,Pt_inc,Fs);
end
sound(ys,Fs);
% sound(y*1000,Fs);


