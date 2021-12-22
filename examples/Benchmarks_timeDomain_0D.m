close all
clear all %#ok

startup

set(0,'defaultTextInterpreter','latex');
startMatlabPool
plotP_inc = 0;
playP_inc = 1;

% applyLoad = 'pointCharge';
applyLoad = 'planeWave';
% applyLoad = 'radialPulsation';

model = 'S35';
% model = 'S5';
f_c = 1500; % (300) source pulse center freq.
ss = 2^5;
% npts = 10;
Fs = 44100;
% Fs = 200000;
dt = 1/Fs;
T = 16; %60/f_c
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
type = 4;
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
    subplot(312), plot(ft,abs(P_incArr))
    subplot(313), plot(ft,atan2(imag(P_incArr),real(P_incArr)));
    drawnow
    %             printResultsToFile([pathToResults 'Pt_inc'], tt.', Pt_inc.', [], 0, 1)
    %             printResultsToFile([pathToResults 'P_inc'], omegat.', abs(P_inc).', [], 0, 1)
    figure(2)
    plot(omegat,abs(P_incArr))
    ylabel('$|P_{\mathrm{inc}}(\omega)|$ [Pa]')
    xlabel('$\omega$ [$\mathrm{s}^{-1}$]')
%     savefig([pathToResults 'Figure17b'])

    figure(3)
    plot(tt,Pt_inc)
    ylabel('$\breve{P}_{\mathrm{inc}}(t)$ [Pa]')
    xlabel('$t$ [s]')
%     savefig([pathToResults 'Figure17a'])
    return
end
if playP_inc
    tt = (0:dt:(T-dt)).';
    Pt_inc = Pt_inc_(tt,0,omega_c,k_c,P_inc,type);
    plot(tt,Pt_inc)
    p_inc = P_inc_(omega(2:end), omega_c,P_inc,type);
    semilogy(f(2:end),abs(p_inc))
    
    sound(Pt_inc,Fs);
%     return
end
%             return
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot dynamic time domain solution in 1D
SSBC = 0;
ESBC = 0; 
SHBC = 0;
switch model
    case 'S15'
        layer = setS15Parameters();
    case 'S35'
        layer = setS35Parameters();
    case 'S5'
        layer = setS5Parameters();
end
layer{1}.calc_p_0 = false;
layer{1}.calc_p = true;
layer{1}.calc_p_inc = false;

defineBCstring

layer{1}.X = [0, 0, 2*layer{1}.R];

if ~exist('options','var')
    options = struct('BC', BC,...
                     'd_vec', d_vec, ...
                     'Display', 'iter', ...
                     'applyLoad',applyLoad);
end
options.omega = omega(2:end);
options.N_max = Inf;
options.d_vec = d_vec;
options.P_inc = @(omega) P_inc_(omega, omega_c,P_inc,type);

if false
    layer{1}.X = [1,0,0];
    layer = e3Dss(layer, options);
    plot(omega(2:end),20*log10(abs(layer{1}.p)))
    plot(omega(2:end),abs(layer{1}.p))
    return
else
    layer = e3Dss(layer, options);
end

startIdx = 2000;
startIdx = 2*round(1900/2);
if N < 100
    startIdx = 1;
end
if type == 4
    startIdx = 1;
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
        totField1(:,n-N/2+1) = layer{1}.p(:,n-N/2);
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

filename = ['../../results/e3Dss/' model '.wav'];
y = real(totFieldTime(end,:));
figure
plot(y)
ys = y/max(abs(y));
if 0
    load handel.mat
    sound(y,Fs);
end
audiowrite(filename,ys,Fs);
% [ys,Fs] = audioread(filename);
sound(Pt_inc,Fs);
sound(ys,Fs);
% sound(y*1000,Fs);


