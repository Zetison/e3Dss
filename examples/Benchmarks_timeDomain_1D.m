close all
clear all %#ok

startup

set(0,'defaultTextInterpreter','latex');
startMatlabPool
plotP_inc = 0;
playP_inc = 1;

% applyLoad = 'pointCharge';
% applyLoad = 'planeWave';
applyLoad = 'radialPulsation';

type = 1;
model = 'S15';
% model = 'S5';
f_c = 1500; % (300) source pulse center freq.
if type == 4
    % npts = 10;
    Fs = 44100;
    % Fs = 200000;
    dt = 1/Fs;
    T = 16; %60/f_c
    N = T*Fs;
else
    ss = 1;
    T = 60/f_c;
    N = 2^11*ss;
end
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
    plot(f(2:end),abs(p_inc))
    
%     sound(Pt_inc,Fs);
%     return
end
%             return
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot dynamic time domain solution in 1D
SSBC = 1;
ESBC = 0; 
switch model
    case 'S15'
        layer = setS15Parameters();
        SHBC = 0;
        layer{2}.calc_sigma_s = [1,0,0,0,0,0];
        layer{1}.calc_p = true;
        layer{3}.calc_p = true;
    case 'S5'
        layer = setS5Parameters();
        SHBC = 0;
        layer{2}.calc_sigma_s = [1,0,0,0,0,0];
        layer{1}.calc_p = true;
        layer{3}.calc_p = true;
end

defineBCstring
R_a = 1.5*layer{1}.R;
npts = 1000;
z1 = linspace(layer{1}.R,4*layer{1}.R,npts).';
z2 = linspace(layer{2}.R,layer{1}.R,npts).';
z3 = linspace(layer{3}.R,layer{2}.R,npts).';
layer{1}.X = [zeros(npts,1), zeros(npts,1), z1];
layer{2}.X = [zeros(npts,1), zeros(npts,1), z2];
layer{3}.X = [zeros(npts,1), zeros(npts,1), z3];

if ~exist('options','var')
    options = struct('BC', BC,...
                     'd_vec', d_vec, ...
                     'Display', 'iter', ...
                     'applyLoad',applyLoad);
end
options.omega = omega(2:end);
options.N_max = Inf;
switch applyLoad
    case 'pointCharge'
        switch model
            case 'S15'
                options.r_s = mean([layer{2}.R,layer{3}.R]);
                P_inc = P_inc*options.r_s;
            case 'S5'
                options.r_s = 0;
        end
        options.d_vec = -d_vec;
    case 'planeWave'
        options.d_vec = d_vec;
end
options.P_inc = @(omega) P_inc_(omega, omega_c,P_inc,type);

if false
    layer{1}.X = [1,0,0];
    layer = e3Dss(layer, options);
    plot(omega(2:end),20*log10(abs(layer{1}.p)))
    return
else
    layer = e3Dss(layer, options);
end

startIdx = 2000;
startIdx = 2*round(1900/2);
if N < 100
    startIdx = 1;
end
totField1 = zeros(npts,N/2);
totField2 = zeros(npts,N/2);
totField3 = zeros(npts,N/2);
PincField = zeros(npts,N/2);

for n = 0:N-1
    f = f_L + (f_R-f_L)/N*n;
    omega = 2*pi*f;
    if n >= N/2+1
        if strcmp(applyLoad,'pointCharge')
            r_s = options.r_s;
            k = omega/layer{3}.c_f;
            x_s = r_s*options.d_vec.';
            r = @(y) norm2(repmat(x_s,size(y,1),1)-y);
            PincField(:,n-N/2+1) = options.P_inc(omega)*exp(1i*k*r(layer{3}.X))./r(layer{3}.X);
            
            totField1(:,n-N/2+1) = layer{1}.p(:,n-N/2);
            totField2(:,n-N/2+1) = layer{2}.sigma_rr(:,n-N/2);
            totField3(:,n-N/2+1) = PincField(:,n-N/2+1) + layer{3}.p(:,n-N/2);
        elseif strcmp(applyLoad,'radialPulsation')
            k = omega/layer{1}.c_f;
            PincField(:,n-N/2+1) = options.P_inc(omega)*layer{1}.R.*exp(-1i*k*(norm2(layer{1}.X)-layer{1}.R))./norm2(layer{1}.X);
            
            totField1(:,n-N/2+1) = PincField(:,n-N/2+1) + layer{1}.p(:,n-N/2);
            totField2(:,n-N/2+1) = layer{2}.sigma_rr(:,n-N/2);
            totField3(:,n-N/2+1) = layer{3}.p(:,n-N/2);
        else
            k = omega/layer{1}.c_f;
            k_vec = options.d_vec*k;
            PincField(:,n-N/2+1) = options.P_inc(omega)*exp(1i*dot3(layer{1}.X, k_vec));
            
            totField1(:,n-N/2+1) = PincField(:,n-N/2+1) + layer{1}.p(:,n-N/2);
            totField2(:,n-N/2+1) = layer{2}.p(:,n-N/2);
            totField3(:,n-N/2+1) = layer{3}.p(:,n-N/2);
        end
    end
end
totFieldTime{1} = 2/T*fft(totField1,N,2);
totFieldTime{2} = 2/T*fft(totField2,N,2);
totFieldTime{3} = 2/T*fft(totField3,N,2);
PincFieldTime = 2/T*fft(PincField,N,2);

for i = 1:3
    temp = totFieldTime{i};
    totFieldTime{i}(:,1:N-startIdx+1) = temp(:,startIdx:end);
    totFieldTime{i}(:,N-startIdx+2:end) = temp(:,1:startIdx-1);
end
temp = PincFieldTime;
PincFieldTime(:,1:N-startIdx+1) = temp(:,startIdx:end);
PincFieldTime(:,N-startIdx+2:end) = temp(:,1:startIdx-1);

figure(2)
m_arr = 0:N-1;
for m = m_arr
    t = dt*m;
    plot(z3,real(totFieldTime{3}(:,m+1)),'blue','DisplayName','p_{tot}');
    hold on
    plot(z2,-real(totFieldTime{2}(:,m+1)),'blue','DisplayName','-\sigma_{rr}');
    plot(z1,real(totFieldTime{1}(:,m+1)),'blue','DisplayName','p_{tot}');
    if strcmp(applyLoad,'pointCharge')
        plot(z3,real(PincFieldTime(:,m+1)),'red','DisplayName','p_{inc}');
    elseif strcmp(applyLoad,'radialPulsation') || strcmp(applyLoad,'planeWave')
        plot(z1,real(PincFieldTime(:,m+1)),'red','DisplayName','p_{inc}');
    end
    ylim([-2,2])
    xlim([z3(1),z1(end)])
    plot(layer{1}.R*[1,1],ylim,'--','color','black')
    plot(layer{2}.R*[1,1],ylim,'--','color','black')
    plot(layer{3}.R*[1,1],ylim,'--','color','black')
    hold off
    titleStr = sprintf('Time step %4d, t = %5.3fs. T = %5.3fs.',m,t,T);
    title(titleStr)
%     legend show % This is really slow...
    drawnow
end

filename = '../../results/e3Dss/'+model+'.wav';
y = real(totFieldTime{1}(end,:));
ys = y/max(abs(y));
if 0
    load handel.mat
    sound(y,Fs);
end
audiowrite(filename,ys,Fs);
% [ys,Fs] = audioread(filename);
sound(ys,Fs);
% plot(psd(spectrum.periodogram,ys,'Fs',Fs,'NFFT',length(signal)));


