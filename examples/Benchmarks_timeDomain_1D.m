close all
clear all %#ok

% pathToResults = '../../../results/e3Dss/';
pathToResults = '../results';
set(0,'defaultTextInterpreter','latex');
startMatlabPool
plotP_inc = 1;

applyLoad = 'pointCharge';
% applyLoad = 'planeWave';
% applyLoad = 'radialPulsation';

% model = 'S15';
model = 'S5';
% f_c = 1500; % (300) source pulse center freq.
f_c = 1500; % (300) source pulse center freq.
N = 2^11;
% N = 2^2;
npts = 1000;
% npts = 10;
T = 120/f_c; %60/f_c
B = N/T; % bandwidth
f_L = -B/2;
f_R = B/2;
df = 1/T;
f = linspace(0,f_R-df,N/2);
omega = 2*pi*f;
type = 1;
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
end
%             return
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot dynamic time domain solution in 1D
SSBC = 0;
ESBC = 0; 
switch model
    case 'S15'
        layer = setS15Parameters();
        SHBC = 1;
        layer{2}.calc_sigma_s = [1,0,0,0,0,0];
        layer{1}.calc_p_0 = false;
        layer{1}.calc_p = true;
        layer{3}.calc_p = true;
    case 'S5'
        layer = setS5Parameters();
        SHBC = 0;
        layer{2}.calc_sigma_s = [1,0,0,0,0,0];
        layer{1}.calc_p_0 = false;
        layer{1}.calc_p = true;
        layer{3}.calc_p = true;
end

defineBCstring
R_a = 1.5*layer{1}.R_i;

z1 = linspace(layer{1}.R_i,4*layer{1}.R_i,npts).';
z2 = linspace(layer{2}.R_i,layer{1}.R_i,npts).';
z3 = linspace(layer{3}.R_i,layer{2}.R_i,npts).';
layer{1}.X = [zeros(npts,1), zeros(npts,1), z1];
layer{2}.X = [zeros(npts,1), zeros(npts,1), z2];
layer{3}.X = [zeros(npts,1), zeros(npts,1), z3];

if ~exist('options','var')
    options = struct('BC', BC,...
                     'd_vec', d_vec, ...
                     'applyLoad',applyLoad);
end
options.omega = omega(2:end);
options.d_vec = d_vec;
options.N_max = Inf;
if strcmp(applyLoad,'pointCharge')
    if false
        options.r_s = mean([layer{2}.R_i,layer{3}.R_i]);
        P_inc = P_inc*options.r_s;
    else
        options.r_s = 0;
    end
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
            x_s = -r_s*options.d_vec.';
            r = @(y) norm2(repmat(x_s,size(y,1),1)-y);
            PincField(:,n-N/2+1) = options.P_inc(omega)*exp(1i*k*r(layer{3}.X))./r(layer{3}.X);
            
            totField1(:,n-N/2+1) = layer{1}.p(:,n-N/2);
            totField2(:,n-N/2+1) = layer{2}.sigma_rr(:,n-N/2);
            totField3(:,n-N/2+1) = PincField(:,n-N/2+1) + layer{3}.p(:,n-N/2);
        elseif strcmp(applyLoad,'radialPulsation')
            k = omega/layer{1}.c_f;
            PincField(:,n-N/2+1) = options.P_inc(omega)*layer{1}.R_i.*exp(-1i*k*(norm2(layer{1}.X)-layer{1}.R_i))./norm2(layer{1}.X);
            
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
dt = T/N;
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
    plot(layer{1}.R_i*[1,1],ylim,'--','color','black')
    plot(layer{2}.R_i*[1,1],ylim,'--','color','black')
    plot(layer{3}.R_i*[1,1],ylim,'--','color','black')
    hold off
    titleStr = sprintf('Time step %4d, t = %5.3fs. T = %5.3fs.',m,t,T);
    title(titleStr)
%     legend show % This is really slow...
    drawnow
end
