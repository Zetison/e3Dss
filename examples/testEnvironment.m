close all
clear all %#ok

addpath ..
addpath ../utils
addpath ../models

%% Test benchmark models
alpha_s = 240*pi/180;
beta_s = 30*pi/180;
alpha_s = 0;
beta_s = 0;

beta_f = beta_s;
beta_f_arr = beta_f;
alpha_f = alpha_s;

d_vec = -[cos(beta_s)*cos(alpha_s);
          cos(beta_s)*sin(alpha_s);
          sin(beta_s)];

% scatteringCase = 'Ray';
scatteringCase = 'BI';
% scatteringCase = 'Sweep';
switch scatteringCase
    case 'Ray'
        alpha_f_arr = alpha_f;
        beta_f_arr = beta_f;
        f = 10e3; % frequencies in Hertz
    case 'BI'
        delta_alpha = 1;
        alpha_f_arr = linspace(0,360,7200)*pi/180;
%         beta_f_arr = 0;
        beta_f_arr = beta_f;
        f = 10e3; % frequencies in Hertz
    case 'Sweep'
        alpha_f_arr = alpha_s;
        delta_f = 1;
        Nfreqs = 3;
        c_f = 1524;
        kstart = 0.01;
        kend = 2;
        f = linspace(kstart*c_f/(2*pi),kend*c_f/(2*pi),Nfreqs);
end
model = 'S1';

switch model
    case 'S1'
        setS1Parameters
    case 'S3'
        setS3Parameters
    case 'S5'
        setS5Parameters
    case 'S13'
        setS13Parameters
    case 'S15'
        setS15Parameters
    case 'S35'
        setS35Parameters
    case 'S135'
        setS135Parameters
    case 'IL'
        setIhlenburgParameters

end
SHBC = true;
ESBC = false;
SSBC = false;
R_i = R_o - t;
omega = 2*pi*f; % Angular frequency

defineBCstring

options = struct('d_vec', d_vec,... 
                 'omega', omega, ...
                 'R_i', R_i, ...
                 'R_o', R_o, ...
                 'P_inc', P_inc, ...
                 'E', E, ...
                 'nu', nu, ...
                 'rho_s', rho_s, ...
                 'rho_f', rho_f, ...
                 'calc_farField', 1, ...
                 'c_f', c_f);
             
v = R_o(1)*[cos(beta_f_arr)*cos(alpha_f_arr); cos(beta_f_arr)*sin(alpha_f_arr); sin(beta_f_arr)*ones(size(alpha_f_arr))]';
data = e3Dss(v, options);
TS = 20*log10(abs(data(1).p));
% plot(180/pi*alpha_f_arr, TS,'DisplayName','Analytic')

% for res = [1,2,3,4]
%     T = readtable(['../../../comsol/models/S1/S1_F' num2str(f) '_res' num2str(res) '_TSVSalpha.txt'],'FileType','text', 'HeaderLines',1);
%     x = T.Var1;
%     y = T.Var2;
% %     plot(x,y,'DisplayName',['resolution = ' num2str(res)])
%     semilogy(x,abs(y-TS),'DisplayName',['resolution = ' num2str(res)])
%     hold on
% end
% 
% legend('off');
% legend('show');
specialValues = [];
% alpha_f_arr = unique(sort([specialValues', linspace(0,2*pi,1000)]));
k = omega./c_f;
p_inc = @(v) P_inc*exp(1i*k*dot3(v,d_vec));
plot(180/pi*alpha_f_arr, real(data(1).p+p_inc(v)), 180/pi*alpha_f_arr, 2*real(p_inc(v)))
xlim([60,240])
figure(2)
plot(180/pi*alpha_f_arr, imag(data(1).p+p_inc(v)), 180/pi*alpha_f_arr, 2*imag(p_inc(v)))
xlim([60,240])
% figure(3)
% plot(180/pi*alpha_f_arr, abs(data(1).p+p_inc(v)))
% xlim([60,240])
% 
% plot(180/pi*alpha_f_arr, real(data(1).p))
% xlim([0,360])
% ylim([-11,-1.8])
