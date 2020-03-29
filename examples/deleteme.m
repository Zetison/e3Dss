%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This study correspond to Fig. 1 and Fig. 4 in Sage1979mri
% Sage1979mri is available at https://doi.org/10.1121/1.382928

close all
clear all %#ok

%% Define parameters
a = 1;
z = 0;
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
% omega = 2*pi*f;
type = 1;
d_vec = -[0,0,1].';
omega_c = 2*pi*f_c;
c_f = 1500;
k_c = omega_c/c_f;
P_inc = 1;
t = linspace(0,1/f_c,1000);
t = [-1e-3,t,3e-3];
P = Pt_inc_(t,z,omega_c,k_c,P_inc,1);
% plot(t,P)
% hold on
% for terms = 1:4
%     P = Pt_inc_(t,z,omega_c,k_c,P_inc,3,terms);
% 
%     plot(t,P)
%     hold on
% end


ft = linspace(0,f_R,2000);
omegat = sort([2*pi*ft,omega_c]);
P = P_inc_(omegat,omega_c,P_inc,1);
plot(omegat,abs(P))
P = P_inc_(omegat,omega_c,P_inc,3);
hold on
plot(omegat,abs(P))