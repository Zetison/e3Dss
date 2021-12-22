close all
clear all %#ok

layer = setS135Parameters();
SHBC = false;
ESBC = false;
SSBC = false;
f = 30e3;
omega = 2*pi*f; % Angular frequency

defineBCstring
d_vec = [0, 0, -1];

options = struct('d_vec', d_vec,... 
                 'omega', omega, ...
                 'BC', BC, ...
                 'P_inc', 1);
             
layer{1}.X = -layer{1}.R*d_vec;
layer{1}.calc_p_0 = true; % Calculate the far field pattern

%% Create spy matrix (requires to be in debug mode in getCoeffs.m)
% uncomment keyboard command in getCoeffs.m to get spy matrix
data = e3Dss(layer, options);