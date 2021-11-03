close all
clear all %#ok
% 
startup
resultsFolder = [folderName '/Benchmarks_NearFieldPlots_1D'];
if ~exist(resultsFolder, 'dir')
    mkdir(resultsFolder);
end

%% Test benchmark models

beta_f = pi/2;
alpha_f = 0;
npts = 1000;
P_inc = 1;

f = 1e3; % frequencies in Hertz
% model = 'S15';
model = 'S1';
% applyLoad = 'planeWave';
applyLoad = 'radialPulsation';
switch model
    case 'S1'
        layer = setS1Parameters();
    case 'S15'
        layer = setS15Parameters();
end
omega = 2*pi*f; % Angular frequency
k = omega./layer{1}.c_f;
switch applyLoad
    case 'planeWave'
        d_vec = [0;0;1];
        p_inc = @(r) P_inc*exp(1i*k*r);
    case 'radialPulsation'
        p_inc = @(r) P_inc*layer{1}.R_i*exp(-1i*k*(r-layer{1}.R_i))./r;
end

r_arr1 = linspace(layer{1}.R_i,2*layer{1}.R_i,npts).';
r_arr2 = linspace(layer{2}.R_i,layer{1}.R_i,npts).';
r_arr3 = linspace(layer{3}.R_i,layer{2}.R_i,npts).';
X{1} = r_arr1*[cos(beta_f)*cos(alpha_f), cos(beta_f)*sin(alpha_f), sin(beta_f)*ones(size(alpha_f))];
X{2} = r_arr2*[cos(beta_f)*cos(alpha_f), cos(beta_f)*sin(alpha_f), sin(beta_f)*ones(size(alpha_f))];
X{3} = r_arr3*[cos(beta_f)*cos(alpha_f), cos(beta_f)*sin(alpha_f), sin(beta_f)*ones(size(alpha_f))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SHBC
switch model
    case 'S1'
        layer = setS1Parameters();
    case 'S15'
        layer = setS15Parameters();
end
layer{1}.X = X{1};
layer{2}.X = X{2};
layer{3}.X = X{3};
layer{2}.calc_sigma_s = [1,0,0,0,0,0];
layer{1}.calc_p = true;
layer{3}.calc_p = true;
SHBC = true;
ESBC = false;
SSBC = false;

defineBCstring

options = struct('BC', BC, ...
                 'omega', omega, ...
                 'P_inc', P_inc, ...
                 'applyLoad',applyLoad);
% data = e3Dss(layer{1}.X, options);
layer = e3Dss(layer, options);
figure(42)
switch numel(layer)
    case 1
        plot(r_arr1, real(p_inc(r_arr1) + layer{1}.p))
    case 2
        plot([r_arr2; r_arr1], [-real(layer{2}.sigma_rr); real(p_inc(r_arr1) + layer{1}.p)])
    otherwise
        plot([r_arr3; r_arr2; r_arr1], [real(layer{3}.p); -real(layer{2}.sigma_rr); real(p_inc(r_arr1) + layer{1}.p)])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ESBC
switch model
    case 'S1'
        layer = setS1Parameters();
    case 'S15'
        layer = setS15Parameters();
end
layer{1}.X = X{1};
layer{2}.X = X{2};
layer{3}.X = X{3};
layer{2}.calc_sigma_s = [1,0,0,0,0,0];
layer{1}.calc_p = true;
layer{3}.calc_p = true;
SHBC = false;
ESBC = false;
SSBC = true;

defineBCstring

options.BC = BC;

layer = e3Dss(layer, options);
hold on
switch numel(layer)
    case 1
        plot(r_arr1, real(p_inc(r_arr1) + layer{1}.p))
    case 2
        plot([r_arr2; r_arr1], [-real(layer{2}.sigma_rr); real(p_inc(r_arr1) + layer{1}.p)])
    otherwise
        plot([r_arr3; r_arr2; r_arr1], [real(layer{3}.p); -real(layer{2}.sigma_rr); real(p_inc(r_arr1) + layer{1}.p)])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NNBC
switch model
    case 'S1'
        layer = setS1Parameters();
    case 'S15'
        layer = setS15Parameters();
end
layer{1}.X = X{1};
layer{2}.X = X{2};
layer{3}.X = X{3};
layer{2}.calc_sigma_s = [1,0,0,0,0,0];
layer{1}.calc_p = true;
layer{3}.calc_p = true;
SHBC = false;
ESBC = false;
SSBC = false;

defineBCstring

options.BC = BC;

layer = e3Dss(layer, options);
hold on
switch numel(layer)
    case 1
        plot(r_arr1, real(p_inc(r_arr1) + layer{1}.p))
    case 2
        plot([r_arr2; r_arr1], [-real(layer{2}.sigma_rr); real(p_inc(r_arr1) + layer{1}.p)])
    otherwise
        plot([r_arr3; r_arr2; r_arr1], [real(layer{3}.p); -real(layer{2}.sigma_rr); real(p_inc(r_arr1) + layer{1}.p)])
end
legend('SHBC','SSBC','NNBC')
legend show
hold off

