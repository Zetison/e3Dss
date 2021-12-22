% This script produces Figure 19 in 

close all
clear all %#ok

startup
resultsFolder = [folderName '/Benchmarks_timeDomain'];
if ~exist(resultsFolder, 'dir')
    mkdir(resultsFolder);
end

startMatlabPool

intermediatePointCharge = 1; % Places the point charge in between the S1 layer and S5 layer
model = 'S15';
% model = 'S5';
SHBC = 0;
SSBC = 0;
ESBC = 0; 
f_c = 3000;
T = 120/f_c;
N = 2^10;
B = N/T; % bandwidth
f_L = -B/2;
f_R = B/2;
df = 1/T;
f = linspace(0,f_R-df,N/2);
f(end)
% N = 2^2;
% applyLoad = 'planeWave';
applyLoad = 'pointCharge';
% applyLoad = 'surfExcitation';
% applyLoad = 'mechExcitation';
% applyLoad = 'radialPulsation';
switch applyLoad
    case 'pointCharge'
        ESBC = 1;
%         SHBC = 1;
    case 'planeWave'
        ESBC = 1;
    case 'surfExcitation'
        SHBC = 1;
    case 'radialPulsation'  
        SHBC = 1; 
end
switch model
    case 'S5'
        layer = setS5Parameters();
    case 'S15'
        layer = setS15Parameters();
end
R = layer{1}.R;
defineBCstring
R_a = 1.5*R;
P_inc = 1;
theta_s = NaN(1,2);
r_s = 2*R; 
if strcmp(applyLoad,'pointCharge')
    if intermediatePointCharge
        d_vec = -[-sqrt(r_s^2-R_a^2), R_a, 0].';  
        r_s = layer{2}.R*1/3 + layer{3}.R*2/3;
    else
        d_vec = [-sqrt(r_s^2-R_a^2), R_a, 0].'; 
    end
    P_inc = P_inc*r_s;
elseif strcmp(applyLoad,'surfExcitation')
    d_vec = -[-sqrt(r_s^2-R_a^2), R_a, 0].';  
    r_s = layer{1}.R;
    theta_s = [40,60]*pi/180;
elseif strcmp(applyLoad,'mechExcitation')
    d_vec = -[-sqrt(r_s^2-R_a^2), R_a, 0].';  
    r_s = layer{1}.R;
else
    d_vec = [1, 0, 0].';  
end
d_vec = d_vec/norm(d_vec);

options = struct('d_vec', d_vec, ...
                 'P_inc', P_inc, ...
                 'applyLoad', applyLoad, ...
                 'plotTimeOscillation', 0, ...
                 'plotInTimeDomain', 1, ...
                 'BC',BC, ...
                 'SHBC', SHBC, ...
                 'ESBC', ESBC, ...
                 'SSBC', SSBC, ...
                 'R_a', R_a, ...
                 'r_s', r_s, ...
                 'theta_s', theta_s, ...
                 'f_c', f_c, ...
                 'N', N, ...
                 'T', T,...
                 'computeForSolidDomain', strcmp(model,'S5'));

extraPts = 8; %40
folderName = [resultsFolder '/paraviewResults/' model '_' BC '_' applyLoad '/'];
if ~exist(folderName, 'dir')
    mkdir(folderName);
end

createParaviewFiles_e3Dss(extraPts, folderName, layer, options)    


