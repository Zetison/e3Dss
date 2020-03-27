% This script produces Figure 19 in 

close all
clear all %#ok

pathToResults = '../../../../../../hugeFiles/e3Dss/';
% pathToResults = '../../../results/e3Dss/';
% pathToResults = '../results/';

startMatlabPool

intermediatePointCharge = true; % Places the point charge in between the S1 layer and S5 layer
% model = 'S15';
model = 'S5';
SHBC = 0;
SSBC = 0;
ESBC = 0; 
f_c = 1500;
T = 120/f_c;
N = 2^10;
applyLoad = 'planeWave';
% applyLoad = 'pointCharge';
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
R_i = layer{1}.R_i;
defineBCstring
R_a = 1.5*R_i;
P_inc = 1;
r_s = 2*R_i; 
theta_s = NaN(1,2);
d_vec = -[-sqrt(r_s^2-R_a^2), R_a, 0].';  
d_vec = d_vec/norm(d_vec);
if strcmp(applyLoad,'pointCharge')
    if intermediatePointCharge
        r_s = layer{2}.R_i*1/3 + layer{3}.R_i*2/3;
    end
    P_inc = P_inc*r_s;
elseif strcmp(applyLoad,'surfExcitation')
    r_s = layer{1}.R_i;
    theta_s = [40,60]*pi/180;
elseif strcmp(applyLoad,'mechExcitation')
    r_s = layer{1}.R_i;
else
    d_vec = [1, 0, 0].';  
end

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

extraPts = 2;
folderName = [pathToResults 'paraviewResults/' model '_' BC '_' applyLoad '/'];
if ~exist(folderName, 'dir')
    mkdir(folderName);
end

createParaviewFiles_e3Dss(extraPts, folderName, layer, options)    


