%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This study is based on Ihlenburg1998fea Figure 5.2
% Ihlenburg1998fea is available at https://books.google.no/books?id=JrMPBwAAQBAJ

close all
clear all %#ok

startup
resultsFolder = [folderName '/Ihlenburg1998fea'];
if ~exist(resultsFolder, 'dir')
    mkdir(resultsFolder);
end

%% Ihlenburg (1998) example
layer = setIhlenburgParameters();
nFreqs = 100000;
k = linspace(2/nFreqs,2,nFreqs)'; % wave number
omega = k*layer{1}.c_f;   % Wave number for outer fluid domain

theta = 180*pi/180;
d_vec = [0,0,1];
SHBC = false;
SSBC = true;
ESBC = false;
defineBCstring
options = struct('BC', BC,...
                 'd_vec', d_vec, ...
                 'omega', omega, ...
                 'P_inc', 1);
layer{1}.calc_p_0 = true; % Calculate the far field pattern
layer{1}.X = -d_vec;
f = @(k)objFunc(k,layer,options);
fminsearchOptions = optimset('TolX',eps, 'TolFun', eps, 'MaxIter', 1000,'MaxFunEvals',1000);
k_extremas = findExtremas(f, 2/nFreqs, 2, 100000,fminsearchOptions)';
k = unique(sort([k; k_extremas]));

R = layer{1}.R;
omega = k*layer{1}.c_f;   % Wave number for outer fluid domain
options.omega = omega;

layer = e3Dss(layer, options);

F = abs(layer{1}.p_0);
plot(k*R, F)
set(0,'defaulttextinterpreter','latex')
hold on
omega = k_extremas*layer{1}.c_f;   % Wave number for outer fluid domain
options.omega = omega;
layer = e3Dss(layer, options);
Fspecial = abs(layer{1}.p_0);
k_resonansFrequencies = k_extremas(Fspecial>10);
F_resonansFrequencies = Fspecial(Fspecial>10);
plot(R*k_resonansFrequencies, F_resonansFrequencies,'x')
xlabel('$$k_1 R_{0,1}$$')
xlim([0, max(k*R)])
ylim([0, 30])
ylabel('$$|F(k)|$$')  
savefig([resultsFolder '/Figure5'])


function F = objFunc(k,layer,options)
options.Display = 'none';
options.omega = k*layer{1}.c_f;
layer = e3Dss(layer, options);
F = abs(layer{1}.p_0);

end
