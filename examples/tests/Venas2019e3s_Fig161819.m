function tasks = Venas2019e3s_Fig161819(plotResults)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This study correspond to Figure 16, 18 and 19 Venas2019e3s
% Venas2019e3s is available at https://doi.org/10.1016/j.jsv.2017.08.006
if nargin < 1
    plotResults = false;
end

startup
resultsFolder = [folderName '/NearFieldPlots'];
if ~exist(resultsFolder, 'dir')
    mkdir(resultsFolder);
end

startMatlabPool

plotTimeOscillation = 0;
% applyLoad = 'planeWave';
applyLoad = 'pointCharge';
% applyLoad = 'surfExcitation';
% applyLoad = 'mechExcitation';
% applyLoad = 'radialPulsation';
extraPts = 1;
N_max = inf;
computeForSolidDomain = 0;
intermediatePointCharge = 0; % Places the point charge in between the S1 layer and S5 layer

f_c = 1500;
T = 120/f_c;
% T = 1;
N = 2^10;
N = 2^4;
% N = 4;
B = N/T; % bandwidth
f_L = -B/2;
f_R = B/2;
df = 1/T;
f = linspace(0,f_R-df,N/2);

% modelCellArr = {'IL'};
modelCellArr = {'S15', 'S5', 'S35', 'S135'};
tasks = {};
for modelCell = modelCellArr
    model = modelCell{1};
    switch model
        case 'S1'
            BCarr = {'ESBC'};
        case 'S5' % Figure 16a,b,c,d or 18
            plotInTimeDomain = true;
            if plotInTimeDomain
                BCarr = {'ESBC'};
                computeForSolidDomain = 1;
            else
                BCarr = {'SHBC','NNBC'};
                computeForSolidDomain = 0;
            end
        case 'S35' % Figure 16e,f
            BCarr = {'SHBC'};
        case 'S135' % Figure 16g,h
            BCarr = {'NNBC'};
        case 'S15'
            BCarr = {'ESBC'};
        case 'IL'
            BCarr = {'NNBC'};
        case 'Skelton1997tao'
            BCarr = {'SSBC'};
    end
    for BC = BCarr
        switch model
            case 'S1'
                f_arr = 1000;
                layer = setS1Parameters();
                extraPts = 1;
                layer = defineBCstring(layer,BC);
                plotInTimeDomain = true;
                computeForSolidDomain = 1;
            case 'S15' % Figure 19
                f_arr = 1000;
                layer = setS15Parameters();
                extraPts = 1; % 40
                layer = defineBCstring(layer,BC);
                plotInTimeDomain = true;
                applyLoad = 'pointCharge';
            case 'S5' % Figure 16a,b,c,d
                layer = setS5Parameters();
                layer = defineBCstring(layer,BC);
                extraPts = 1; % 40
                k_arr = 6;
                omega_arr = k_arr*layer{1}.c;
                f_arr = omega_arr/(2*pi);
                applyLoad = 'planeWave';
            case 'S35' % Figure 16e,f
                layer = setS35Parameters();
                layer = defineBCstring(layer,BC);
                plotInTimeDomain = false;
                extraPts = 1;
                k_arr = 6;
                omega_arr = k_arr*layer{1}.c;
                f_arr = omega_arr/(2*pi);
                applyLoad = 'planeWave';
            case 'S135' % Figure 16g,h
                layer = setS135Parameters();
                layer = defineBCstring(layer,BC);
                plotInTimeDomain = false;
                extraPts = 1;
                k_arr = 6;
                omega_arr = k_arr*layer{1}.c;
                f_arr = omega_arr/(2*pi);
                applyLoad = 'planeWave';
        end
        R = layer{1}.R;
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
%         R_a = 2*layer{1}.R;

        options = struct('BC', BC{1}, ...
                         'd_vec', d_vec, ...
                         'P_inc', P_inc, ...
                         'N_max', N_max, ...
                         'computeForSolidDomain', computeForSolidDomain, ...
                         'plotTimeOscillation', plotTimeOscillation, ...
                         'plotInTimeDomain', plotInTimeDomain, ...
                         'plotDisplacementVectors', false, ...
                         'compDisplacementDers', false, ...
                         'applyLoad', applyLoad, ...
                         'r_s', r_s, ...
                         'Display','none', ...
                         'theta_s', theta_s, ...
                         'f_c', f_c, ...
                         'N', N, ...
                         'T', T,...
                         'R_a', R_a);

        folderName = [resultsFolder '/paraviewResults/' model '/'];
        if ~exist(folderName, 'dir')
            mkdir(folderName);
        end
        if plotInTimeDomain
            omega = 2*pi*f;
            options.omega = omega;
            vtfFileName = [folderName '_' BC{1}];
            tasks_i = createParaviewFiles_e3Dss(extraPts, vtfFileName, layer, options);
            tasks = [tasks, tasks_i];
        else
            for f = f_arr
    %         parfor f = f_arr
                omega = 2*pi*f;
                options.omega = omega;
                vtfFileName = [folderName '_' BC{1} '_f_' num2str(f)];
                tasks_i = createParaviewFiles_e3Dss(extraPts, vtfFileName, layer, options);
                tasks = [tasks, tasks_i];
            end
        end
    end
end

        