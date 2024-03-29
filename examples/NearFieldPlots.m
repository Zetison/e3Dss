%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This study correspond to Figure 16, 18 and 19 Venas2019e3s
% Venas2019e3s is available at https://doi.org/10.1016/j.jsv.2017.08.006
close all
clear all %#ok

startup
resultsFolder = [folderName '/NearFieldPlots'];
if ~exist(resultsFolder, 'dir')
    mkdir(resultsFolder);
end

startMatlabPool

plotTimeOscillation = 0;
% applyLoad = 'planeWave';
% applyLoad = 'pointCharge';
applyLoad = 'surfExcitation';
% applyLoad = 'mechExcitation';
% applyLoad = 'radialPulsation';
extraPts = 8;
N_max = inf;
computeForSolidDomain = 0;
intermediatePointCharge = 0; % Places the point charge in between the S1 layer and S5 layer

f_c = 1500;
T = 120/f_c;
% T = 1;
N = 2^10;
% N = 4;
B = N/T; % bandwidth
f_L = -B/2;
f_R = B/2;
df = 1/T;
f = linspace(0,f_R-df,N/2);

% modelCellArr = {'IL'};
modelCellArr = {'S5', 'S35', 'S135'};
modelCellArr = {'S135'};
modelCellArr = {'S15'};
modelCellArr = {'S45'};
% modelCellArr = {'S5'};
% modelCellArr = {'S1'};
% modelCellArr = {'Skelton1997tao'};
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
            BCarr = {'NNBC'};
        case 'IL'
            BCarr = {'NNBC'};
        case 'Skelton1997tao'
            BCarr = {'SSBC'};
        otherwise
            BCarr = {'NNBC'};
    end
    theta_s = NaN(1,2);
    for BC = BCarr
        switch model
            case 'S1'
                f_arr = 10000;
                layer = setS1Parameters();
                applyLoad = 'planeWave';
                extraPts = 100;
                layer = defineBCstring(layer,BC);
                plotInTimeDomain = false;
                computeForSolidDomain = 1;
                plotTimeOscillation = 1;
            case 'S15' % Figure 19
                f_arr = 1000;
                layer = setS15Parameters();
                extraPts = 2; % 40
                layer = defineBCstring(layer,BC);
                plotInTimeDomain = true;
%                 applyLoad = 'pointCharge';
                applyLoad = 'surfExcitation';
                theta_s = [40,60]*pi/180;
            case 'S5' % Figure 16a,b,c,d
                layer = setS5Parameters();
                layer = defineBCstring(layer,BC);
                extraPts = 16; % 40
                k_arr = 6;
                omega_arr = k_arr*layer{1}.c_f;
                f_arr = omega_arr/(2*pi);
                applyLoad = 'planeWave';
            case 'S35' % Figure 16e,f
                layer = setS35Parameters();
                layer = defineBCstring(layer,BC);
                plotInTimeDomain = false;
                extraPts = 40;
                k_arr = 6;
                omega_arr = k_arr*layer{1}.c_f;
                f_arr = omega_arr/(2*pi);
                applyLoad = 'planeWave';
            case 'S135' % Figure 16g,h
                layer = setS135Parameters();
                layer = defineBCstring(layer,BC);
                plotInTimeDomain = false;
                extraPts = 40;
                k_arr = 6;
                omega_arr = k_arr*layer{1}.c_f;
                f_arr = omega_arr/(2*pi);
                applyLoad = 'planeWave';
            case 'S45'
                layer = setS45Parameters();
                layer = defineBCstring(layer,BC);
                plotInTimeDomain = 1;
                extraPts = 10;
                k_arr = 6;
                omega_arr = k_arr*layer{1}.c_f;
                f_arr = omega_arr/(2*pi);
                applyLoad = 'surfExcitation';
                theta_s = [0,10]*pi/180;
            case 'IL'
                plotTimeOscillation = 0;
                layer = setIhlenburgParameters();
                plotInTimeDomain = false;
                % k_arr = [];
                k_arr = linspace(0.001, 2, 3000).';
                specialValues = [0.250621182794531 %
                               0.320300579445871 %
                               0.370671479527136 %
                               0.412992731010227 %
                               0.454191270410376 %
                               0.499088778976889 %
                               0.551286412239756 %
                               0.613370456080303 %
                               0.687008309336546 %
                               0.773084257718564 %
                               0.871890313908958 %
                               0.983323027396819 %
                               1.107045032953710 %
                               1.242597693362253 %
                               1.389470517759271 %
                               1.547139666101034 %
                               1.715087015757087 %
                               1.892808062465205 %
                               ];
                k_arr = unique(sort([k_arr;specialValues]));
                omega_arr = k_arr*layer{1}.c_f;
                f_arr = omega_arr/(2*pi);
                extraPts = 15;
            case 'Skelton1997tao'
                layer = setSkelton1997taoParameters();
                plotInTimeDomain = false;
                withCoating = false;
                if ~withCoating
                    layerSSBC = layer([1,3]);
                    layerSSBC{1}.R = layer{2}.R;
                    layer = layerSSBC;
                end
                extraPts = 200;
                computeForSolidDomain = 0;
                plotTimeOscillation = 0;
                npts = 10;
                npts = 200;
                f_max = 25e4;
                f_arr = linspace(f_max/npts, f_max, npts).';
                f_arr = 140e3;
        end
        R = 5;
        R_a = 1.5*R;
        P_inc = 1;
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
            if strcmp(model,'S45')
                r_s = layer{end-2}.R;
            else
                r_s = layer{1}.R;
            end
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
                         'plotDisplacementVectors', true, ...
                         'compDisplacementDers', false, ...
                         'applyLoad', applyLoad, ...
                         'r_s', r_s, ...
                         'type', 1, ...
                         'Display','iter', ...
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
            createParaviewFiles_e3Dss(extraPts, vtfFileName, layer, options)
        else
            for f = f_arr
    %         parfor f = f_arr
                omega = 2*pi*f;
                options.omega = omega;
                vtfFileName = [folderName '_' BC{1} '_f_' num2str(f)];
                createParaviewFiles_e3Dss(extraPts, vtfFileName, layer, options)
            end
        end
    end
end

        