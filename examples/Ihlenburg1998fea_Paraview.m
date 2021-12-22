close all
clear all %#ok

startup
resultsFolder = [folderName '/Ihlenburg1998fea'];
if ~exist(resultsFolder, 'dir')
    mkdir(resultsFolder);
end

startMatlabPool
applyLoad = 'planeWave';
ESBC = 0;
SSBC = 0;
k_arr = linspace(0.001, 2, 3000).';
% k_arr = [];
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
% k_arr = [1, 2];
for SHBC = 0 %[0, 1]
    for modelCell = {'IL'} %{'IL', 'S5', 'S35', 'S135'}
        model = modelCell{1};
        switch model
            case 'S1'
                layer = setS1Parameters();
            case 'S5'
                layer = setS5Parameters();
            case 'S35'
                layer = setS35Parameters();
            case 'S135'
                layer = setS135Parameters();
            case 'IL'
                layer = setIhlenburgParameters();
        end
        
        defineBCstring
        resultsFolderParaview = [resultsFolder '/paraviewResults/' model '/'];
        if ~exist(resultsFolderParaview, 'dir')
            mkdir(resultsFolderParaview);
        end
%         for i = 1:length(k_arr)
        parfor i = 1:length(k_arr)
            k = k_arr(i);
            omega = k*layer{1}.c_f;
            R_a = 2*layer{1}.R;
            
            alpha_s = pi;
            beta_s = 0;
            beta_f = beta_s;
            beta_f_arr = beta_s;
            d_vec = -[cos(beta_s)*cos(alpha_s);
                      cos(beta_s)*sin(alpha_s);
                      sin(beta_s)]; 
                  
            options = struct('BC', BC, ...
                             'd_vec', d_vec, ...
                             'omega', omega, ...
                             'P_inc', 1, ...
                             'SHBC', SHBC, ...
                             'ESBC', ESBC, ...
                             'SSBC', SSBC, ...
                             'computeForSolidDomain', 0, ...
                             'plotTimeOscillation', 0, ...
                             'plotInTimeDomain', 0, ...
                             'applyLoad', applyLoad, ...
                             'R_a', R_a);

            vtfFileName = [resultsFolderParaview '/' BC '_' num2str(i)];

            extraPts = 15;

            createParaviewFiles_e3Dss(extraPts, vtfFileName, layer, options)
        end
    end
end