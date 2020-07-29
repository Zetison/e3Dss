close all
clear all %#ok

pathToResults = '../../../../../../hugeFiles/e3Dss/';
% pathToResults = '../../../results/e3Dss/';
% pathToResults = '../results';

startMatlabPool
for extraPts = [2,4,8,16,32]
    if false
        alpha_s = 240*pi/180;
        beta_s = 30*pi/180;
        beta_f = beta_s;
        beta_f_arr = beta_s;
        d_vec = -[cos(beta_s)*cos(alpha_s);
                  cos(beta_s)*sin(alpha_s);
                  sin(beta_s)]; 
    else
        d_vec = [1, 0, 0]'; 
    end
    N_max = inf;
    ESBC = 1;
    SSBC = 0;
    for SHBC = 0 %[0, 1]
        for modelCell = {'S5'} %{'IL', 'S5', 'S35', 'S135'}
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
            for f = 1000
                omega = 2*pi*f;
                R_a = 2*layer{1}.R_i;

                options = struct('BC', BC, ...
                                 'd_vec', d_vec, ...
                                 'omega', omega, ...
                                 'P_inc', 1, ...
                                 'SHBC', SHBC, ...
                                 'ESBC', ESBC, ...
                                 'SSBC', SSBC, ...
                                 'N_max', N_max, ...
                                 'computeForSolidDomain', 1, ...
                                 'plotTimeOscillation', 0, ...
                                 'plotInTimeDomain', 0, ...
                                 'applyLoad', 'planeWave', ...
                                 'R_a', R_a);

                folderName = [pathToResults 'nearfields/paraviewResults/' model '/'];
                if ~exist(folderName, 'dir')
                    mkdir(folderName);
                end

                vtfFileName = [folderName '_' BC '_' num2str(extraPts)];
                createParaviewFiles_e3Dss(extraPts, vtfFileName, layer, options)
            end
        end
    end
end

        