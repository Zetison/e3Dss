function tasks = Ayres1987ars_Fig1(plotResults)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This study correspond to Figure 1 in Ayres1987ars
% Ayres1987ars is available at https://doi.org/10.1121/1.394950

if nargin < 1
    plotResults = false;
end

startup
folderName = [homeDir '/Dropbox/Apps/Overleaf/createFigures/data/e3Dss_article2'];
resultsFolder = [folderName '/Ayres1987ars'];
if ~exist(resultsFolder, 'dir')
    mkdir(resultsFolder);
end

alpha = [0,0.5e7,0.5e8];
beta = [0,0.1e6,0.1e7];
% alpha = [0,5e7,5e7];
% beta = [0,0.1e6,1e6];
fluids = {'air'};
for j = 1:numel(fluids)
    fluid = fluids{j};
    if plotResults
        figure
        title(sprintf('Figure 1-3 in Ayres1987ars (%s)', fluid))
    end
    for i = 1:numel(alpha)
        layer = setAyres1987arsParameters(fluid,alpha(i),beta(i));
        a = layer{1}.R;
        k1a_max = 20;
        npts = 3000;
        k1a = linspace(k1a_max/npts,k1a_max,npts).';
        k = k1a/a;
        omega = k*layer{1}.c_f;
        d_vec = [0,0,1].';
        options = struct('applyLoad', 'planeWave', ...
                         'Display', 'none', ...
                         'd_vec', d_vec, ...
                         'BC', 'NNBC', ...
                         'N_max', 50-1, ...
                         'omega', omega);
        layer{1}.X = layer{1}.R*[0,0,-1];
        layer{1}.calc_p_0 = true;

        if 0
            startMatlabPool
            f = @(k1a)-objFunc(k1a,layer,options);
            specialValues = findExtremas(f, k1a(1), k1a(end), 1e7)';
            delta = 1e-5*k1a(end);
            specialValues = sort([specialValues; (specialValues-delta); (specialValues+delta)]);
            save(['miscellaneous/Ayres_' fluid '_extremas_' num2str(i)], 'specialValues')
        else
            load(['miscellaneous/Ayres_' fluid '_extremas_' num2str(i)])
        end
        k1a = unique(sort([k1a; specialValues]));
        options.omega = k1a/layer{1}.R*layer{1}.c_f;
        layer = e3Dss(layer, options);

        if plotResults
            plot(k1a, abs(layer{1}.p_0)*2/a,'DisplayName',sprintf('e3Dss: Rubber sphere $\\alpha$ = %.1e, $\\beta$ = %.1e',alpha(i),beta(i)))
            hold on
            printResultsToFile([resultsFolder '/Figure1_' fluid '_NNBC_' num2str(i)], {'x', k1a, 'y', abs(layer{1}.p_0).'*2/a, 'xlabel','k1a', 'ylabel','fs'})
            xlabel('$k_1 a$')
        end
        tasks(i).p_0_2a = abs(layer{1}.p_0)*2/a;
    end
    if strcmp(fluid,'air')
        form_Ayres = importdata('models/Ayres1987ars/Figure1a.csv');
        if plotResults
            plot(form_Ayres(:,1), form_Ayres(:,2),'DisplayName','Ayres1987ars')
        end
        options.BC = 'SHBC';
        layerSHBC = e3Dss(layer(1), options);
        if plotResults
            plot(k1a, abs(layerSHBC{1}.p_0)*2/a,'DisplayName','SHBC')
            printResultsToFile([resultsFolder '/Figure1_SHBC'], {'x', k1a, 'y', abs(layerSHBC{1}.p_0).'*2/a, 'xlabel','k1a', 'ylabel','fs'})
            printResultsToFile([resultsFolder '/Figure1_Ayres'], {'x', form_Ayres(:,1), 'y', form_Ayres(:,2), 'xlabel','k1a', 'ylabel','fs'})
            hleg = legend('show','interpreter','latex');
            plot([k1a(1),k1a(end)],[1,1],'black','DisplayName','')
            hleg.String(end) = [];
        end
        tasks(4).p_0_2a = abs(layerSHBC{1}.p_0)*2/a;
    end
    if plotResults
        ylim([0,1.2])
    end
end

function fs = objFunc(k1a,layer,options)
options.Display = 'none';
options.omega = k1a/layer{1}.R*layer{1}.c_f;
layer = e3Dss(layer, options);
fs = abs(layer{1}.p_0)*2/layer{1}.R;


