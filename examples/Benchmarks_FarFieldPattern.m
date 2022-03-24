close all
clear all %#ok


startup
resultsFolder = [folderName '/Benchmark_FarFieldPattern'];
if ~exist(resultsFolder, 'dir')
    mkdir(resultsFolder);
end

% Calculate Far-Field pattern (BI and Sweep)
tic
alpha_s = 240*pi/180;
beta_s = 30*pi/180;
beta_f = beta_s;
beta_f_arr = beta_s;
d_vec = -[cos(beta_s)*cos(alpha_s);
          cos(beta_s)*sin(alpha_s);
          sin(beta_s)];
models = {'S1','S3','S5','S13','S15','S35','S135'};
% models = {'tripleShell'};
models = {'S135'};
BCarr = {'SHBC','SSBC','ESBC','NNBC'};
% BCarr = {'NNBC'};

for scatteringCase = {'Sweep'} %,'Sweep'}
    switch scatteringCase{1}
        case 'BI'
            delta_alpha = 0.1;
            alpha_f_arr = (0:delta_alpha:360)*pi/180;
%             f = [1e3, 3e3, 10e3, 30e3]; % frequencies in Hertz
            f = 1e3; % frequencies in Hertz
        case 'Sweep'
            alpha_f_arr = alpha_s;
            f = linspace(1e3,10e3,3000);
    end
    for ii = 1:length(models)
        model = models{ii};
        close all
        for BC = BCarr
            omega = 2*pi*f; % Angular frequency
            switch model
                case 'S1'
                    layer = setS1Parameters();
                case 'S3'
                    layer = setS3Parameters();
                case 'S5'
                    layer = setS5Parameters();
                case 'S13'
                    layer = setS13Parameters();
                case 'S15'
                    layer = setS15Parameters();
                case 'S35'
                    layer = setS35Parameters();
                case 'S135'
                    layer = setS135Parameters();
                case 'tripleShell'
                    layer = setTripleShellParameters();
            end
            layer{1}.calc_p_0 = true; % Calculate the far field pattern
            layer = defineBCstring(layer,BC{1});
            options = struct('BC', BC{1},...
                             'd_vec', d_vec, ...
                             'omega', omega, ...
                             'P_inc', 1);
            if strcmp(scatteringCase{1},'BI')
                if strcmp(model,'S5') && exist('miscellaneous/S5_SHBC_F30_specialValues.mat', 'file')
                    if 0
                        objf = @(alpha)-objFunc(alpha,beta_f,layer,options);
                        specialValues = findExtremas(objf, 0, 2*pi, 100000)';
                        delta = 1e-5*2*pi;
                        specialValues = sort([specialValues; (specialValues-delta); (specialValues+delta)]);
                        save('miscellaneous/S5_SHBC_F30_extremas', 'specialValues')
                    else
                        load('miscellaneous/S5_SHBC_F30_extremas.mat')
                    end
                    alpha_f_arr = unique(sort([specialValues', alpha_f_arr]));
                else
                    alpha_f_arr = (0:delta_alpha:360)*pi/180;
                end
            end
            X = [cos(beta_f)*cos(alpha_f_arr); cos(beta_f)*sin(alpha_f_arr); sin(beta_f)*ones(size(alpha_f_arr))]';

            layer{1}.X = X;

            if length(alpha_s) > 1
                aspect = 'S';
            else
                aspect = num2str(round(alpha_s*180/pi, 15, 'significant'));
            end
            if length(beta_s) > 1
                elevation = 'S';
            else
                elevation = num2str(round(beta_s*180/pi, 15, 'significant'));
            end
            if strcmp(scatteringCase{1}, 'Sweep')
                frequency = 'S';
            else
                frequency = num2str(f/1000);
            end
            [layer,~,flag] = e3Dss(layer, options);
            legendArr = cell(0,1);

            task.alpha_s = alpha_s;
            task.beta_s = beta_s;
            task.scatteringCase = scatteringCase{1};
            if strcmp(BC,'SHBC')
                colorArr = [0,70,147]/255;
            elseif strcmp(BC,'SSBC')
                colorArr = [178,0,0]/255;
            elseif strcmp(BC,'ESBC')
                colorArr = [59,124,37]/255;
            else
                colorArr = [149,49,157]/255;
            end
            task.comments = 'Exact solution';
            if strcmp(scatteringCase{1}, 'Sweep')
                saveName = [model '_' BC{1} '_' scatteringCase{1} '_A'  aspect '_E' elevation '_F' frequency];
                resultsFolderName = [resultsFolder '/' saveName];
                legendArr = {saveName};
                TS = 20*log10(abs(layer{1}.p_0)).'; 
                filename = [resultsFolder '/' saveName];
                task.saveName = saveName;
                task.f_arr = f;
                if any(~flag)
                    figure(30+ii)
                    TS_plot = TS(logical(~flag));
                    kR_01 = 2*pi*f(logical(~flag))/layer{1}.c*layer{1}.R;
                    printResultsToFile(filename, {'x',kR_01.', 'y', TS_plot, 'task',task})
                    plot(kR_01, TS_plot,'color',colorArr,'DisplayName', [model ' with ' BC{1}]);
                    hold on
                    xlabel('$kR_{0,1}$','interpreter','latex')
                    xlim([kR_01(1), kR_01(end)])
                    legend('off');
                    legend('show','Location','northeast');
                end
            elseif strcmp(scatteringCase{1}, 'BI')
                for i = 1:length(f)
                    figure(i)
                    saveName = [model '_' BC{1} '_' scatteringCase{1} '_A'  aspect '_E' elevation '_F' num2str(f(i)/1000)];
                    filename = [resultsFolder '/' saveName];
                    task.f = f(i);
                    task.saveName = saveName;
                    TS = 20*log10(abs(layer{1}.p_0(:,i))); 
                    if ~flag(i)
                        printResultsToFile(filename, {'x',180/pi*alpha_f_arr.', 'y', TS, 'task',task})
                        plot(alpha_f_arr*180/pi, TS,'color',colorArr,'DisplayName', [model ' with ' BC{1}]);
                        ylabel('TS','interpreter','latex')
                        xlabel('$\alpha_f$','interpreter','latex')
                        xlim([0, 360])
                        legend('off');
                        legend('show','Location','northeast');
                        hold on
                    end
                end
            end
        end
        if strcmp(scatteringCase{1}, 'Sweep')
            saveName = [model '_' scatteringCase{1} '_A'  aspect '_E' elevation '_F' frequency];
            filename = [resultsFolder '/' saveName];
            savefig([filename '.fig'])
        elseif strcmp(scatteringCase{1}, 'BI')
            for i = 1:length(f)
                figure(i)
                saveName = [model '_' scatteringCase{1} '_A'  aspect '_E' elevation '_F' num2str(f(i)/1000)];
                title([saveName(1:end-2) num2str(f(i)/1000)], 'interpreter', 'none')
                filename = [resultsFolder '/' saveName];
                savefig([filename(1:end-2) num2str(f(i)/1000) '.fig'])
            end
        end
    end
end
toc

function TS = objFunc(alpha,beta,layer,options)

options.Display = 'none';
layer{1}.X = [cos(beta)*cos(alpha); cos(beta)*sin(alpha); sin(beta)*ones(size(alpha))]';
layer = e3Dss(layer, options);
TS = 20*log10(abs(layer{1}.p_0)).';

end