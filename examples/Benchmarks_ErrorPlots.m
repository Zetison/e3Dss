close all
clear all %#ok

startup
folderName = [homeDir '/Dropbox/Apps/Overleaf/createFigures/data/e3Dss_article2'];
resultsFolder = [folderName '/Benchmarks_ErrorPlots'];
if ~exist(resultsFolder, 'dir')
    mkdir(resultsFolder);
end

% startMatlabPool
% mp = NaN;
%% Calculate errors
for nu_a = [-1, 100]
    for useSymbolicPrecision = 0 %[0,1]
        if useSymbolicPrecision
            prec = 'mp';
    %         prec = 'sym';
            mpstartup
        else
            prec = 'double';
        end
        models = {'S1','S3','S5','S13','S15','S35','S135'};
    %     models = {'S13','S15','S35','S135'};
%         models = {'Skelton1997tao'};
%         models = {'Hetmaniuk2012raa'};
%         models = {'Sage1979mri'};
%         models = {'S135'};
%         models = {'S35'};
%         models = {'S1'};
        counter = 1;
        for i_model = 1:length(models)
            for ESBC = [0, 1]
                for SHBC = [0, 1]
                    for SSBC = [0, 1]
                        if ~(ESBC + SHBC + SSBC > 1)
                            tasks(counter).model = models{i_model};
                            tasks(counter).ESBC = ESBC;
                            tasks(counter).SHBC = SHBC;
                            tasks(counter).SSBC = SSBC;
                            counter = counter + 1;
                        end
                    end
                end
            end
        end
        for i = 1:length(tasks)
    %     parfor i = 1:length(tasks)
            nosymdigits = 40;
            switch prec
                case 'single'
                    Eps = 1e-7;
                    O = 0;
                case 'double'
                    Eps = eps;
                    O = 0;
                case 'sym'
                    Eps = vpa('1e-40');
    %                 digits(2000);
                    digits(nosymdigits);
                    O = vpa('0');
                case 'mp'
                    Eps = 10^(-mp(nosymdigits));
    %                 Eps = 1e-30;
                    mp.Digits(nosymdigits);
    %                 mp.Digits(2000);
                    O = mp('0');
            end
            PI = getPI(prec);
            model = tasks(i).model;
            ESBC = tasks(i).ESBC;
            SHBC = tasks(i).SHBC;
            SSBC = tasks(i).SSBC;
            tic
            r_s = NaN;
            switch model
                case 'S1'
                    layer = setS1Parameters(prec);
                case 'S3'
                    layer = setS3Parameters(prec);
                case 'S5'
                    layer = setS5Parameters(prec);
                case 'S13'
                    layer = setS13Parameters(prec);
                case 'S15'
                    layer = setS15Parameters(prec);
                case 'S35'
                    layer = setS35Parameters(prec);
                case 'S135'
                    layer = setS135Parameters(prec);
                case 'Skelton1997tao'
                    layer = setSkelton1997taoParameters();
                    if ~SSBC
                        continue
                    end
                case 'Hetmaniuk2012raa'
                    layer = setHetmaniukParameters();
                    if ~SSBC
                        continue
                    end
                case 'Sage1979mri'
                    layer = setSage1979mriParameters();
                    if any([SHBC,SSBC,ESBC])
                        continue
                    end
            end

            alpha_s = 240*PI/180;
            beta_s = 30*PI/180;
            d_vec = zeros(3,1,class(PI));
            d_vec(1) = -cos(beta_s)*cos(alpha_s);
            d_vec(2) = -cos(beta_s)*sin(alpha_s);
            d_vec(3) = -sin(beta_s);

            switch model
                case {'S1','S3','S5','S13','S15','S35','S135'}
                    if SHBC
                        layer = layer(1:end-2);
                        BC = 'SHBC';
                    elseif ESBC
                        layer = layer(1:end-1);
                        layer{end}.R_i = 0;
                        BC = 'ESBC';
                    elseif SSBC
                        layer = layer(1:end-1);
                        BC = 'SSBC';
                    else
                        BC = 'NNBC';
                    end
                case 'Hetmaniuk2012raa'
                    BC = 'SSBC';
                case 'Skelton1997tao'
                    BC = 'SSBC';
                case 'Sage1979mri'
                    BC = 'NNBC';
            end
            M = length(layer);

            npts_r = 4;
            npts_theta = 4;
            npts_phi = 4;
%             npts_r = 2;
%             npts_theta = 2;
%             npts_phi = 1;
            for m = 1:M
                if m == 1
                    r = linspaceHP(layer{m}.R_i, 2*layer{m}.R_i, npts_r);
                else
                    r = linspaceHP(layer{m}.R_i, layer{m-1}.R_i, npts_r);
                end
%                 if m == M
%                     r(1) = 0.001;
%                 end
                theta = linspaceHP(O,PI,npts_theta);
                phi = linspaceHP(O,2*PI,npts_phi);
                pts = zeros(length(r)*length(theta)*length(phi),3,class(PI));

                counter = 1;
                for ii = 1:length(r)
                    for jj = 1:length(theta)
                        for ll = 1:length(phi)
                            pts(counter,:) = r(ii)*[sin(theta(jj))*cos(phi(ll)), sin(theta(jj))*sin(phi(ll)), cos(theta(jj))];
                            counter = counter + 1;
                        end
                    end
                end
                [~, I, ~] = uniquetol(double(pts),10*eps,'ByRows',true, 'DataScale',max(max(abs(double(pts)))));
                layer{m}.X = pts(I,:);
            end
            R_1 = layer{1}.R_i;
            if 1
                nFreqs = 100;
%                 nFreqs = 4;
                kR_start = 1e-1;
                kR_end = 1e4;
%                 kR_end = 1e0;
                kR = 10.^linspaceHP(log10(kR_start),log10(kR_end),nFreqs);
                k = kR/R_1;
                omega = k*layer{1}.c_f;
                f = omega/(2*PI);
            else
                nFreqs = 100;
                f_max = 25e4;
                f_max = 10000;
                f = linspace(f_max/nFreqs,f_max,nFreqs);
                omega = 2*PI*f;
                k = omega/layer{1}.c_f;
                kR = k*R_1;
            end
%             N_max = 1122;
            N_max = Inf;
            options = struct('BC', BC, ...
                             'd_vec', d_vec, ...
                             'omega', omega, ...
                             'P_inc', ones(1,class(O)), ...
                             'prec', prec, ...
                             'Display','final',...
                             'debug',0, ...
                             'nu_a',nu_a, ...
                             'N_max', N_max, ...
                             'Eps', Eps);
    %         options.N_max = 2;
            for m = 1:M
                layer{m}.calc_err_pc = true; 
                layer{m}.calc_err_dc = true;  
                switch layer{m}.media
                    case 'fluid'
                        layer{m}.calc_err_helmholtz = true;  
                    case 'solid'
                        layer{m}.calc_err_navier = true(1,2);
                end
            end
            [layer,N_eps,flag,relTermMaxArr] = e3Dss(layer, options);

            err_navier1 = zeros(1,nFreqs,class(PI));
            err_navier2 = zeros(1,nFreqs,class(PI));
            err_helmholtz = zeros(1,nFreqs,class(PI));
            err_pc = zeros(1,nFreqs,class(PI));
            err_dc = zeros(1,nFreqs,class(PI));
            for m = 1:M
                isSphere = layer{m}.R_i == 0;
                if ~isSphere && ~(m == M && SHBC)
                    err_pc = max([err_pc; layer{m}.err_pc],[],1);
                end
                if ~isSphere && ~(m == M && SSBC)
                    err_dc = max([err_dc; layer{m}.err_dc],[],1);
                end
                switch layer{m}.media
                    case 'fluid'
                        err_helmholtz = max([err_helmholtz; layer{m}.err_helmholtz],[],1);
                    case 'solid'
                        err_navier1 = max([err_navier1; layer{m}.err_navier{1}],[],1);
                        err_navier2 = max([err_navier2; layer{m}.err_navier{2}],[],1);
                end
            end

            figure
            if M == 1 && SHBC
                loglog(kR, err_helmholtz,'color',[0,70,147]/255,'DisplayName','Helmholtz')
                hold on
                loglog(kR, err_dc,'color',[149,49,157]/255,'DisplayName','Displacenemnt condition')
                legendArr = {'Helmholtz', 'DisplacementCond'};
                err = [err_helmholtz; err_dc];
            else
                loglog(kR, err_helmholtz,'color',[0,70,147]/255,'DisplayName','Helmholtz')
                hold on
                loglog(kR, err_navier1,'color',[178,0,0]/255,'DisplayName', 'Navier - $1^{\mathrm{st}}$ component')
                loglog(kR, err_navier2,'color',[59,124,37]/255,'DisplayName','Navier - $2^{\mathrm{nd}}$ component')
                loglog(kR, err_dc,'color',[149,49,157]/255,'DisplayName','Displacenemnt condition')
                loglog(kR, err_pc,'color',[247, 158,30]/255,'DisplayName','Pressure condition')
                legendArr = {'Helmholtz', 'Navier1', 'Navier2', 'DisplacementCond', 'PressureCond'};
                err = [err_helmholtz; err_navier1; err_navier2; err_dc; err_pc];
            end
            leg1 = legend('show','Location','northwest');
            set(leg1,'Interpreter','latex');
            useScaling = nu_a ~= -1;
            
            filename = [resultsFolder '/errors_' model '_' BC '_Symbolic' num2str(useSymbolicPrecision) '_Scaling' num2str(useScaling)];

            printResultsToFile(filename, {'x',double(kR.'), 'y', double(err.'), 'xlabel','kR', 'ylabel',legendArr})
            xlabel('$k_1 R_1$','interpreter','latex')
            ylabel('Relative residual error')
            title(['Errors for model ' model '_' BC '_Symbolic' num2str(useSymbolicPrecision) '_Scaling' num2str(useScaling)], 'interpreter', 'none')
            hold off
            if ~useSymbolicPrecision
                ylim([0.1*eps 1e2])
            end
            xlim([double(kR(1)), double(kR(end))])
            drawnow
            savefig([filename '.fig'])
            fprintf('Finished a case in %f seconds!\n\n', toc)
            flags = find(flag);
            if any(flags)
                fprintf('Flags at f > %f kHz (kR > %f)\n\n', f(flags(1))/1000, kR(flags(1)))
            end
        end
    end
end
