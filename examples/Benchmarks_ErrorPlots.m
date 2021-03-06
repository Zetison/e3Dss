close all
clear all %#ok

startup
resultsFolder = [folderName '/Benchmarks_NearFieldPlots'];
if ~exist(resultsFolder, 'dir')
    mkdir(resultsFolder);
end

mpstartup
startMatlabPool
% mp = NaN;
%% Calculate errors
for useSymbolicPrecision = 1 %[0,1]
    if useSymbolicPrecision
        prec = 'mp';
%         prec = 'sym';
    else
        prec = 'double';
    end
    models = {'S1','S3','S5','S13','S15','S35','S135'};
%     models = {'S13','S15','S35','S135'};
%     models = {'S35'};
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
%     for i = 1:length(tasks)
    parfor i = 1:length(tasks)
        switch prec
            case 'single'
                Eps = 1e-7;
                PI = pi;b
                O = 0;
            case 'double'
                Eps = eps;
                PI = pi;
                O = 0;
            case 'sym'
                Eps = 1e-40;
%                 digits(2000);
                digits(40);
                PI = vpa('pi');
                O = vpa('0');
            case 'mp'
                Eps = 1e-40;
%                 Eps = 1e-30;
                mp.Digits(40);
%                 mp.Digits(2000);
                PI = mp('pi');
                O = mp('0');
        end
        model = tasks(i).model;
        ESBC = tasks(i).ESBC;
        SHBC = tasks(i).SHBC;
        SSBC = tasks(i).SSBC;
        tic
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
        end

        alpha_s = 240*PI/180;
        beta_s = 30*PI/180;
        d_vec = zeros(3,1,class(PI));
        d_vec(1) = -cos(beta_s)*cos(alpha_s);
        d_vec(2) = -cos(beta_s)*sin(alpha_s);
        d_vec(3) = -sin(beta_s);


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
        M = length(layer);
        
        E = [];
        nu = [];
        rho_s = [];
        c_f = [];
        R_i = [];
        R_o = [];
        for m = 1:M
            switch layer{m}.media
                case 'fluid'
                    c_f = [c_f, layer{m}.c_f];
                case 'solid'
                    E = [E, layer{m}.E];
                    nu = [nu, layer{m}.nu];
                    rho_s = [rho_s, layer{m}.rho];
                    if ~(strcmp(BC,'ESBC') && m == M)
                        R_i = [R_i,layer{m}.R_i];
                    end
                    R_o = [R_o,layer{m-1}.R_i];
            end
        end
        if strcmp(BC,'SHBC')
            R_o = [R_o, layer{end}.R_i];
        end
        
        npts_r = 4;
        npts_theta = 4;
        npts_phi = 4;
        for m = 1:M
            if m == 1
                r = linspaceHP(layer{m}.R_i, 2*layer{m}.R_i, npts_r);
            else
                r = linspaceHP(layer{m}.R_i, layer{m-1}.R_i, npts_r);
            end
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
        K = E./(3*(1-2*nu));
        G = E./(2*(1+nu));
        c_s_1 = sqrt((3*K+4*G)./(3*rho_s));
        c_s_2 = sqrt(G./rho_s);

        switch BC
            case {'ESBC', 'SHBC'}
                if isempty(R_i)
                    Upsilon = R_o./c_f;
                else
                    Upsilon = min([R_i./c_s_1(1:end-1), R_i./c_s_2(1:end-1), R_o./c_f]);
                end
            case 'SSBC'
                Upsilon = min([R_i./c_s_1, R_i./c_s_2, R_o./c_f]);
            case 'NNBC'
                Upsilon = min([R_i./c_s_1, R_i./c_s_2, R_o./c_f(1:end-1)]);
        end

        C = (layer{1}.R_i./layer{1}.c_f)^(3/2)/Upsilon^(1/2);
        if 0
            nFreqs = 2; %
%             f = 10.^linspaceHP(-log10(1e3*C),-log10(5e2*C),nFreqs);
            f = 10.^linspaceHP(-log10(1e3*C),log10(4e2/C),nFreqs);
        else
            nFreqs = 100; %
            f = 10.^linspaceHP(-log10(1e3*C),log10(4e2/C),nFreqs);
        end

        omega = 2*PI*f;
        options = struct('BC', BC, ...
                         'd_vec', d_vec, ...
                         'omega', omega, ...
                         'P_inc', ones(1,class(O)), ...
                         'prec', prec, ...
                         'Display','none',...
                         'Eps', Eps);
%         options.N_max = 2;
        for m = 1:M
            layer{m}.calc_errPresCond = true; 
            layer{m}.calc_errDispCond = true;  
            switch layer{m}.media
                case 'fluid'
                    layer{m}.calc_errHelm = true;  
                case 'solid'
                    layer{m}.calc_errNav = true;
            end
        end
        layer = e3Dss(layer, options);
        
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
                    err_navier1 = max([err_navier1; layer{m}.err_navier1],[],1);
                    err_navier2 = max([err_navier2; layer{m}.err_navier2],[],1);
            end
        end

        sc = f*C;
        figure
        if M == 1 && SHBC
            loglog(sc, err_helmholtz,'color',[0,70,147]/255,'DisplayName','Helmholtz')
            hold on
            loglog(sc, err_dc,'color',[149,49,157]/255,'DisplayName','Displacenemnt condition')
            legendArr = {'Helmholtz', 'DisplacementCond'};
            err = [err_helmholtz; err_dc];
        else
            loglog(sc, err_helmholtz,'color',[0,70,147]/255,'DisplayName','Helmholtz')
            hold on
            loglog(sc, err_navier1,'color',[178,0,0]/255,'DisplayName', 'Navier - $1^{\mathrm{st}}$ component')
            loglog(sc, err_navier2,'color',[59,124,37]/255,'DisplayName','Navier - $2^{\mathrm{nd}}$ component')
            loglog(sc, err_dc,'color',[149,49,157]/255,'DisplayName','Displacenemnt condition')
            loglog(sc, err_pc,'color',[247, 158,30]/255,'DisplayName','Pressure condition')
            legendArr = {'Helmholtz', 'Navier1', 'Navier2', 'DisplacementCond', 'PressureCond'};
            err = [err_helmholtz; err_navier1; err_navier2; err_dc; err_pc];
        end
        leg1 = legend('show','Location','northwest');
        set(leg1,'Interpreter','latex');
        filename = [resultsFolder '/errors_' model '_' BC '_Symbolic' num2str(useSymbolicPrecision)];

%         printResultsToFile(filename, {'x',double(sc.'), 'y', double(err.')})
        xlabel('$C f$','interpreter','latex')
        ylabel('Relative residual error')
        title(['Errors for model ' model '_' BC], 'interpreter', 'none')
        hold off
        if ~useSymbolicPrecision
            ylim([0.1*eps 1e2])
        end
        xlim([double(sc(1)), double(sc(end))])
        drawnow
        savefig([filename '.fig'])
        fprintf('Finished a case in %f seconds!\n\n', toc)
    end
end

