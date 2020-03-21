close all
clear all %#ok

addpath ..
addpath ../utils
addpath ../models

% pathToResults = '../../../results/e3Dss/';
pathToResults = '../results/';
% mpstartup
startMatlabPool
mp = NaN;
%% Calculate errors
for useSymbolicPrecision = 0 %[0,1]
    if useSymbolicPrecision
%         prec = 'mp';
        prec = 'sym';
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
    for i = 1:length(tasks)
%     parfor i = 1:length(tasks)
        switch prec
            case 'single'
                Eps = 1e-7;
                PI = pi;
                O = 0;
            case 'double'
                Eps = eps;
                PI = pi;
                O = 0;
            case 'sym'
                Eps = 1e-50;
                digits(2000);
                PI = vpa('pi');
                O = vpa('0');
            case 'mp'
                Eps = 1e-50;
                mp.Digits(2000);
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
                Upsilon = min([R_i./c_s_1(1:end-1), R_i./c_s_2(1:end-1), R_o./c_f]);
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
                         'Eps', Eps);
        for m = 1:M
            layer{m}.calc_errPresCond = true; 
            layer{m}.calc_errDispCond = true;  
            switch layer{m}.media
                case 'fluid'
                    layer{m}.calc_errHelm = true;  
                    layer{m}.calc_p_0 = false;
                case 'solid'
                    layer{m}.calc_du = true(3,3);
                    layer{m}.calc_sigma = true(1,6);
            end
        end
        layer = e3Dss(layer, options);
        
        err_stress_xx = zeros(M,nFreqs);
        err_stress_yy = zeros(M,nFreqs);
        err_stress_zz = zeros(M,nFreqs);
        err_stress_yz = zeros(M,nFreqs);
        err_stress_xz = zeros(M,nFreqs);
        err_stress_xy = zeros(M,nFreqs);
        for m = 1:M
            n_X = size(layer{m}.X,1);
            if strcmp(layer{m}.media, 'solid')
                strain = cell(6,1);
                for ii = 1:6
                    strain{ii} = zeros(n_X,nFreqs,class(PI));
                end                         
                stress = strain;             
                du_X = cell(3,3);
                du_X{1,1} = layer{m}.du_xdx;
                du_X{1,2} = layer{m}.du_xdy;
                du_X{1,3} = layer{m}.du_xdz;
                du_X{2,1} = layer{m}.du_ydx;
                du_X{2,2} = layer{m}.du_ydy;
                du_X{2,3} = layer{m}.du_ydz;
                du_X{3,1} = layer{m}.du_zdx;
                du_X{3,2} = layer{m}.du_zdy;
                du_X{3,3} = layer{m}.du_zdz;
                vgtinv = [1 6 5;
                          6 2 4;
                          5 4 3];
                for ii = 1:3
                    for jj = ii:3
                        strain{vgtinv(ii,jj)} = 0.5*(du_X{ii,jj}+du_X{jj,ii});
                    end
                end
                for ii = 1:3
                    for jj = 1:3
                        if ii == jj
                            stress{ii} = stress{ii} + (layer{m}.K+4*layer{m}.G/3)*strain{jj};
                        else
                            stress{ii} = stress{ii} + (layer{m}.K-2*layer{m}.G/3)*strain{jj};
                        end
                    end
                end
                for ii = 4:6
                    stress{ii} = 2*layer{m}.G*strain{ii};
                end
                err_stress_xx(m,:) = max(abs(stress{1} - layer{m}.sigma_xx),[],1)./max(abs(layer{m}.sigma_xx),[],1);
                err_stress_yy(m,:) = max(abs(stress{2} - layer{m}.sigma_yy),[],1)./max(abs(layer{m}.sigma_yy),[],1);
                err_stress_zz(m,:) = max(abs(stress{3} - layer{m}.sigma_zz),[],1)./max(abs(layer{m}.sigma_zz),[],1);
                err_stress_yz(m,:) = max(abs(stress{4} - layer{m}.sigma_yz),[],1)./max(abs(layer{m}.sigma_yz),[],1);
                err_stress_xz(m,:) = max(abs(stress{5} - layer{m}.sigma_xz),[],1)./max(abs(layer{m}.sigma_xz),[],1);
                err_stress_xy(m,:) = max(abs(stress{6} - layer{m}.sigma_xy),[],1)./max(abs(layer{m}.sigma_xy),[],1);
            end
        end
        err_stress_xx = max(err_stress_xx,[],1);
        err_stress_yy = max(err_stress_yy,[],1);
        err_stress_zz = max(err_stress_zz,[],1);
        err_stress_yz = max(err_stress_yz,[],1);
        err_stress_xz = max(err_stress_xz,[],1);
        err_stress_xy = max(err_stress_xy,[],1);

        sc = f*C;
        figure(i)
        if ~(M == 1 && SHBC)
            loglog(sc, err_stress_xx,'color',[0,70,147]/255,'DisplayName', '$\sigma_{xx}$')
            hold on
            loglog(sc, err_stress_yy,'color',[178,0,0]/255,'DisplayName', 'Error in $\sigma_{yy}$')
            loglog(sc, err_stress_zz,'color',[59,124,37]/255,'DisplayName', 'Error in $\sigma_{zz}$')
            loglog(sc, err_stress_yz,'color',[149,49,157]/255,'DisplayName', 'Error in $\sigma_{yz}$')
            loglog(sc, err_stress_xz,'color',[247, 158,30]/255,'DisplayName', 'Error in $\sigma_{xz}$')
            loglog(sc, err_stress_xy,'color',[0,172,239]/255,'DisplayName', 'Error in $\sigma_{xy}$')
            err = [err_stress_xx; err_stress_yy; err_stress_zz; err_stress_yz; err_stress_xz; err_stress_xy];
        else
            continue
        end
        leg1 = legend('show','Location','northwest');
        set(leg1,'Interpreter','latex');
        filename = [pathToResults 'sigmaErrors_' model '_' BC '_Symbolic' num2str(useSymbolicPrecision)];
        
        xlabel('$C f$','interpreter','latex')
        ylabel('Relative error')
        title(['Errors for model ' model '_' BC], 'interpreter', 'none')
        hold off
        if ~useSymbolicPrecision
            ylim([0.1*eps 1e2])
        end
        xlim([double(sc(1)), double(sc(end))])
        drawnow
%                     savefig([filename '.fig'])
        fprintf('Finished a case in %f seconds!\n\n', toc)
    end
end




