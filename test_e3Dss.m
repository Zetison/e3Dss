close all
clear all %#ok

addpath utils
addpath models

example = 'testEnvironment'; % use this example to play around with the code
% example = 'create Far-Field Pattern Plots'; % for ALL benchmarks
% example = 'create error plots'; % for ALL benchmarks
% example = 'Chang (1994)';
% example = 'Ihlenburg (1998)';
% example = 'Fender (1972)';
% example = 'Bessel functions for large n';
% example = 'Bessel functions upper bound';
% example = 'Plot in Paraview';
% example = 'Solid sphere';
% example = 'Solid sphere inside spherical shell';
% example = 'create spy matrix';
% example = 'time domain solution - 1D Visualization';

switch example
    case 'testEnvironment'
        %% Test benchmark models
        alpha_s = 240*pi/180;
%         alpha_s = 90*pi/180;
        beta_s = 30*pi/180;
%         beta_s = 0;
        beta_f = beta_s;
%         beta_f_arr = beta_s;
%         beta_f = 0*pi/180;
%         beta_f = 0;
        beta_f_arr = beta_f;
        alpha_f = alpha_s;

        d_vec = -[cos(beta_s)*cos(alpha_s);
                  cos(beta_s)*sin(alpha_s);
                  sin(beta_s)];

%         scatteringCase = 'Ray';
        scatteringCase = 'BI';
%         scatteringCase = 'Sweep';
        switch scatteringCase
            case 'Ray'
                alpha_f_arr = alpha_f;
                beta_f_arr = beta_f;
                f = 10e3; % frequencies in Hertz
            case 'BI'
                delta_alpha = 1;
                alpha_f_arr = (0:delta_alpha:360)*pi/180;
                alpha_f_arr = 0;
                f = 30e3; % frequencies in Hertz
            case 'Sweep'
                alpha_f_arr = alpha_s;
                delta_f = 1;
                Nfreqs = 3;
                c_f = 1524;
                kstart = 0.01;
                kend = 2;
                f = linspace(kstart*c_f/(2*pi),kend*c_f/(2*pi),Nfreqs);
        end
        model = 'S5';
        % model = 'test';
        
        switch model
            case 'S1'
                setS1Parameters
            case 'S3'
                setS3Parameters
            case 'S5'
                setS5Parameters
            case 'S13'
                setS13Parameters
            case 'S15'
                setS15Parameters
            case 'S35'
                setS35Parameters
            case 'S135'
                setS135Parameters
            case 'test'
                P_inc = 1; % Amplitude of incident wave
                rho_s = [7850, 7850, 7850]; % Density of solid
                rho_f = [1000, 1000, 1.2, 1.2]; % Density of fluids
                c_f = [1500, 1500, 340, 340];  % Speed of sound in fluid domains
                t = [0.008, 0.02, 0.05]; % The thickness of the sphere
                E = [210e9, 210e9, 210e9]; % Youngs modulus of elastic material
                nu = [0.3, 0.3, 0.3]; % Poisson ratio of elastic material
                R_o = [5, 1, 0.5];
            case 'IL'
                setIhlenburgParameters

        end
        SHBC = true;
        ESBC = false;
        SSBC = false;
        R_i = R_o - t;
        omega = 2*pi*f; % Angular frequency
        
        defineBCstring
        
        options = struct('d_vec', d_vec,... 
                         'omega', omega, ...
                         'R_i', R_i, ...
                         'R_o', R_o, ...
                         'P_inc', P_inc, ...
                         'E', E, ...
                         'nu', nu, ...
                         'rho_s', rho_s, ...
                         'rho_f', rho_f, ...
                         'calc_farField', 1, ...
                         'c_f', c_f);
                 
%         r_vec = linspace(R_o(1),2*R_o(1),1000).';
%         v = r_vec*[cos(beta_f)*cos(alpha_f_arr); cos(beta_f)*sin(alpha_f_arr); sin(beta_f)*ones(size(alpha_f_arr))]'; % points valid for far-field pattern calculation
        load('results/S5_SHBC_F30_specialValues')
        alpha_f_arr = unique(sort([specialValues', linspace(0,2*pi,1000)]));
        v = R_o(1)*[cos(beta_f)*cos(alpha_f_arr); cos(beta_f)*sin(alpha_f_arr); sin(beta_f)*ones(size(alpha_f_arr))]';
        data = e3Dss(v, options);
        
        plot(180/pi*alpha_f_arr, 20*log10(abs(data(1).p)))
        xlim([0,360])
%         switch scatteringCase
%             case 'BI'
%                 TS = 20*log10(abs(data(1).p)); 
%                 plot(alpha_f_arr*180/pi, TS);
%                 ylabel('TS')
%                 xlabel('\alpha_f')
%                 xlim([0, 360])
%             case 'Sweep'
%                 indices = ~data(1).flag;
%                 TS = 20*log10(abs(data(1).p(indices))); 
%                 plot(f(indices), TS);
%                 ylabel('TS')
%                 xlabel('f')
%                 xlim([f(1), f(end)])
%         end
%         
%         k = omega./c_f;
%         p_inc = @(v) P_inc*exp(1i*k*dot3(v,d_vec));
%         gp_inc = @(v) 1i*p_inc(v)*k(1);
%         plot(f, 20*log10(abs(data(1).p)))
%         plot(f, abs(data(1).p+p_inc(v)))
%         plot(k, 20*log10(abs(data(1).p)))
%         plot(r_vec, real(data(1).p+p_inc(v)), r_vec, real(data(1).dpdz+gp_inc(v)))
%         legend('Pressure field', 'Velocity field in z-direction')
        legend('S5 with SHBC f = 30Hz')
        ylabel('TS')
        xlabel('k')
        
    case 'create Far-Field Pattern Plots'
        % Calculate Far-Field pattern (BI and Sweep)
        tic
        alpha_s = 240*pi/180;
        beta_s = 30*pi/180;
        beta_f = beta_s;
        beta_f_arr = beta_s;
        d_vec = -[cos(beta_s)*cos(alpha_s);
                  cos(beta_s)*sin(alpha_s);
                  sin(beta_s)];
        models = {'S1','S3','S5','S35','S15','S13','S135'};
        
        for scatteringCase = {'BI'}
            switch scatteringCase{1}
                case 'BI'
                    delta_alpha = 0.1;
                    alpha_f_arr = (0:delta_alpha:360)*pi/180;
                    f = [1e3, 3e3, 10e3, 30e3]; % frequencies in Hertz
%                     f = [1e3, 10e3]; % frequencies in Hertz
                case 'Sweep'
                    alpha_f_arr = alpha_s;
                    delta_f = 0.1;
                    f = linspace(1e3,10e3,3000);
%                     f = 1e3:delta_f:30e3;
%                     f = 1e3:delta_f:5e3;
            end
            for ii = 1:length(models)
                model = models{ii};
                close all
                for SHBC = [0, 1]
                    for ESBC = [0, 1]
                        for SSBC = [0, 1]
                            if ~(ESBC + SHBC + SSBC > 1)
                                if strcmp(model,'S5')
                                    load('results/S5_SHBC_F30_specialValues')
                                    alpha_f_arr = unique(sort([specialValues', alpha_f_arr]));
                                else
                                    alpha_f_arr = (0:delta_alpha:360)*pi/180;
                                end
                                v = [cos(beta_f)*cos(alpha_f_arr); cos(beta_f)*sin(alpha_f_arr); sin(beta_f)*ones(size(alpha_f_arr))]'; 

                                switch model
                                    case 'S1'
                                        setS1Parameters
                                    case 'S3'
                                        setS3Parameters
                                    case 'S5'
                                        setS5Parameters
                                    case 'S13'
                                        setS13Parameters
                                    case 'S15'
                                        setS15Parameters
                                    case 'S35'
                                        setS35Parameters
                                    case 'S135'
                                        setS135Parameters
                                    case 'tripleShell'
                                        setTripleShellParameters
                                end
                                R_i = R_o - t;
                                omega = 2*pi*f; % Angular frequency
                                
                                defineBCstring

                                options = struct('d_vec', d_vec, ...
                                                 'omega', omega, ...
                                                 'R_i', R_i, ...
                                                 'R_o', R_o, ...
                                                 'P_inc', P_inc, ...
                                                 'E', E, ...
                                                 'nu', nu, ...
                                                 'rho_s', rho_s, ...
                                                 'rho_f', rho_f, ...
                                                 'c_f', c_f, ...
                                                 'calc_farField', 1);

                                % Create folders
                                folderName = 'results';
                                if ~exist(folderName, 'dir')
                                    mkdir(folderName);
                                end

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
                                data = e3Dss(v, options);
                                legendArr = cell(0,1);

                                varCol.alpha_s = alpha_s;
                                varCol.beta_s = beta_s;
                                varCol.scatteringCase = scatteringCase{1};
                                if strcmp(scatteringCase{1}, 'Sweep')
                                    saveName = [model '_' BC '_' scatteringCase{1} '_A'  aspect '_E' elevation '_F' frequency];
                                    resultsFolderName = [folderName '/' saveName];
                                    legendArr = {saveName};
                                    TS = 20*log10(abs(data(1).p)).'; 
                                    filename = [folderName '/' saveName];
                                    varCol.saveName = saveName;
                                    varCol.f_arr = f;
                                    if any(~data(1).flag)
                                        TS_plot = TS(logical(~data(1).flag));
                                        kR_01 = 2*pi*f(logical(~data(1).flag))/c_f(1)*R_o(1);
                                        printResultsToFile(filename, kR_01.', TS_plot, varCol, 1, 0, 'NTNU_FFI', 'Exact solution')
                                        plot(kR_01, TS_plot);
                                        hold on
                                        xlabel('kR_{0,1}')
                                        xlim([kR_01(1), kR_01(end)])
                                    end
                                elseif strcmp(scatteringCase{1}, 'BI')
                                    for i = 1:length(f)
                                            
                                        figure(i)
                                        saveName = [model '_' BC '_' scatteringCase{1} '_A'  aspect '_E' elevation '_F' num2str(f(i)/1000)];
                                        filename = [folderName '/' saveName];
                                        varCol.f = f(i);
                                        varCol.saveName = saveName;
                                        TS = 20*log10(abs(data(1).p(:,i))); 
                                        if ~data(1).flag(i)
                                            printResultsToFile(filename, 180/pi*alpha_f_arr.', TS, varCol, 1, 0, 'NTNU_FFI', 'Exact solution')
                                            plot(alpha_f_arr*180/pi, TS);
                                            ylabel('TS')
                                            xlabel('\alpha_f')
                                            xlim([0, 360])
                                            hold on
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
                if strcmp(scatteringCase{1}, 'Sweep')
                    title(filename(9:end), 'interpreter', 'none')
                    legend({'NNBC','SSBC','ESBC','SHBC'}, 'interpreter', 'none')
                    savefig([filename '.fig'])
                elseif strcmp(scatteringCase{1}, 'BI')
                    for i = 1:length(f)
                        figure(i)
                        legend({'NNBC','SSBC','ESBC','SHBC'}, 'interpreter', 'none')
                        title([filename(9:end-2) num2str(f(i)/1000)], 'interpreter', 'none')
                        savefig([filename(1:end-2) num2str(f(i)/1000) '.fig'])
                    end
                end
            end
        end
        toc
    case 'create error plots'
        %% Calculate errors
        useSymbolicPrecision = 0;
        if useSymbolicPrecision
            Eps = 1e-40;
        else
            Eps = eps;
        end
        
        models = {'S1','S3','S5','S35','S15','S13','S135'};
%         models = {'S135'};
        
        poolobj = gcp('nocreate');
        noCoresToUse = feature('numCores');
        if noCoresToUse > 7
            noCoresToUse = 7;
        end
        if isempty(poolobj)
            parpool(noCoresToUse, 'IdleTimeout', Inf)
        end
        tic
        noTestCases = 1;
        for i = 1:length(models)
            model = models{i};
            for ESBC = [0, 1]
                for SHBC = [0, 1]
                    for SSBC = [0, 1]
                        if ~(ESBC + SHBC + SSBC > 1)
                            tic
                            switch model
                                case 'S1'
%                                     setS1Parameters
                                    P_inc = 1; % Amplitude of incident wave
                                    rho_s = 7850; % Density of solid
                                    rho_f = [1000, 1.2]; % Density of fluids
                                    c_f = [1500, 340];  % Speed of sound in fluid domains
                                    t = 0.05; % The thickness of the sphere
                                    E = 210e9; % Youngs modulus of elastic material
                                    nu = 0.3; % Poisson ratio of elastic material

                                    R_o = 1; % Distance to far field point

                                case 'S3'
%                                     setS3Parameters
                                    P_inc = 1; % Amplitude of incident wave
                                    rho_s = 7850; % Density of solid
                                    rho_f = [1000, 1.2]; % Density of fluids
                                    c_f = [1500, 340];  % Speed of sound in fluid domains
                                    t = 0.02; % The thickness of the sphere
                                    E = 210e9; % Youngs modulus of elastic material
                                    nu = 0.3; % Poisson ratio of elastic material

                                    R_o = 3; % Distance to far field point
                                case 'S5'
%                                     setS5Parameters
                                    P_inc = 1; % Amplitude of incident wave
                                    rho_s = 7850; % Density of solid
                                    rho_f = [1000, 1000]; % Density of fluids
                                    c_f = [1500, 1500];  % Speed of sound in fluid domains
                                    t = 0.008; % The thickness of the sphere
                                    E = 210e9; % Youngs modulus of elastic material
                                    nu = 0.3; % Poisson ratio of elastic material

                                    R_o = 5; % Distance to far field point
                                case 'S13'
%                                     setS13Parameters
                                    P_inc = 1; % Amplitude of incident wave
                                    rho_s = [7850, 7850]; % Density of solid
                                    rho_f = [1000, 1.2, 1.2]; % Density of fluids
                                    c_f = [1500, 340, 340];  % Speed of sound in fluid domains
                                    t = [0.02 0.05]; % The thickness of the sphere
                                    E = [210e9, 210e9]; % Youngs modulus of elastic material
                                    nu = [0.3, 0.3]; % Poisson ratio of elastic material

                                    R_o = [3, 1];
                                case 'S15'
%                                     setS15Parameters
                                    P_inc = 1; % Amplitude of incident wave
                                    rho_s = [7850, 7850]; % Density of solid
                                    rho_f = [1000, 1000, 1.2]; % Density of fluids
                                    c_f = [1500, 1500, 340];  % Speed of sound in fluid domains
                                    t = [0.008, 0.05]; % The thickness of the sphere
                                    E = [210e9, 210e9]; % Youngs modulus of elastic material
                                    nu = [0.3, 0.3]; % Poisson ratio of elastic material

                                    R_o = [5, 1];

                                case 'S35'
%                                     setS35Parameters
                                    P_inc = 1; % Amplitude of incident wave
                                    rho_s = [7850, 7850]; % Density of solid
                                    rho_f = [1000, 1000, 1.2]; % Density of fluids
                                    c_f = [1500, 1500, 340];  % Speed of sound in fluid domains
                                    t = [0.008, 0.02]; % The thickness of the sphere
                                    E = [210e9, 210e9]; % Youngs modulus of elastic material
                                    nu = [0.3, 0.3]; % Poisson ratio of elastic material

                                    R_o = [5, 3]; 

                                case 'S135'
%                                     setS135Parameters
                                    P_inc = 1; % Amplitude of incident wave
                                    rho_s = [7850, 7850, 7850]; % Density of solid
                                    rho_f = [1000, 1000, 1.2, 1.2]; % Density of fluids
                                    c_f = [1500, 1500, 340, 340];  % Speed of sound in fluid domains
                                    t = [0.008, 0.02, 0.05]; % The thickness of the sphere
                                    E = [210e9, 210e9, 210e9]; % Youngs modulus of elastic material
                                    nu = [0.3, 0.3, 0.3]; % Poisson ratio of elastic material

                                    R_o = [5, 1, 0.5];

                            end
                            R_i = R_o - t; % Inner radius of shell
                            alpha_s = 240*pi/180;
                            beta_s = 30*pi/180;
        
                            if useSymbolicPrecision
                                P_inc = vpa(P_inc);
                                rho_s = vpa(rho_s);
                                rho_f = vpa(rho_f);
                                c_f = vpa(c_f);
                                t = vpa(t);
                                E = vpa(E);
                                nu = vpa(nu);
                                R_o = vpa(R_o); % Density of solid
                                R_i = R_o - t;
                                alpha_s = 240*vpa(pi)/180;
                                beta_s = 30*vpa(pi)/180;
                            end
                            d_vec = -[cos(beta_s)*cos(alpha_s);
                                      cos(beta_s)*sin(alpha_s);
                                      sin(beta_s)];

%                             defineBCstring
                            M = length(R_o);
                            if SHBC
                                E = E(1:end-1);
                                rho_s = rho_s(1:end-1);
                                rho_f = rho_f(1:end-1);
                                c_f = c_f(1:end-1);
                                nu = nu(1:end-1);
                                R_i = R_i(1:end-1);
                                BC = 'SHBC';
                            elseif ESBC
                                rho_f = rho_f(1:end-1);
                                c_f = c_f(1:end-1);
                                R_i = R_i(1:end-1);
                                BC = 'ESBC';
                            elseif SSBC
                                rho_f = rho_f(1:end-1);
                                c_f = c_f(1:end-1);
                                BC = 'SSBC';
                            else
                                BC = 'NNBC';
                            end

                            M = length(R_o);
                            noDomains = 2*M+1-ESBC-2*SHBC-SSBC;
                            vv = cell(noDomains,1);
                            m = 1;
                            npts_r = 4;
                            npts_theta = 4;
                            npts_phi = 4;
                            for j = 1:noDomains
                                if j == 1
                                    r = linspace(R_o(1), 2*R_o(1), npts_r);
                                elseif mod(j,2)
                                    if m == M+1
                                        r = linspace(0, R_i(m-1), npts_r-1);
                                    else
                                        r = linspace(R_o(m), R_i(m-1), npts_r);
                                    end
                                else
                                    if ESBC && m == M
                                        r = linspace(0, R_o(m), npts_r-1);    
                                    else
                                        r = linspace(R_i(m), R_o(m), npts_r);     
                                    end
                                    m = m + 1;
                                end
                                r = r.';
                                theta = linspace(0,pi,npts_theta);
                                phi = zeros(1,1,npts_phi);
                                phi(1,1,1:npts_phi) = linspace(0,2*pi,npts_phi);
                                if useSymbolicPrecision
                                    r = vpa(r);
                                    theta = vpa(theta);
                                    phi = vpa(phi);
                                end
                                R = repmat(r,1,length(theta),length(phi));
                                Theta = repmat(theta,length(r),1,length(phi));
                                Phi = repmat(phi,length(r),length(theta),1);
                                v_x = R.*sin(Theta).*cos(Phi);
                                v_y = R.*sin(Theta).*sin(Phi);
                                v_z = R.*cos(Theta);
                                v_x = v_x(:);
                                v_y = v_y(:);
                                v_z = v_z(:);
                                pts = double(abs([v_x, v_y, v_z]));
                                [~, I, ~] = uniquetol(pts,10*eps,'ByRows',true, 'DataScale',max(max(abs(pts))));
                                vv{j} = [v_x(I), v_y(I), v_z(I)];
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
                            
                            
                            
                            
                            C_eps = (R_o(1)./c_f(1))^(3/2) ...
                                    /Upsilon^(1/2);
                            nFreqs = 100;
                            f = 10.^linspace(log10(1e-2/C_eps),log10(2e3/C_eps),nFreqs);
                            if useSymbolicPrecision
                                f = vpa(f);
                            end
                            omega = 2*pi*f;
                            options = struct('d_vec', d_vec, ...
                                             'omega', omega, ...
                                             'R_i', R_i, ...
                                             'R_o', R_o, ...
                                             'P_inc', P_inc, ...
                                             'E', E, ...
                                             'nu', nu, ...
                                             'rho_s', rho_s, ...
                                             'rho_f', rho_f, ...
                                             'c_f', c_f, ...
                                             'Eps', Eps);
        
                            options.calc_u_x = 1;
                            options.calc_u_y = 1;
                            options.calc_u_z = 1;
                            options.calc_dpdx = 1;
                            options.calc_dpdy = 1;
                            options.calc_dpdz = 1;
                            options.calc_p_laplace = 1;
                            options.calc_sigma_xx = 1;
                            options.calc_sigma_yy = 1;
                            options.calc_sigma_zz = 1;
                            options.calc_sigma_yz = 1;
                            options.calc_sigma_xz = 1;
                            options.calc_sigma_xy = 1;
                            options.calc_sigma_rr = 1;
                            options.calc_navier = 1;
                            options.calc_errorsNavier = 1;
                            options.calc_errorsDisplacementCondition = 1;
                            options.calc_errorsPressureCondition = 1;
                            options.calc_errorsHelmholtz = 1;
                            
                            options.useSymbolicPrecision = useSymbolicPrecision;
                            data = e3Dss(vv, options);
                            err_navier1 = zeros(M,nFreqs);
                            err_navier2 = zeros(M,nFreqs);
                            err_helmholtz = zeros(M+1,nFreqs);
                            err_pc = zeros(M,nFreqs);
                            err_dc = zeros(M,nFreqs);
                            for m = 1:M+1
                                if m ~= M+1
                                    err_dc(m,:) = data(m).err_dc;
                                    if ~(SHBC && m == M)
                                        err_pc(m,:) = data(m).err_pc;
                                        err_navier1(m,:) = data(m).err_navier1;
                                        err_navier2(m,:) = data(m).err_navier2;
                                    end
                                end
                                if ~(m == M+1 && (SHBC || SSBC || ESBC))
                                    err_helmholtz(m,:) = data(m).err_helmholtz;
                                end
                            end
                            err_navier1 = max(err_navier1,[],1);
                            err_navier2 = max(err_navier2,[],1);
                            err_helmholtz = max(err_helmholtz,[],1);
                            err_pc = max(err_pc,[],1);
                            err_dc = max(err_dc,[],1);
                            
                            sc = f*C_eps;
                            if M == 1 && SHBC
                                loglog(sc, err_helmholtz, sc, err_dc)
                                legendArr = {'Helmholtz', 'DispCondition'};
                                legend(legendArr, 'interpreter', 'none')
                                err = [err_helmholtz; err_dc];
                            else
                                loglog(sc, err_pc, sc, err_navier1, sc, err_navier2, sc, err_helmholtz, sc, err_dc)
                                legendArr = {'PressureCondition', 'Navier1', 'Navier2', 'Helmholtz', 'DispCondition'};
                                legend(legendArr, 'interpreter', 'none','Location','north')
                                err = [err_pc; err_navier1; err_navier2; err_helmholtz; err_dc];
                            end
                            filename = ['results/errors_' model '_' BC '_Symbolic' num2str(useSymbolicPrecision)];
                            
                            printResultsToFile(filename, double(sc.'), double(err.'), [], 0, 1, [], [], {'C_ef'},legendArr)
                            xlabel('C_\epsilon f')
                            ylabel('Relative residual error')
                            title(['Errors for model ' model '_' BC], 'interpreter', 'none')
                            hold off
                            if ~useSymbolicPrecision
                                ylim([0.1*eps 1])
                            end
                            xlim([double(sc(1)), double(sc(end))])
                            drawnow
                            savefig([filename '.fig'])
                            noTestCases = noTestCases + 1;
                            fprintf('Finished a case in %f seconds!\n\n', toc)
                        end
                    end
                end
            end
        end
        toc
    case 'Chang (1994)'
        %% Chang and Demkowiz (1994) example
        P_inc = 1; % Amplitude of incident wave
        rho_f = 1000; % Density of fluids
        c_f = 1460; % Speed of sound in fluid domains
        rho_s = 7800; % Density of solid
        E = 2.0e11;
        nu = 0.3;
        a = 1; % mid surface radius
        h = 0.01*a;
        R_i = a - h/2; % Inner radius of shell
        R_o = a + h/2;  % Outer radius of shell
        
        %%%%%%%%%
        k = [15, 20];
        omega = k*c_f(1);
        
        d_vec = [0,0,1].';
        p_inc = @(v) P_inc*exp(1i*dot3(v,d_vec)*k);
        p_tot_Chang15 = importdata('models/Chang1994voa2/Figure16.csv');
        p_tot_Chang20 = importdata('models/Chang1994voa2/Figure17.csv');
        theta = linspace(0,pi,2000).';
        phi = 0;
        v = R_o*[sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)];
        options = struct('d_vec', d_vec, ...
                         'omega', omega, ...
                         'R_i', R_i, ...
                         'R_o', R_o, ...
                         'P_inc', P_inc, ...
                         'E', E, ...
                         'nu', nu, ...
                         'rho_s', rho_s, ...
                         'rho_f', rho_f, ...
                         'c_f', c_f);
        
        data = e3Dss(v, options);
        p_tot = data(1).p + p_inc(v);
        
        figure(16)
        plot(theta*180/pi, real(p_tot(:,1)), p_tot_Chang15(:,1), p_tot_Chang15(:,2))
        title(sprintf('$$h/a = %.2f, k = %.1f, N_{eps} = %d$$', h/a, k(1), data(1).N_eps(1)),'interpreter','latex')
        xlabel('theta (degree)')
        ylabel('real part of pressure')
        ylim([-2 2])
        legend({'Present work', 'Reference Solution from Chang (1994) - Figure 16'})
        xlim([0 180])
        
        figure(17)
        plot(theta*180/pi, real(p_tot(:,2)), p_tot_Chang20(:,1), p_tot_Chang20(:,2))
        title(sprintf('$$h/a = %.2f, k = %.1f, N_{eps} = %d$$', h/a, k(2), data(1).N_eps(2)),'interpreter','latex')
        xlabel('theta (degree)')
        ylabel('real part of pressure')
        ylim([-2 2])
        legend({'Present work', 'Reference Solution from Chang (1994) - Figure 17'})
        xlim([0 180])
        
        
        folderName = 'results';
        if ~exist(folderName, 'dir')
            mkdir(folderName);
        end
        
        aspect = '0';
        elevation = 'm90';
        waveNumber = num2str(15);
        BC = 'SSBC';
        
        varCol.alpha_s = 0;
        varCol.beta_s = -90;
        scatteringCase = 'BI';
        model = 'Chang';
        varCol.scatteringCase = scatteringCase;
        
        
        
        varCol.f = k(1)*c_f(1)/(2*pi);
        saveName = [model '_' BC '_' scatteringCase '_A'  aspect '_E' elevation '_k' waveNumber '_Chang1'];
        varCol.saveName = saveName;
        filename = [folderName '/' saveName];
        printResultsToFile(filename, p_tot_Chang15(:,1), p_tot_Chang15(:,2), varCol, 1, 0, 'Chang', 'Results using WebPlotDigitizer to scan the results from Chang (1994)')
        
        
        varCol.f = k(2)*c_f(1)/(2*pi);
        waveNumber = num2str(20);
        saveName = [model '_' BC '_' scatteringCase '_A'  aspect '_E' elevation '_k' waveNumber '_Chang2'];
        varCol.saveName = saveName;
        filename = [folderName '/' saveName];
        printResultsToFile(filename, p_tot_Chang20(:,1), p_tot_Chang20(:,2), varCol, 1, 0, 'Chang', 'Results using WebPlotDigitizer to scan the results from Chang (1994)')
        
        
        varCol.f = k(1)*c_f(1)/(2*pi);
        waveNumber = num2str(15);
        saveName = [model '_' BC '_' scatteringCase '_A'  aspect '_E' elevation '_k' waveNumber '_1'];
        varCol.saveName = saveName;
        filename = [folderName '/' saveName];
        printResultsToFile(filename, theta*180/pi, real(p_tot(:,1)), varCol, 1, 0, 'NTNU_FFI', 'Analytic solution')
        
        varCol.f = k(2)*c_f(1)/(2*pi);
        waveNumber = num2str(20);
        saveName = [model '_' BC '_' scatteringCase '_A'  aspect '_E' elevation '_k' waveNumber '_2'];
        varCol.saveName = saveName;
        filename = [folderName '/' saveName];
        printResultsToFile(filename, theta*180/pi, real(p_tot(:,2)), varCol, 1, 0, 'NTNU_FFI', 'Analytic solution')
        
        saveName = [model '_' BC '_' scatteringCase '_A'  aspect '_E' elevation];
        filename = [folderName '/' saveName];
        figure(42)
        createConvergencePlot('2D',options,v,60,filename)
        
    case 'Ihlenburg (1998)'
        %% Ihlenburg (1998) example
        Ieye = [1, 0;
                0, 1;
                0, 0];
        ESBC = 0;
        for i = 1:3
            if i == 1
                nFreqs = 2000;
            elseif i == 2
                nFreqs = 5000;
            else
                nFreqs = 5000;
            end
            P_inc = 1; % Amplitude of incident wave
            rho_f = [1000, 1000]; % Density of outer fluid
            rho_s = 7669; % Density of solid
            c_f = [1524, 1524];   % Speed of sound in outer (exterior) fluid domain
            t = 0.15; % The thickness of the sphere
            R = 5; % Midsurface radius
            R_o = R+t/2; % Outer radius of shell
            R_i = R-t/2; % Inner radius of shell
            E = 207e9; % Youngs modulus of elastic material
            nu = 0.3; % Poisson ratio of elastic material

            SHBC = Ieye(i,1);
            SSBC = Ieye(i,2);
            defineBCstring
            if SSBC
                specialValues = [0.317054564603519 %
                                   0.392367003924396 %
                                   0.450625070843336 %
                                   0.495802583768421 %
                                   0.538194539627076 %
                                   0.584225864137528 %
                                   0.638340507928227 %
                                   0.703448659735311 %
                                   0.781278648263076 %
                                   0.872642563605135 %
                                   0.977700430550731 %
                                   1.096200471571476 %
                                   1.227664074628201 %
                                   1.371508743070356 %
                                   1.527120122444878 %
                                   1.693889192704899 %
                                   1.871228417747961];
            elseif ~SHBC
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
                               1.892808062465205]; %
            else
                specialValues = [];
            end
            
            k = linspace(2/nFreqs,2,nFreqs)'; % wave number
            omega = k*c_f(1);   % Wave number for outer fluid domain

            theta = 180*pi/180;
            d_vec = [1,0,0];
            options = struct('d_vec', d_vec, ...
                             'omega', omega, ...
                             'R_i', R_i, ...
                             'R_o', R_o, ...
                             'P_inc', P_inc, ...
                             'E', E, ...
                             'nu', nu, ...
                             'rho_s', rho_s, ...
                             'rho_f', rho_f, ...
                             'c_f', c_f,...
                             'calc_farField', 1);
                         
            if SSBC || ~SHBC
                if false
                    v = R_o*[cos(0),0,0];
                    f = @(k)-objFunc(k,options,v,c_f(1),0);
                    specialValues = findExtremas(f, 2/nFreqs, 2, 100000)';
                    v = R_o*[cos(pi),0,0];
                    f = @(k)-objFunc(k,options,v,c_f(1),0);
                    specialValues = [specialValues; findExtremas(f, 2/nFreqs, 2, 100000)'];
                    save(['results/Ihlenburg_' BC '_extremas'], 'specialValues')
                else
                    load(['results/Ihlenburg_' BC '_extremas'])
                end
                delta = 1e-4;
                specialValues = sort([specialValues; (specialValues-delta); (specialValues+delta)]);
            end
            k = unique(sort([k; specialValues]));
            omega = k*c_f(1);   % Wave number for outer fluid domain
            options.omega = omega;
            
            v = R_o*[cos(pi),0,0;
                     cos(0),0,0];
            data = e3Dss(v, options);

            figure(3)
            F = data(1).p;
            TS = 20*log10(abs(F));
            plot(k*R_o, TS(1,:))
            set(0,'defaulttextinterpreter','latex')
            hold on
            title('Ihlenburg (1998) example, $$\theta = 180^\circ$$')
            xlabel('$$kR_0$$')
            xlim([0, max(k*R_o)])
            ylabel('TS')   
            
            figure(4)
            F = data(1).p;
            TS = 20*log10(abs(F));
            plot(k*R_o, TS(2,:))
            set(0,'defaulttextinterpreter','latex')
            hold on
            title('Ihlenburg (1998) example - $$\theta = 0^\circ$$')
            xlabel('$$kR_0$$')
            xlim([0, max(k*R_o)])
            ylabel('TS')  
%             return
            folderName = 'results';
            if ~exist(folderName, 'dir')
                mkdir(folderName);
            end
%             
%             elevation = '0';
%             frequency = 'S';
%             
%             varCol.alpha_s = pi;
%             varCol.beta_s = 0;
%             scatteringCase = 'Sweep';
%             model = 'IL';
%             varCol.scatteringCase = scatteringCase;
%             
%             varCol.f_arr = omega/(2*pi);
%             
%             saveName = [model '_' BC '_' scatteringCase '_A180_E' elevation '_F' frequency];
%             varCol.saveName = saveName;
%             filename = [folderName '/' saveName];
%             printResultsToFile(filename, k*R_o, TS(1,:).', varCol, 1, 0, 'NTNU_FFI', 'Analytic solution')
% %             
%             saveName = [model '_' BC '_' scatteringCase '_A0_E' elevation '_F' frequency];
%             varCol.saveName = saveName;
%             filename = [folderName '/' saveName];
%             printResultsToFile(filename, k*R_o, TS(2,:).', varCol, 1, 0, 'NTNU_FFI', 'Analytic solution')
% %             
            if 0
                figure(40+i)
                nFreqs = 500;
                k = linspace(2/nFreqs,2,nFreqs)'; % wave number
                k = unique(sort([k; specialValues]));
                omega = k*c_f(1);   % Wave number for outer fluid domain
                options.omega = omega;

                createConvergencePlot('3D',options,v,35, ['results/IhlenburgError_' num2str(i)])
            end
        end
        
    case 'Fender (1972)'
        %% Fender (1972) example
        setFenderParameters
        
        a = R_o; % Outer radius of shell
        delta = t;
        b = a - delta; % Inner radius of shell
        nFreqs = 3000;
        k = linspace(32/nFreqs,32,nFreqs)';
        omega = k*c_f(1);
        
        d_vec = -[0,0,1].';

        %%%%%%%%%
        SPL_Fender0 = importdata('models/Fender1972sfa/Fig2.csv');
        SPL_Fender180 = importdata('models/Fender1972sfa/Fig3.csv');
        theta = 0;
        options = struct('d_vec', d_vec, ...
                         'omega', omega, ...
                         'R_i', b, ...
                         'R_o', a, ...
                         'P_inc', P_inc, ...
                         'E', E, ...
                         'nu', nu, ...
                         'rho_s', rho_s, ...
                         'rho_f', rho_f, ...
                         'c_f', c_f);
        if true
            v = [0,0,a*cos(0)];
            f = @(k)-objFunc(k,options,v,c_f(1),1);
%             specialValues = findExtremas(f, 2/nFreqs, 32, 100000)';
            v = [0,0,a*cos(pi)];
            f = @(k)-objFunc(k,options,v,c_f(1),1);
%             specialValues = [specialValues; findExtremas(f, 2/nFreqs, 32, 100000)'];
%             save('results/Fender_extremas', 'specialValues')
            load('results/Fender_extremas')
            delta = 1e-4;
            specialValues = sort([specialValues; (specialValues-delta); (specialValues+delta)]);
        end
        k = unique(sort([k; specialValues]));
        omega = k*c_f(1);   % Wave number for outer fluid domain
        options.omega = omega;
        v = [0,0,a*cos(0);
             0,0,a*cos(pi)];
        data = e3Dss(v, options);
        p_inc = @(v) P_inc*exp(1i*dot3(v,d_vec)*k.');
        p_tot = data(1).p + p_inc(v);
        SPL = 20*log10(abs(p_tot)); % sound pressure level
        figure(2)
        plot(k*a, SPL(1,:), SPL_Fender0(:,1), SPL_Fender0(:,2))
        set(0,'defaulttextinterpreter','latex')
        title('Predicted total pressure as a function of $$ka$$ at the surface of the shell, $$\theta = 0^\circ$$')
        xlabel('$$k_Oa$$')
        ylabel('Surface sound pressure level')
        ylim([-80 120])
        legend({'Present work', 'Reference Solution from Fender (1972) - Figure 2'})
        xlim(a*[k(1) k(end)])
        
        %%%%%%%%
        figure(3)
        plot(k*a, SPL(2,:), SPL_Fender180(:,1), SPL_Fender180(:,2))
        set(0,'defaulttextinterpreter','latex')
        title('Predicted total pressure as a function of $$ka$$ at the surface of the shell, $$\theta = 180^\circ$$')
        xlabel('$$k_Oa$$')
        ylabel('Surface sound pressure level')
        ylim([-60 120])
        legend({'Present work', 'Reference Solution from Fender (1972) - Figure 2'})
        xlim(a*[k(1) k(end)])


        
        folderName = 'results';
        if ~exist(folderName, 'dir')
            mkdir(folderName);
        end
        
        aspect = '0';
        elevation = '90';
        frequency = 'S';
        BC = 'NNBC';
        
        varCol.alpha_s = 0;
        varCol.beta_s = pi/2;
        scatteringCase = 'Sweep';
        model = 'Fender';
        varCol.scatteringCase = scatteringCase;
        
        varCol.f_arr = 2*pi*omega;
        
        saveName = [model '_' BC '_' scatteringCase '_A'  aspect '_E' elevation '_F' frequency '_Fender1'];
        varCol.saveName = saveName;
        filename = [folderName '/' saveName];
        printResultsToFile(filename, SPL_Fender0(:,1), SPL_Fender0(:,2), varCol, 1, 0, 'Fender', 'Results using WebPlotDigitizer to scan the results from Fender (1972)')
        
        
        saveName = [model '_' BC '_' scatteringCase '_A'  aspect '_E' elevation '_F' frequency '_Fender2'];
        varCol.saveName = saveName;
        filename = [folderName '/' saveName];
        printResultsToFile(filename, SPL_Fender180(:,1), SPL_Fender180(:,2), varCol, 1, 0, 'Fender', 'Results using WebPlotDigitizer to scan the results from Fender (1972)')
        
        
        saveName = [model '_' BC '_' scatteringCase '_A'  aspect '_E' elevation '_F' frequency '_1'];
        varCol.saveName = saveName;
        filename = [folderName '/' saveName];
        printResultsToFile(filename, k*a, SPL(1,:).', varCol, 1, 0, 'NTNU_FFI', 'Analytic solution')
        
        saveName = [model '_' BC '_' scatteringCase '_A'  aspect '_E' elevation '_F' frequency '_2'];
        varCol.saveName = saveName;
        filename = [folderName '/' saveName];
        printResultsToFile(filename, k*a, SPL(2,:).', varCol, 1, 0, 'NTNU_FFI', 'Analytic solution')

        figure(42)
        nFreqs = 500;
        k = linspace(2/nFreqs,32,nFreqs)'; % wave number
        k = unique(sort([k; specialValues]));
        omega = k*c_f(1);   % Wave number for outer fluid domain
        options.omega = omega;
        createConvergencePlot('3D',options,v,75, '../../../LaTeX/createFigures/contents/e3Dss/FenderError')
        savefig('results/FenderError.fig')
        
    case 'Bessel functions for large n'
        %% Create plot of bessel functinos for large n
        xi = 500;
%         xi = vpa(xi); % Symbolic precision was not used in the article
%         for this plot
%         digits(100)
        type = 1;
        switch type
            case 1
                N = 550;
            case 2
                N = 1050;
        end
        npts = N;
        B = zeros(npts,1);
        B1 = zeros(npts,1);
        B2 = zeros(npts,1);
        N_arr = linspace(1,N,npts).';
        parfor i = 1:length(N_arr)
            n = N_arr(i);
%             B(i) = bessel_s(n,xi,2).*bessel_s(n,xi,1);
        %     B(n) = bessel_s(n,xi,1);
            B1(i) = bessel_s(n,xi,1);
            B2(i) = bessel_s(n,xi,2);
        end
        semilogy(N_arr,abs(B1),N_arr,abs(B2))
        xlim([0,N])
%         semilogy(N_arr,abs(B))
        switch type
            case 1
                ylim([1e-10 1e4])
                printResultsToFile('results/besseljForLargeN1', N_arr, abs(B1), [], 0, 1)
                printResultsToFile('results/besselyForLargeN1', N_arr, abs(B2), [], 0, 1)
            case 2
                ylim([1e-220 1e220])
                printResultsToFile('results/besseljForLargeN2', N_arr, abs(B1), [], 0, 1)
                printResultsToFile('results/besselyForLargeN2', N_arr, abs(B2), [], 0, 1)
        end
    case 'Bessel functions upper bound'
        m = 0.953675026078925;
        b = -31.956879721039723;
        for ii = 100
            N = 1000;
            npts = N;
            X = zeros(npts,2);
            x_j = @(N) 1.0072*N+11.99;
%             x_y = @(N) m*N+b;
            x_y = @(N) 0.953675026078925*N-31.956879721039723;
            N_arr = linspace(100,N,npts).';
            for i = 1:length(N_arr)
                n = N_arr(i);
    %             
%                 y = @(x) bessel_s(n,x,1);
%                 f = @(x) abs(y(x))-10^-ii;
%                 X00 = x_j(n);
%     %             X0 = fminbnd(f,X00*[0.5, 1.5]);
%     %             X(i,1) = fminbnd(f,X00*0.5, X00*1.5);
%                 X(i,1) = fzero(f,X00);

                y = @(x) bessel_s(n,x,2);
                f = @(x) abs(y(x))-10^ii;
                X00 = x_y(n);
    %             X(i,2) = fminbnd(f,X00*0.5, X00*1.5);
                X(i,2) = fzero(f,X00);
    %             X(i,2) = bisection(f, X00*0.5, X00*1.5, 100, 1e-1);

            end
            [~,m,b] = regression(N_arr.',X(:,2).');
        end
        x_y = @(N) m*N+b;
        [~,m,b] = regression(N_arr.',X(:,1).');
%         [~,m,b] = regression(N_arr(500:end).',X(500:end,2).');
        x_j = @(N) m*N+b;
%         plot(N_arr,X(:,1),N_arr,X(:,2)) %, N_arr, x_j(N_arr), N_arr, x_y(N_arr))
%         plot(N_arr,X(:,2), N_arr, x_y(N_arr))
        plot(X(:,2), N_arr)
        xlabel('x')
        ylabel('N')
%         legend({'Plot of $$x$$ as a function of $$N$$ for when $$j_N(x) = 10^{-16}$$', 'Plot of $$x$$ as a function of $$N$$ for when $$y_N(x) = 10^{200}$$', 'Reference line for j', 'Reference line for y'}, 'interpreter','latex')
%         legend({'Plot of $$x$$ as a function of $$N$$ for when $$y_N(x) = 10^{200}$$', 'Reference line for y'}, 'interpreter','latex')
        
    case 'Plot in Paraview'
        ESBC = 0;
        SSBC = 0;
        for SHBC = 0 %[0, 1]
            for modelCell = {'IL'} %{'S5', 'S35', 'S135'}
                model = modelCell{1};
                switch model
                    case 'S5'
                        setS5Parameters
                    case 'S35'
                        setS35Parameters
                    case 'S135'
                        setS135Parameters
                    case 'IL'
                        setIhlenburgParameters
                end
                R_i = R_o - t;
                defineBCstring
                k_arr = [1.892808062465205, 1.8929];
                plotRealPart = 1;
                for k = 2 % 30/R_o(1)
                    omega = k*c_f(1);
                    % rho_f(end) = 0;
                    R_a = 2*5; %R_o(1);
                    d_vec = [1, 0, 0]';  
                    usePointChargeWave = false;
                    options = struct('d_vec', d_vec, ...
                                     'omega', omega, ...
                                     'R_i', R_i, ...
                                     'R_o', R_o, ...
                                     'P_inc', P_inc, ...
                                     'E', E, ...
                                     'nu', nu, ...
                                     'rho_s', rho_s, ...
                                     'rho_f', rho_f, ...
                                     'c_f', c_f, ...
                                     'SHBC', SHBC, ...
                                     'ESBC', ESBC, ...
                                     'SSBC', SSBC, ...
                                     'plotTimeOscillation', 0, ...
                                     'plotInTimeDomain', 0, ...
                                     'usePointChargeWave', usePointChargeWave, ...
                                     'usePlaneWave', ~usePointChargeWave, ...
                                     'R_a', R_a);
                    
%                     folderName = ['results/paraviewResults/' model '/'];
                    folderName = ['otherFunctions/e3Dss/results/paraviewResults/' model '/'];
                    if ~exist(folderName, 'dir')
                        mkdir(folderName);
                    end
                    vtfFileName = [folderName BC '_ka_' num2str(k*R_o(1))];

                    extraPts = 40;

                    createParaviewFiles_exact3(extraPts, vtfFileName, options)
                end
            end
        end

% %         solid = getSphericalShellData(5, 4.992,'Zaxis');
% %         solid = getSphericalShellData(3, 2.98,'Zaxis');
%         solid = getSphericalShellData(1, 0.95,'Zaxis');
%         [nodes, ~, visElements] = buildVisualization3dMesh_new3(solid.knots{1}, solid.knots{2}, solid.knots{3}, 200, 200, 0, solid);
%         VTKoptions = struct('name','otherFunctions/e3Dss/results/sources/S1_solid', 'celltype', 'VTK_HEXAHEDRON'); 
%         VTKdata.nodes = nodes;
%         VTKdata.visElements = visElements;
%         VTKdata.omega = 1;
%         makeVTKfile(VTKdata, VTKoptions);
    case 'Solid sphere'
        model = 'S5';
        ESBC = 1;
        SHBC = 0;
        SSBC = 0;
        setS5Parameters
        
        R_i = R_o - t;
        defineBCstring
        R_a = 1.5*R_o(1);
        r_s = 2*R_o(1);
        
        usePointChargeWave = 0;
        if usePointChargeWave
            d_vec = -[-sqrt(r_s^2-R_a^2), R_a, 0]';  
            d_vec = d_vec/norm(d_vec); 
        else
            d_vec = [1, 0, 0]';   
        end
        
        options = struct('d_vec', d_vec, ...
                         'R_i', R_i, ...
                         'R_o', R_o, ...
                         'P_inc', P_inc, ...
                         'E', E, ...
                         'nu', nu, ...
                         'rho_s', rho_s, ...
                         'rho_f', rho_f, ...
                         'c_f', c_f, ...
                         'usePointChargeWave', usePointChargeWave, ...
                         'usePlaneWave', ~usePointChargeWave, ...
                         'plotTimeOscillation', 0, ...
                         'plotInTimeDomain', 1, ...
                         'SHBC', SHBC, ...
                         'ESBC', ESBC, ...
                         'SSBC', SSBC, ...
                         'R_a', R_a, ...
                         'r_s', r_s, ...
                         'f_c', 1500, ...
                         'N', 2^9, ...
                         'computeForSolidDomain', 1);

        options.Eps = 1e-8;
        extraPts = 30;
        folderName = ['results/paraviewResults/' model '_' BC '/'];
        if ~exist(folderName, 'dir')
            mkdir(folderName);
        end

        createParaviewFiles_exact3(extraPts, folderName, options)  
    case 'Solid sphere inside spherical shell'
        model = 'S15';
        ESBC = 1;
        SHBC = 0;
        SSBC = 0;
        setS15Parameters
        
        R_i = R_o - t;
        defineBCstring
        R_a = 1.5*R_o(1);
        r_s = 2*R_o(1);
        usePointChargeWave = 1;
        if usePointChargeWave
            d_vec = -[-sqrt(r_s^2-R_a^2), R_a, 0]';  
            d_vec = d_vec/norm(d_vec); 
        else
            d_vec = [1, 0, 0]';   
        end
        
        options = struct('d_vec', d_vec, ...
                         'R_i', R_i, ...
                         'R_o', R_o, ...
                         'P_inc', P_inc, ...
                         'E', E, ...
                         'nu', nu, ...
                         'rho_s', rho_s, ...
                         'rho_f', rho_f, ...
                         'c_f', c_f, ...
                         'usePointChargeWave', usePointChargeWave, ...
                         'usePlaneWave', ~usePointChargeWave, ...
                         'plotTimeOscillation', 0, ...
                         'plotInTimeDomain', 1, ...
                         'SHBC', SHBC, ...
                         'ESBC', ESBC, ...
                         'SSBC', SSBC, ...
                         'R_a', R_a, ...
                         'r_s', r_s, ...
                         'f_c', 1500, ...
                         'N', 2^9);

        options.Eps = 1e-8;
        extraPts = 40;
        folderName = ['results/paraviewResults/' model '_' BC '/'];
        if ~exist(folderName, 'dir')
            mkdir(folderName);
        end

        createParaviewFiles_exact3(extraPts, folderName, options)        
    case 'create spy matrix'
        setS135Parameters
        SHBC = false;
        ESBC = false;
        SSBC = false;
        R_i = R_o - t;
        f = 30e3;
        omega = 2*pi*f; % Angular frequency
        
        defineBCstring
        d_vec = [0, 0, -1];
        
        options = struct('d_vec', d_vec,... 
                         'omega', omega, ...
                         'R_i', R_i, ...
                         'R_o', R_o, ...
                         'P_inc', P_inc, ...
                         'E', E, ...
                         'nu', nu, ...
                         'rho_s', rho_s, ...
                         'rho_f', rho_f, ...
                         'c_f', c_f);
        v = -R_o(1)*d_vec;
        data = e3Dss(v, options);
        %% Create spy matrix (requires to be in debug mode in getCoeffs.m)
        % uncomment keyboard command in getCoeffs.m to get spy matrix
        % Remember to use model = 'S135' and scatteringCase = 'BI' and f = 30e3
        % to reproduce plot
    case 'time domain solution - 1D Visualization'
        if 0
            f_c = 300; % (300) source pulse center freq.
            N = 2^9;
            M = 2*N; % corresponds to fs = sampling rate
            T = 30/f_c;
            B = N/T; % bandwidth
            f_L = -B/2;
            f_R = B/2;
            df = 1/T;
            type = 1;
            useNegativeFreqsInScatPres = 0;
            if useNegativeFreqsInScatPres
                f = linspace(f_L,f_R-df,N);
            else
                f = linspace(0,f_R-df,N/2);
            end
            omega = 2*pi*f;

            d_vec = [0,0,1].';
            R_o = 2;
            omega_c = 2*pi*f_c;
            P_inc = @(omega) P_inc_(omega, omega_c,type);
            c_f = 1500;
            k_c = omega_c/c_f;
            npts = 400;
            z = linspace(-40,-R_o,npts).';

            %%%%%%%%%%%%%%%%%%%%%%%%%%
            % Plots
            figure(1); clf;
            t = linspace(0,30/f_c,1000);
            subplot(311), plot(t,Pt_inc_(t,0,omega_c,k_c,type))
            Sw1 = P_inc_(omega,omega_c,type);
            subplot(312), plot(f,abs(Sw1))
            subplot(313), plot(f,atan2(imag(Sw1),real(Sw1)));
            drawnow
    %         return
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            % Plot dynamic time domain solution in 1D
            vv = [zeros(npts,1), zeros(npts,1), z];
            options = struct('d_vec', d_vec, ...
                             'omega', omega(2:end), ...
                             'R_o', R_o, ...
                             'P_inc', P_inc, ...
                             'c_f', c_f);
            data = e3Dss(vv, options);
        else
            f_c = 1500; % (300) source pulse center freq.
            N = 2^11;
            M = N; % corresponds to fs = sampling rate
            T = 240/f_c; %60/f_c
            B = N/T; % bandwidth
            f_L = -B/2;
            f_R = B/2;
            df = 1/T;
            f = linspace(0,f_R-df,N/2);
            omega = 2*pi*f;
            type = 1;
            d_vec = [0,0,1].';
            omega_c = 2*pi*f_c;
            c_f = 1500;
            k_c = omega_c/c_f;

            %%%%%%%%%%%%%%%%%%%%%%%%%%
            % Plots
            figure(1); clf;
%             t = linspace(-1/f_c,6/f_c,1000);
            tt = linspace(0,1/f_c,1000);
            tt = [-1e-3,tt,3e-3];
            Pt_inc = Pt_inc_(tt,0,omega_c,k_c,type);
            subplot(311), plot(tt,Pt_inc)
            ft = linspace(0,f_R,2000);
            omegat = 2*pi*ft;
            P_inc = P_inc_(omegat,omega_c,type);
            subplot(312), plot(ft,abs(P_inc))
            subplot(313), plot(ft,atan2(imag(P_inc),real(P_inc)));
            drawnow
%             printResultsToFile('results/Pt_inc', tt.', Pt_inc.', [], 0, 1)
%             printResultsToFile('results/P_inc', omegat.', abs(P_inc).', [], 0, 1)
%             return
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            % Plot dynamic time domain solution in 1D
%             setS15Parameters
            setS5Parameters
            npts = 1000;
            R_a = 1.5*R_o(1);
%             z = linspace(-1.3*R_a,-R_o(1),npts).';
            z = linspace(R_o(1),20,npts).';
            R_i = R_o - t;
            ESBC = 1;
            vv = [zeros(npts,1), zeros(npts,1), z];
            rho_f = rho_f(1:end-1);
            c_f = c_f(1:end-1);
            R_i = R_i(1:end-1);
            P_inc = @(omega) P_inc_(omega, omega_c,type);
            
            options = struct('d_vec', d_vec, ...
                             'omega', omega(2:end), ...
                             'R_i', R_i, ...
                             'R_o', R_o, ...
                             'P_inc', P_inc, ...
                             'E', E, ...
                             'nu', nu, ...
                             'rho_s', rho_s, ...
                             'rho_f', rho_f, ...
                             'c_f', c_f);
            data = e3Dss(vv, options);
        end
        startIdx = 2000;
        totField = zeros(npts,N/2);
        PincField = zeros(npts,N/2);
        
        for n = 0:N-1
            f = f_L + (f_R-f_L)/N*n;
            omega = 2*pi*f;
            k = omega/c_f(1);
            k_vec = d_vec*k;
            if n >= N/2+1
                PincField(:,n-N/2+1) = P_inc_(omega,omega_c,type).*exp(1i*dot3(vv, k_vec));
                totField(:,n-N/2+1) = PincField(:,n-N/2+1) + data(1).p(:,n-N/2);
            end
        end
        dt = T/M;
        totFieldTime = 2/T*fft(totField,M,2);
        PincFieldTime = 2/T*fft(PincField,M,2);
        
        temp = totFieldTime;
        totFieldTime(:,1:M-startIdx+1) = temp(:,startIdx:end);
        totFieldTime(:,M-startIdx+2:end) = temp(:,1:startIdx-1);
        temp = PincFieldTime;
        PincFieldTime(:,1:M-startIdx+1) = temp(:,startIdx:end);
        PincFieldTime(:,M-startIdx+2:end) = temp(:,1:startIdx-1);
        
        figure(2)
        m_arr = 0:M-1;
        for m = m_arr
            t = dt*m;
            plot(z,real(totFieldTime(:,m+1)),z,real(PincFieldTime(:,m+1)));
            titleStr = sprintf('Time step %4d, t = %5.3fs. T = %5.3fs.',m,t,T);
            title(titleStr)
            ylim([-1.3,1.3])
            legend('p_{tot}', 'p_{inc}')
            drawnow
            pause(0.1)
            if m == 120
                pause
            end
        end
        figure(3)
        plot(dt*m_arr,real(totFieldTime(1,:)),dt*m_arr,real(PincFieldTime(1,:)));
        xlim([0, 120*dt])
        
        
        
end