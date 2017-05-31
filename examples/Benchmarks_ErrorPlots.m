close all
clear all %#ok

addpath ..
addpath ../utils
addpath ../models

startMatlabPool

%% Calculate errors
useSymbolicPrecision = 0;
if useSymbolicPrecision
    Eps = 1e-40;
else
    Eps = eps;
end

models = {'S1','S3','S5','S35','S15','S13','S135'};

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
parfor i = 1:length(models)
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
                    f = 10.^linspace(log10(1e-3/C_eps),log10(4e2/C_eps),nFreqs);
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
                    filename = ['../results/errors_' model '_' BC '_Symbolic' num2str(useSymbolicPrecision)];

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