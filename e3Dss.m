function [layer,N_eps,flag] = e3Dss(newLayers, newOptions)

% The function e3Dss (exact 3D scattering solutions) computes the solution
% to scattering problems on multilayered spherical shells impinged by a 
% plane wave, a wave due to a point source or a radially pulsating wave.
%
% The function handles the case in which both X and omega are vectors.
% 
% The following improvement will be implemented in the future:
%   - Implement a scaling in the linear system of equations such that
%   evaluations of the spherical bessel functions j_n(x) and y_n(x) are
%   replaced by j_n(x)*y_{n-1}(x) and y_n(x)*j_{n+1}(x), respectively. This will in
%   turn solve the problem of evaluations extending the range of the given
%   precision (10^290 is implemented for double), and thus enable
%   evaluations at even higher frequencies.
%
% Author: Jon Vegard Venås
% E-mail: JonVegard.Venas@sintef.no
% Institute: SINTEF Digital
% Release: 2
% Release date: 21/03/2020

options = struct('d_vec',   [0;0;1],  ... 	% Direction of the incident wave
                'omega',   2*pi*1e3, ...    % Angular frequency
                'P_inc',   1,    ...     	% Amplitude of incident wave
                'Eps',     eps,  ...        % Small parameter for series truncation
                'N_max',   inf,  ...        % Upper limit for the number of terms in the series
                'prec',    'double', ...    % Precision of the calculations
                'Display', 'final', ...     % Print status during computations
                'BC',      'SHBC', ...      % Boundary condition on the inermost layer 'SSBC' (Sound soft boundary condition), 'NNBC' (Neumann-Neumann boundary condition) 
                'applyLoad', 'planeWave');  % Incident wave type: I.e. planeWave, pointCharge, mechExcitation, surfExcitation, radialPulsation
if nargin > 1
    options = updateOptions(options,newOptions);
end
layer = getDefaultParameters(newLayers);

layer = updateOptions(layer,newLayers);

% Supress warning for ill-conditioned matrices in getCoeffs.m
poolobj = gcp('nocreate');
if isempty(poolobj)
    warning('off', 'MATLAB:nearlySingularMatrix')
    warning('off', 'MATLAB:singularMatrix')
else
    pctRunOnAll warning('off', 'MATLAB:singularMatrix')
    pctRunOnAll warning('off', 'MATLAB:nearlySingularMatrix')
end
prec = options.prec;

omega = options.omega;
if isrow(omega)
    options.omega = omega.';
end
if strcmp(layer{1}.media,'solid') % the outermost unbounded domain is solid
    error('This case is not implemented')
end
M = numel(layer);
for m = 1:M
    calc_errPresCond = layer{m}.calc_errPresCond;
    calc_errDispCond = layer{m}.calc_errDispCond;
    switch layer{m}.media
        case 'fluid'
            calc_errHelm = layer{m}.calc_errHelm;
            calc_dp = or(layer{m}.calc_dp, calc_errDispCond);

            calc_p_laplace = layer{m}.calc_p_laplace;

            layer{m}.calc_p = layer{m}.calc_p || calc_errPresCond;
            layer{m}.calc_dpdr = any(calc_dp) || calc_p_laplace;
            layer{m}.calc_dp = calc_dp;
            layer{m}.calc_p_laplace = calc_p_laplace || calc_errHelm;
            layer{m}.calc_dpdt = layer{m}.calc_dpdr;
            layer{m}.calc_d2pdr2 = calc_p_laplace || calc_errHelm;
            layer{m}.calc_d2pdt2 = calc_p_laplace || calc_errHelm;
            layer{m}.calc_errors = calc_errDispCond || calc_errPresCond || calc_errHelm;
            layer{m}.calc_farFieldOnly = ~any([layer{m}.calc_p, layer{m}.calc_dpdr, layer{m}.calc_dpdt, layer{m}.calc_d2pdr2, layer{m}.calc_d2pdt2]);
        case {'solid','viscoelastic'}
            
            calc_u = layer{m}.calc_u;
            calcCartesianDispDerivatives = any(layer{m}.calc_du(:));
            calc_errNav = layer{m}.calc_errNav;
            layer{m}.calc_u = or(calc_u, calc_errDispCond);
            layer{m}.calc_u_r = calc_errDispCond || calcCartesianDispDerivatives || any(calc_u) || calc_errNav;
            layer{m}.calc_u_t = calc_errDispCond || calcCartesianDispDerivatives || any(calc_u) || calc_errNav;
            layer{m}.calc_du_rdr = calcCartesianDispDerivatives;
            layer{m}.calc_du_rdt = calcCartesianDispDerivatives;
            layer{m}.calc_du_tdr = calcCartesianDispDerivatives;
            layer{m}.calc_du_tdt = calcCartesianDispDerivatives;
            calc_sigma = layer{m}.calc_sigma;
            calc_sigma_s = layer{m}.calc_sigma_s;
            if any(calc_sigma_s(2:end))
                error('Not yet implemented')
            end
            calcStresses = any(calc_sigma) || any(calc_sigma_s) || calc_errNav || calc_errPresCond;
            layer{m}.calc_sigma_rr = calcStresses;
            layer{m}.calc_sigma_tt = calcStresses;
            layer{m}.calc_sigma_pp = calcStresses;
            layer{m}.calc_sigma_rt = calcStresses;
            layer{m}.calcStresses = calcStresses;
            layer{m}.calc_navier1 = calc_errNav;
            layer{m}.calc_navier2 = calc_errNav;
            layer{m}.calc_sigma = or(calc_sigma, calcStresses);
            layer{m}.calcCartesianDispDerivatives = calcCartesianDispDerivatives;
            layer{m}.calc_errors = calc_errDispCond || calc_errPresCond || calc_errNav;
    end
    if isfield(layer{m},'X')
        %% Coordinate transformation
        % Due to symmetry, we can do a coordinate transformation such that we only
        % need to compute the solution for the special case k_vec = k*[0, 0, 1].
        [r, theta, phi, A] = coordTransform(layer{m}.X, options.d_vec);
        layer{m}.r = r.';
        layer{m}.theta = theta.';
        layer{m}.phi = phi.';
    end
end
options.SHBC = strcmp(options.BC,'SHBC');
options.SSBC = strcmp(options.BC,'SSBC');
options.IBC = strcmp(options.BC,'IBC');
if options.IBC && ~isfield(options,'z')
    options.z = 1;
end

if any(omega == 0) && ( strcmp(options.applyLoad, 'mechExcitation') || ...
                       (strcmp(options.applyLoad, 'pointCharge') && options.r_s ~= 0) || ...
                       (strcmp(options.applyLoad, 'surfExcitation') && abs(options.theta_s - pi) > options.Eps))
    warning('This case has no unique solution for omega=0, these values will be removed from the computation.')
    omega(omega == 0) = [];
    options.omega = omega;
end
if any([options.SHBC,options.SSBC,options.IBC]) && layer{end}.R_i == 0
    error('Boundary conditions may not be used at the origin')
end
switch options.applyLoad
    case {'planeWave','radialPulsation'}
        m_s = 1;
    case {'pointCharge','mechExcitation','surfExcitation','custom'}
        if ~isfield(options,'r_s')
            if strcmp(options.applyLoad, 'pointCharge')
                options.r_s = 2*layer{1}.R_i;
            elseif strcmp(options.applyLoad, 'surfExcitation') || strcmp(options.applyLoad, 'custom')
                options.r_s = layer{1}.R_i;
            end
        end
        if strcmp(options.applyLoad, 'surfExcitation') && ~isfield(options,'theta_s')
            options.theta_s = [0,pi];
        end
        if strcmp(options.applyLoad, 'custom') && ~isfield(options,'theta_s')
            options.theta_s = 1;
        end
        r_s = options.r_s;
        m_s = 1;
        for m = 1:M
            if r_s < layer{m}.R_i
                m_s = m_s + 1;
            else
                break
            end
        end
end
options.m_s = m_s;

%% Compute the solution with d_vec = [0, 0, 1]
[layer,N_eps,flag] = e3Dss_0(layer, options);

%% Coordinate transformation (Transform back to original Cartesian coordinate system)

% Allocate memory
nFreqs = length(omega);

for m = 1:M
    if isfield(layer{m},'X')
        Theta = repmat(layer{m}.theta,nFreqs,1);
        Phi = repmat(layer{m}.phi,nFreqs,1);
        R = repmat(layer{m}.r,nFreqs,1);
        indices = logical(R < options.Eps);
        n_X = size(layer{m}.X,1);
        switch layer{m}.media
            case 'fluid'
                if any(layer{m}.calc_dp) || layer{m}.calc_p_laplace
                    dpdr = layer{m}.dpdr;
                    dpdt = sin(Theta).*layer{m}.dpdt; % rescale dpdt
                    dpdX_m = cell(3,1);

                    dpdX_m{1} = dpdr.*sin(Theta).*cos(Phi) + dpdt.*cos(Theta).*cos(Phi)./R;
                    dpdX_m{1}(indices) = 0;

                    dpdX_m{2} = dpdr.*sin(Theta).*sin(Phi) + dpdt.*cos(Theta).*sin(Phi)./R;
                    dpdX_m{2}(indices) = 0;

                    dpdX_m{3} = dpdr.*cos(Theta) - dpdt.*sin(Theta)./R;
                    dpdX_m{3}(indices) = dpdr(indices);

                    dpdX = cell(3,1);
                    dpdX{1} = zeros(n_X,nFreqs);
                    dpdX{2} = zeros(n_X,nFreqs);
                    dpdX{3} = zeros(n_X,nFreqs);
                    for ii = 1:3
                        for jj = 1:3
                            dpdX{ii} = dpdX{ii} + A(ii,jj)*dpdX_m{jj}.';
                        end
                    end
                    if layer{m}.calc_dp(1)
                        layer{m}.dpdx = dpdX{1};
                    end
                    if layer{m}.calc_dp(2)
                        layer{m}.dpdy = dpdX{2};
                    end
                    if layer{m}.calc_dp(3)
                        layer{m}.dpdz = dpdX{3};
                    end
                    dpdt = layer{m}.dpdt; % use scaled version of dpdt for the laplace operator
                    if layer{m}.calc_p_laplace
                        d2pdr2 = layer{m}.d2pdr2;
                        d2pdt2 = layer{m}.d2pdt2;

                        temp = 2./R.*dpdr + 1./R.^2.*cos(Theta).*dpdt + 1./R.^2.*d2pdt2;                    
                        temp(indices) = 0;
                        layer{m}.p_laplace = d2pdr2 + temp;
                        layer{m}.p_laplace = layer{m}.p_laplace.';
                    end
                end
                if layer{m}.calc_p
                    layer{m}.p = layer{m}.p.';
                end
                if layer{m}.calc_p_0
                    layer{m}.p_0 = layer{m}.p_0.';
                end
            case {'solid','viscoelastic'}
                if layer{m}.calcStresses
                    % Transform the stresses in the spherical coordinate system to the 
                    % Cartesian coordinate system
                    sigma_m = cell(6,1);
                    sigma_m{1} = layer{m}.sigma_rr;
                    sigma_m{2} = layer{m}.sigma_tt;
                    sigma_m{3} = layer{m}.sigma_pp;
                    sigma_m{4} = zeros(nFreqs,n_X,prec);
                    sigma_m{5} = zeros(nFreqs,n_X,prec);
                    sigma_m{6} = layer{m}.sigma_rt;

                    D = getStressTransformationMatrix(layer{m}.theta,layer{m}.phi,2);
                    sigma_X_m = cell(6,1);
                    for ii = 1:6
                        sigma_X_m{ii} = zeros(nFreqs,n_X,prec);
                    end
                    sigma_X = sigma_X_m;
                    for ii = 1:6
                        for jj = 1:6
                            D_kl = D(ii, jj, :);
                            D_kl = D_kl(:);
                            sigma_X_m{ii} = sigma_X_m{ii} + repmat(D_kl.',nFreqs,1).*sigma_m{jj};
                        end
                    end
                    if m == M
                        sigma_X_m{1}(indices) = sigma_m{1}(indices);
                        sigma_X_m{2}(indices) = sigma_m{2}(indices);
                        sigma_X_m{3}(indices) = sigma_m{3}(indices);
                        sigma_X_m{4}(indices) = 0;
                        sigma_X_m{5}(indices) = 0;
                        sigma_X_m{6}(indices) = 0;
                    end

                    alpha = zeros(3,prec);
                    I = eye(3);
                    for ii = 1:3
                        for jj = 1:3
                            alpha(ii,jj) = dot(I(:,ii), A(:,jj));
                        end
                    end

                    vgt = [1 1;
                           2 2;
                           3 3;
                           2 3;
                           1 3;
                           1 2];
                    vgtinv = [1 6 5;
                              6 2 4;
                              5 4 3];
                    for vgtIdx = 1:6
                        for ii = 1:3
                            for jj = 1:3
                                sigma_X{vgtIdx} = sigma_X{vgtIdx} + alpha(vgt(vgtIdx,1),ii)*alpha(vgt(vgtIdx,2),jj)*sigma_X_m{vgtinv(ii,jj)};
                            end
                        end
                    end
                    for ii = 1:6
                        sigma_X{ii} = sigma_X{ii}.';
                    end
                    if layer{m}.calc_sigma(1)
                        layer{m}.sigma_xx = sigma_X{1};
                    end
                    if layer{m}.calc_sigma(2)
                        layer{m}.sigma_yy = sigma_X{2};
                    end
                    if layer{m}.calc_sigma(3)
                        layer{m}.sigma_zz = sigma_X{3};
                    end
                    if layer{m}.calc_sigma(4)
                        layer{m}.sigma_yz = sigma_X{4};
                    end
                    if layer{m}.calc_sigma(5)
                        layer{m}.sigma_xz = sigma_X{5};
                    end
                    if layer{m}.calc_sigma(6)
                        layer{m}.sigma_xy = sigma_X{6};
                    end
                    if layer{m}.calc_sigma_s(1)
                        layer{m}.sigma_rr = layer{m}.sigma_rr.';
                    end
                end
                if layer{m}.calcCartesianDispDerivatives
                    % Transform the derivatives in the spherical coordinate system to the 
                    % Cartesian coordinate system
                    u_r = layer{m}.u_r;
                    u_t = layer{m}.u_t;

                    du_m = cell(3,3);
                    du_m{1,1} = layer{m}.du_rdr;
                    du_m{1,2} = layer{m}.du_rdt.*sin(Theta); % rescale du_rdt
                    du_m{1,3} = zeros(nFreqs,n_X,prec);
                    du_m{2,1} = layer{m}.du_tdr.*sin(Theta); % rescale du_tdr
                    du_m{2,2} = layer{m}.du_tdt;
                    du_m{2,3} = zeros(nFreqs,n_X,prec);
                    du_m{3,1} = zeros(nFreqs,n_X,prec);
                    du_m{3,2} = zeros(nFreqs,n_X,prec);
                    du_m{3,3} = zeros(nFreqs,n_X,prec);

                    [J_e,~,~,J_si,J_1,J_2] = getDerivativeTransformationMatrices(layer{m}.theta,layer{m}.phi,layer{m}.r);
                    du_X_m = cell(3,3);
                    for ii = 1:3
                        for jj = 1:3
                            du_X_m{ii,jj} = zeros(nFreqs,n_X,prec);
                        end
                    end
                    du_X = du_X_m;
                    temp = du_X_m;
                    for ii = 1:3
                        for jj = 1:3
                            for ll = 1:3
                                J_si_lj = J_si(ll, jj, :);
                                J_si_lj = J_si_lj(:);
                                temp{ii,jj} = temp{ii,jj} + du_m{ii,ll}.*repmat(J_si_lj.',nFreqs,1);
                            end
                        end
                    end
                    for ii = 1:3
                        for jj = 1:3
                            for ll = 1:3
                                J_e_lj = J_e(ll, ii, :); % transpose of J_e ...
                                J_e_lj = J_e_lj(:);
                                du_X_m{ii,jj} = du_X_m{ii,jj} + repmat(J_e_lj.',nFreqs,1).*temp{ll,jj};
                            end
                        end
                    end
                    for ii = 1:3
                        for jj = 1:3
                            J_1_ij = J_1(ii, jj, :);
                            J_1_ij = J_1_ij(:);
                            du_X_m{ii,jj} = du_X_m{ii,jj} + repmat(J_1_ij.',nFreqs,1).*u_r;

                            J_2_ij = J_2(ii, jj, :);
                            J_2_ij = J_2_ij(:);
                            du_X_m{ii,jj} = du_X_m{ii,jj} + repmat(J_2_ij.',nFreqs,1).*u_t;
                        end
                    end
                    if m == M
                        du_X_m{1,1}(indices) = layer{m}.du_rdr(indices);
                        du_X_m{2,2}(indices) = layer{m}.du_rdt(indices);
                        du_X_m{3,3}(indices) = layer{m}.du_tdr(indices);
                        for ii = 1:3
                            for jj = 1:3
                                if ii ~= jj
                                    du_X_m{ii,jj}(indices) = 0;
                                end
                            end
                        end
                    end

                    for ii = 1:3
                        for jj = 1:3
                            temp{ii,jj}(:) = 0;
                        end
                    end
                    Ainv = inv(A);
                    for ii = 1:3
                        for jj = 1:3
                            for ll = 1:3
                                temp{ii,jj} = temp{ii,jj} + A(ii, ll)*du_X_m{ll,jj};
                            end
                        end
                    end
                    for ii = 1:3
                        for jj = 1:3
                            for ll = 1:3
                                du_X{ii,jj} = du_X{ii,jj} + temp{ii, ll}*Ainv(ll,jj);
                            end
                        end
                    end
                    for ii = 1:3
                        for jj = 1:3
                            du_X{ii,jj} = du_X{ii,jj}.';
                        end
                    end

                    if layer{m}.calc_du(1,1)
                        layer{m}.du_xdx = du_X{1,1};
                    end
                    if layer{m}.calc_du(1,2)
                        layer{m}.du_xdy = du_X{1,2};
                    end
                    if layer{m}.calc_du(1,3)
                        layer{m}.du_xdz = du_X{1,3};
                    end
                    if layer{m}.calc_du(2,1)
                        layer{m}.du_ydx = du_X{2,1};
                    end
                    if layer{m}.calc_du(2,2)
                        layer{m}.du_ydy = du_X{2,2};
                    end
                    if layer{m}.calc_du(2,3)
                        layer{m}.du_ydz = du_X{2,3};
                    end
                    if layer{m}.calc_du(3,1)
                        layer{m}.du_zdx = du_X{3,1};
                    end
                    if layer{m}.calc_du(3,2)
                        layer{m}.du_zdy = du_X{3,2};
                    end
                    if layer{m}.calc_du(3,3)
                        layer{m}.du_zdz = du_X{3,3};
                    end
                    if strcmp(layer{m}.media,'viscoelastic') && layer{m}.calc_p
                        Omega = repmat(omega,n_X,1);
                        E = layer{m}.E;
                        nu = layer{m}.nu;
                        rho = layer{m}.rho;
                        % Compute derived quantities
                        K = E/(3*(1-2*nu));
                        G = E/(2*(1+nu));

                        c_s_1 = sqrt((3*K+4*G)/(3*rho)); % longitudinal wave velocity 
                        c_s_2 = sqrt(G/rho); % shear wave velocity
                        cL1 = c_s_1;
                        cT1 = c_s_2; 
                        kT = wv/cT1;
                        kL = wv./cL1;
                        muv = 1i*rho*wv/kT^2;
                        c = repmat(4*1i*wv*muv/(3*rho) + ev.^2./kL.^2,n_X,1);

                        layer{m}.p = -1i*rho*c^2./Omega.*(du_X{1,1}+du_X{2,2}+du_X{3,3});
                    end
                end
                if any(layer{m}.calc_u)
                    u_X_m = cell(3,1);
                    u_r = layer{m}.u_r;
                    u_t = layer{m}.u_t.*sin(Theta); % rescale u_t
                    u_X_m{1} = u_r.*sin(Theta).*cos(Phi) + u_t.*cos(Theta).*cos(Phi);
                    u_X_m{2} = u_r.*sin(Theta).*sin(Phi) + u_t.*cos(Theta).*sin(Phi);
                    u_X_m{3} = u_r.*cos(Theta) - u_t.*sin(Theta);
                    if m == M
                        u_X_m{1}(indices) = 0;
                        u_X_m{2}(indices) = 0;
                        u_X_m{3}(indices) = u_r(indices);
                    end
                    u_X = cell(3,1);
                    u_X{1} = zeros(n_X,nFreqs,prec);
                    u_X{2} = zeros(n_X,nFreqs,prec);
                    u_X{3} = zeros(n_X,nFreqs,prec);

                    for ii = 1:3
                        for jj = 1:3
                            u_X{ii} = u_X{ii} + A(ii,jj)*u_X_m{jj}.';
                        end
                    end
                    if layer{m}.calc_u(1)
                        layer{m}.u_x = u_X{1};
                    end
                    if layer{m}.calc_u(2)
                        layer{m}.u_y = u_X{2};
                    end
                    if layer{m}.calc_u(3)
                        layer{m}.u_z = u_X{3};
                    end
                end
                if layer{m}.calc_errNav
                    layer{m}.navier1 = layer{m}.navier1.';
                    layer{m}.navier2 = layer{m}.navier2.';
                end
        end
    end
end
SHBC = options.SHBC;
SSBC = options.SSBC;
if M > 1
    nextMedia = layer{2}.media;
else
    nextMedia = 'void';
end
for m = 1:M
    if layer{m}.calc_errors
        Eps = options.Eps;
        R_i = layer{m}.R_i;
        isSphere = R_i == 0;
        n_X = size(layer{m}.X,1);
        X = layer{m}.X;
        switch layer{m}.media
            case 'fluid'
                if layer{m}.calc_errHelm
                    p_laplace = layer{m}.p_laplace;
                    k = omega/layer{m}.c_f;
                    K = repmat(k,n_X,1);

                    layer{m}.err_helmholtz = max(abs(p_laplace+K.^2.*layer{m}.p),[],1)./max(abs(K.^2.*layer{m}.p),[],1);
                end
            case {'solid','viscoelastic'}
                if layer{m}.calc_errNav
                    Theta = repmat(layer{m}.theta,nFreqs,1);
                    rho = layer{m}.rho;
                    u_r = layer{m}.u_r;
                    u_t = layer{m}.u_t.*sin(Theta); % rescale u_t
                    u_r = u_r.';
                    u_t = u_t.';
                    Omega = repmat(omega,n_X,1);
                    layer{m}.err_navier1 = max(abs(layer{m}.navier1 + rho*Omega.^2.*u_r),[],1)./max(abs(rho*Omega.^2.*u_r),[],1);
                    layer{m}.err_navier2 = max(abs(layer{m}.navier2 + rho*Omega.^2.*u_t),[],1)./max(abs(rho*Omega.^2.*u_t),[],1);
                end
        end
        P_inc = options.P_inc;
        if layer{m}.calc_errPresCond && ~isSphere
            if m < M
                X2 = layer{m+1}.X;
                [indices, indices2] = findMatchingPoints(X,X2,Eps);
            else
                indices = abs(norm2(X) - R_i) < 10*Eps;
            end
            X_s = X(indices,:);
            switch layer{m}.media
                case 'fluid'
                    p_tot = layer{m}.p(indices,:);
                    if m == m_s
                        k = omega/layer{m}.c_f;
                        switch options.applyLoad
                            case 'planeWave'
                                k_vec = options.d_vec*k;
                                p_inc = P_inc*exp(1i*X_s*k_vec);
                                p_tot = p_tot + p_inc;
                            otherwise
                                error('Not implemented')
                        end
                    end
                    switch nextMedia
                        case {'solid','viscoelastic'}
                            sigma_cart = cell(6,1);
                            sigma_cart{1} = layer{m+1}.sigma_xx(indices2,:);
                            sigma_cart{2} = layer{m+1}.sigma_yy(indices2,:);
                            sigma_cart{3} = layer{m+1}.sigma_zz(indices2,:);
                            sigma_cart{4} = layer{m+1}.sigma_yz(indices2,:);
                            sigma_cart{5} = layer{m+1}.sigma_xz(indices2,:);
                            sigma_cart{6} = layer{m+1}.sigma_xy(indices2,:);

                            sigma_rr = zeros(size(sigma_cart{1}),prec);

                            phi_temp = atan2(X_s(:,2),X_s(:,1));
                            r_temp = sqrt(X_s(:,1).^2+X_s(:,2).^2+X_s(:,3).^2);
                            theta_temp = acos(X_s(:,3)./r_temp);
                            D = getStressTransformationMatrix(theta_temp,phi_temp,1);
                            for l = 1:6
                                D_kl = D(1, l, :);
                                D_kl = repmat(D_kl(:),1,length(omega));
                                sigma_rr = sigma_rr + D_kl.*sigma_cart{l};
                            end
                            layer{m}.err_pc = max(abs(p_tot+sigma_rr),[],1)./max(abs(p_tot),[],1);
                    end
                case {'solid','viscoelastic'}
                    sigma_cart = cell(6,1);
                    sigma_cart{1} = layer{m}.sigma_xx(indices,:);
                    sigma_cart{2} = layer{m}.sigma_yy(indices,:);
                    sigma_cart{3} = layer{m}.sigma_zz(indices,:);
                    sigma_cart{4} = layer{m}.sigma_yz(indices,:);
                    sigma_cart{5} = layer{m}.sigma_xz(indices,:);
                    sigma_cart{6} = layer{m}.sigma_xy(indices,:);

                    sigma_rr = zeros(size(sigma_cart{1}),prec);

                    phi_temp = atan2(X_s(:,2),X_s(:,1));
                    r_temp = sqrt(X_s(:,1).^2+X_s(:,2).^2+X_s(:,3).^2);
                    theta_temp = acos(X_s(:,3)./r_temp);
                    D = getStressTransformationMatrix(theta_temp,phi_temp,1);
                    for l = 1:6
                        D_kl = D(1, l, :);
                        D_kl = repmat(D_kl(:),1,length(omega));
                        sigma_rr = sigma_rr + D_kl.*sigma_cart{l};
                    end
                    switch nextMedia
                        case 'void'
                            if SSBC
                                layer{m}.err_pc = max(abs(sigma_rr),[],1)/P_inc;
                            end
                        case 'fluid'
                            p_tot = layer{m+1}.p(indices2,:);
                            if m+1 == 1
                                k = omega/layer{m+1}.c_f;
                                k_vec = options.d_vec*k;
                                p_inc = P_inc*exp(1i*X_s*k_vec);
                                p_tot = p_tot + p_inc;
                            end
                            layer{m}.err_pc = max(abs(p_tot+sigma_rr),[],1)./max(abs(p_tot),[],1);
                    end
            end
        end
        if layer{m}.calc_errDispCond && ~isSphere
            if m < M
                X2 = layer{m+1}.X;
                [indices, indices2] = findMatchingPoints(X,X2,Eps);
            else
                indices = abs(norm2(X) - R_i) < 10*Eps;
            end
            switch layer{m}.media
                case 'fluid'
                    P_inc = options.P_inc;

                    rho = layer{m}.rho;

                    n_x = X_s(:,1)./norm2(X_s);
                    n_y = X_s(:,2)./norm2(X_s);
                    n_z = X_s(:,3)./norm2(X_s);
                    n_x = repmat(n_x,1,length(omega));
                    n_y = repmat(n_y,1,length(omega));
                    n_z = repmat(n_z,1,length(omega));

                    dpdx = layer{m}.dpdx(indices,:);
                    dpdy = layer{m}.dpdy(indices,:);
                    dpdz = layer{m}.dpdz(indices,:);
                    if m == m_s
                        k = omega/layer{m}.c_f;
                        switch options.applyLoad
                            case 'planeWave'
                                k_vec = options.d_vec*k;
                                p_inc = P_inc*exp(1i*X_s*k_vec);
                                dpdx = dpdx + 1i*p_inc.*repmat(k_vec(1,:),size(p_inc,1),1);
                                dpdy = dpdy + 1i*p_inc.*repmat(k_vec(2,:),size(p_inc,1),1);
                                dpdz = dpdz + 1i*p_inc.*repmat(k_vec(3,:),size(p_inc,1),1);
                            otherwise
                                error('Not implemented')
                        end
                    end
                    switch nextMedia
                        case 'void'
                            if SHBC
                                layer{m}.err_dc = max(abs(dpdx.*n_x + dpdy.*n_y + dpdz.*n_z),[],1)/P_inc;
                            end
                        case {'solid','viscoelastic'}
                            u_x = layer{m+1}.u_x(indices2,:);
                            u_y = layer{m+1}.u_y(indices2,:);
                            u_z = layer{m+1}.u_z(indices2,:);
                            Omega = repmat(omega,size(u_x,1),1);
                            layer{m}.err_dc = max(abs(    (dpdx-rho*Omega.^2.*u_x).*n_x ...
                                                        + (dpdy-rho*Omega.^2.*u_y).*n_y ...
                                                        + (dpdz-rho*Omega.^2.*u_z).*n_z),[],1)./max(abs(dpdx.*n_x + dpdy.*n_y + dpdz.*n_z),[],1);
                    end
                case {'solid','viscoelastic'}
                    u_x = layer{m}.u_x(indices,:);
                    u_y = layer{m}.u_y(indices,:);
                    u_z = layer{m}.u_z(indices,:);
                    Omega = repmat(omega,size(u_x,1),1);
                    switch nextMedia
                        case 'fluid'
                            n_x = X_s(:,1)./norm2(X_s);
                            n_y = X_s(:,2)./norm2(X_s);
                            n_z = X_s(:,3)./norm2(X_s);
                            n_x = repmat(n_x,1,length(omega));
                            n_y = repmat(n_y,1,length(omega));
                            n_z = repmat(n_z,1,length(omega));
                            rho = layer{m+1}.rho;

                            dpdx = layer{m+1}.dpdx(indices2,:);
                            dpdy = layer{m+1}.dpdy(indices2,:);
                            dpdz = layer{m+1}.dpdz(indices2,:);
                            layer{m}.err_dc = max(abs(    (dpdx-rho*Omega.^2.*u_x).*n_x ...
                                                        + (dpdy-rho*Omega.^2.*u_y).*n_y ...
                                                        + (dpdz-rho*Omega.^2.*u_z).*n_z),[],1)./max(abs(dpdx.*n_x + dpdy.*n_y + dpdz.*n_z),[],1);
                    end
            end
        end
    end
    if M > 1 && ~(strcmp(nextMedia,'void') || strcmp(nextMedia,'origin'))
        if m+1 == M
            if layer{end}.R_i == 0
                nextMedia = 'origin';
            else
                nextMedia = 'void';
            end
        else
            nextMedia = layer{m+2}.media;
        end
    end
end
fieldsToRemove = {'r','theta','phi','P','dP','d2P','Z_zeta','Z_zeta_i','Z_zeta_o','Z_xi_i','Z_eta_i','Z_xi_o','Z_eta_o','a_temp','a','b','b_temp','k','k_temp','Z_r_s','K','G', ...
                  'calc_sigma_rr','calc_sigma_tt','calc_sigma_pp','calc_sigma_rt','calcStresses','calcCartesianDispDerivatives','calc_errors', ...
                  'calc_u_r','calc_u_t','calc_du_rdr','calc_du_rdt','calc_du_tdr','calc_du_tdt','calc_navier1','calc_navier2', ...
                  'calc_dpdr','calc_dpdt','calc_d2pdr2','calc_d2pdt2','calc_farFieldOnly'};
for m = 1:M
    for field = fieldsToRemove
        if isfield(layer{m},field{1})
            layer{m} = rmfield(layer{m},field{1});
        end
    end
end

function layer = getDefaultParameters(layer)

for m = 1:numel(layer)
    layer{m}.R_i              = 1;       % Inner radius of layer
    layer{m}.rho              = 1000;    % Mass density
    layer{m}.calc_errDispCond = false;   % Calculate the errors for the displacement conditions
    layer{m}.calc_errPresCond = false;   % Calculate the errors for the pressure conditions
    switch layer{m}.media
        case 'fluid'
            layer{m}.c_f          	= 1500;       % Speed of sound
            layer{m}.calc_p_0       = false;      % Toggle calculation of the far field pattern
            layer{m}.calc_p       	= false;      % Toggle calculation of the scattered pressure
            layer{m}.calc_dp      	= false(1,3); % Toggle calculation of the three components of the gradient of the pressure
            layer{m}.calc_p_laplace	= false;      % Toggle calculation of the Laplace operator of the scattered pressure fields
            layer{m}.calc_errHelm	= false;      % Toggle calculation of the errors for the Helmholtz equation
        case {'solid','viscoelastic'}
            layer{m}.E            = 200e9;      % Youngs modulus for solid layers
            layer{m}.nu           = 0.3;        % Poisson ratio for solid layers
            layer{m}.calc_u       = false(1,3); % Toggle calculation of the three components of the displacement
            layer{m}.calc_du      = false(3,3); % Toggle calculation of the three cartesian derivatives of the three components of the displacement [du_xdx du_xdy du_xdz; 
                                                %                                                                                                    du_ydx du_ydy du_ydz; 
                                                %                                                                                                    du_zdx du_zdy du_zdz]
            layer{m}.calc_sigma   = false(1,6); % Toggle calculation of the six components of the stress field (cartesian coordinates) [sigma_xx sigma_yy sigma_zz sigma_yz sigma_xz sigma_xy]
            layer{m}.calc_sigma_s = false(1,6); % Toggle calculation of the six components of the stress field (spherical coordinates) [sigma_rr sigma_tt sigma_pp sigma_tp sigma_rp sigma_rt]
            layer{m}.calc_errNav  = false;      % Toggle calculation of the errors for the Navier equation
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [layer,N_eps,flag] = e3Dss_0(layer, options)
% This function computes exact 3D scattering solutions when the axis of
% symmetry is the z-axis
%
%
% Note that dpdt, u_t, du_tdr, and du_rdt are scaled by csc(theta)
displayIter = strcmp(options.Display, 'iter');
prec = options.prec;
omega = options.omega;
nFreqs = length(omega);
M = numel(layer);
Eps = options.Eps;
N_max = options.N_max;
fluidFieldNames = {'p_0','p','dpdr','dpdt','d2pdr2','d2pdt2'};
solidFieldNames = {'u_r','u_t','du_rdr','du_rdt','du_tdr','du_tdt','sigma_rr','sigma_tt','sigma_pp','sigma_rt','navier1','navier2'};
if any(omega == 0)
    computeForStaticCase = true;
else
    computeForStaticCase = false;
end

R_i = inf;
for m = 1:M
    R_o = R_i;
    R_i = layer{m}.R_i;
    isSphere = R_i == 0;
    if isfield(layer{m} ,'X')
        n_X = size(layer{m}.X,1);
    else
        n_X = 0;
    end
    if n_X > 0
        layer{m}.P = zeros(2,n_X,prec); 
        layer{m}.dP = zeros(2,n_X,prec); 
        layer{m}.d2P = zeros(2,n_X,prec);
    end
    switch layer{m}.media
        case 'fluid'
            layer{m}.k = omega*(1./layer{m}.c_f);

            for i = 1:2
                for j = 1:2
                    if n_X > 0
                        layer{m}.Z_zeta{i,j} = zeros(nFreqs-computeForStaticCase,n_X,prec);
                    end
                    if ~isSphere
                        layer{m}.Z_zeta_i{i,j} = zeros(nFreqs,1,prec);
                    end
                    if ~isinf(R_o)
                        layer{m}.Z_zeta_o{i,j} = zeros(nFreqs,1,prec);
                    end
                    if strcmp(options.applyLoad,'pointCharge') && m == options.m_s
                        layer{m}.Z_r_s{i,j} = zeros(nFreqs,1,prec);
                    end
                end
            end
            if n_X > 0
                for fieldName = fluidFieldNames
                    if layer{m}.(['calc_' fieldName{1}])
                        layer{m}.(fieldName{1}) = zeros(nFreqs,n_X,prec);
                    end
                end
            end
        case {'solid','viscoelastic'}
            E = layer{m}.E;
            nu = layer{m}.nu;
            rho = layer{m}.rho;
            % Compute derived quantities
            K = E/(3*(1-2*nu));
            G = E/(2*(1+nu));

            c_s_1 = sqrt((3*K+4*G)/(3*rho)); % longitudinal wave velocity 
            c_s_2 = sqrt(G/rho); % shear wave velocity

            layer{m}.K = K;
            layer{m}.G = G;
            layer{m}.a = omega/c_s_1;
            layer{m}.b = omega/c_s_2;

            for i = 1:2
                for j = 1:2
                    if n_X > 0
                        layer{m}.Z_xi{i,j} = zeros(nFreqs-computeForStaticCase,n_X,prec);
                        layer{m}.Z_eta{i,j} = zeros(nFreqs-computeForStaticCase,n_X,prec);
                    end
                    if ~isSphere
                        layer{m}.Z_xi_i{i,j} = zeros(nFreqs,1,prec);
                        layer{m}.Z_eta_i{i,j} = zeros(nFreqs,1,prec);
                    end
                    if ~isinf(R_o)
                        layer{m}.Z_xi_o{i,j} = zeros(nFreqs,1,prec);
                        layer{m}.Z_eta_o{i,j} = zeros(nFreqs,1,prec);
                    end
                end
            end
            if n_X > 0
                for fieldName = solidFieldNames
                    if layer{m}.(['calc_' fieldName{1}])
                        layer{m}.(fieldName{1}) = zeros(nFreqs,n_X,prec);
                    end
                end
            end
    end
end

indices = (1:nFreqs).';
Zindices = indices;
if computeForStaticCase
    if isa(options.P_inc,'function_handle')
        P_inc = options.P_inc(0);
    else
        P_inc = options.P_inc;
    end
    staticIdx = find(omega == 0);
    for m = 1:M
        if isfield(layer{m} ,'X')
            n_X = size(layer{m}.X,1);
        else
            n_X = 0;
        end
        if n_X > 0
            switch layer{m}.media
                case 'fluid'
                    if m > 1
                        layer{m}.p(staticIdx,:) = P_inc*ones(1,n_X,prec);
                    end
                case {'solid','viscoelastic'}
                    A = -P_inc/(3*layer{m}.K);
                    if layer{m}.calc_u_r
                        layer{m}.u_r(staticIdx,:) = A*repmat(layer{m}.r,nFreqs,1);
                    end
                    if layer{m}.calc_sigma_rr
                        layer{m}.sigma_rr(staticIdx,:) = A*ones(1,n_X,prec);
                    end
                    if layer{m}.calc_sigma_tt
                        layer{m}.sigma_tt(staticIdx,:) = A*ones(1,n_X,prec);
                    end
                    if layer{m}.calc_sigma_pp
                        layer{m}.sigma_pp(staticIdx,:) = A*ones(1,n_X,prec);
                    end
            end
        end
    end
    indices(staticIdx) = [];
    Zindices = Zindices(1:end-1);
end

% In order to avoid premeture termination of the series summation, the
% vector wasNonZero{m} contains the information about the magnitude of the
% previous contribution to the sum. In particular, if the previous
% nExtraTerms terms have relative magnitude less than Eps compared to the
% total sum, the summation will terminate (for the particular omega)
nExtraTerms = 2;
hasCnvrgd = zeros(nFreqs,nExtraTerms); % matrix of element that "has converged"
if computeForStaticCase
    hasCnvrgd(staticIdx,:) = 1;
end
if strcmp(prec,'sym')
    tiny = vpa('1e-1000'); % To avoid dividing by zero.
else
    tiny = realmin(prec); % To avoid dividing by zero.
end
n = zeros(1,prec);
N_eps = NaN(nFreqs,1);
flag = zeros(size(hasCnvrgd,1),1); % Program terminated successfully unless error occurs (for each frequency)
singleModeSolution = (strcmp(options.applyLoad,'pointCharge') && options.r_s == 0) || strcmp(options.applyLoad,'radialPulsation') ...
                     || (strcmp(options.applyLoad,'surfExcitation') && options.theta_s(1) == 0 && abs(options.theta_s(2) - pi) < Eps);

while n <= N_max && ~(singleModeSolution && n > 0)
    try % and hope that no spherical Bessel functions are evaluated to be too large
        if displayIter
            tic
        end
        omega_temp = omega(indices);
        R_i = inf;
        for m = 1:M
            R_o = R_i;
            R_i = layer{m}.R_i;
            isSphere = R_i == 0;
            switch layer{m}.media
                case 'fluid'
                    evalLoad = strcmp(options.applyLoad,'pointCharge') && m == options.m_s;
                    layer{m}.k_temp = layer{m}.k(indices,:);
                    if ~isSphere
                        layer{m}.Z_zeta_i = iterateZ(n,layer{m}.k_temp*R_i,layer{m}.Z_zeta_i,Zindices,~isSphere);
                    end
                    if ~isinf(R_o)
                        layer{m}.Z_zeta_o = iterateZ(n,layer{m}.k_temp*R_o,layer{m}.Z_zeta_o,Zindices,~isSphere || evalLoad);
                    end
                    if evalLoad
                        layer{m}.Z_r_s = iterateZ(n,layer{m}.k_temp*options.r_s,layer{m}.Z_r_s,Zindices,~isSphere);
                    end
                case {'solid','viscoelastic'}
                    layer{m}.a_temp = layer{m}.a(indices,:);
                    layer{m}.b_temp = layer{m}.b(indices,:);
                    if ~isSphere
                        layer{m}.Z_xi_i  = iterateZ(n,layer{m}.a_temp*R_i,layer{m}.Z_xi_i, Zindices,~isSphere);
                        layer{m}.Z_eta_i = iterateZ(n,layer{m}.b_temp*R_i,layer{m}.Z_eta_i,Zindices,~isSphere);
                    end
                    if ~isinf(R_o)
                        layer{m}.Z_xi_o  = iterateZ(n,layer{m}.a_temp*R_o,layer{m}.Z_xi_o, Zindices,~isSphere);
                        layer{m}.Z_eta_o = iterateZ(n,layer{m}.b_temp*R_o,layer{m}.Z_eta_o,Zindices,~isSphere);
                    end
            end
        end
        C = getCoeffs(n, omega_temp, layer, options);

        hasCnvrgdTmp = zeros(length(indices),M); % temporary hasCnvrgd matrix
        for m = 1:M
            isSphere = layer{m}.R_i == 0;
            hasCnvrgdTmp2 = ones(size(indices)); % temporary hasCnvrgd vector    
            if isfield(layer{m},'r')
                [layer{m}.P, layer{m}.dP, layer{m}.d2P] = legendreDerivs(n, cos(layer{m}.theta), layer{m}.P, layer{m}.dP, layer{m}.d2P);
                switch layer{m}.media
                    case 'fluid'
                        zeta = layer{m}.k_temp*layer{m}.r;
                        if ~layer{m}.calc_farFieldOnly
                            layer{m}.Z_zeta = iterateZ(n,zeta,layer{m}.Z_zeta,Zindices,~isSphere);
                        end
                        media = p_(m,n,zeta,layer{m}.theta,C{m},layer{m}.k_temp,layer{m}.P(2,:), layer{m}.dP(2,:), layer{m}.d2P(2,:), ...
                                    layer{m}.Z_zeta, layer{m},isSphere);
                        fieldNames = fluidFieldNames;
                    case {'solid','viscoelastic'}
                        xi = layer{m}.a_temp*layer{m}.r;
                        eta = layer{m}.b_temp*layer{m}.r;
                        layer{m}.Z_xi = iterateZ(n,xi,layer{m}.Z_xi,Zindices,~isSphere);
                        layer{m}.Z_eta = iterateZ(n,eta,layer{m}.Z_eta,Zindices,~isSphere);
                        if isSphere
                            A = C{m}(:,1);
                            if n > 0
                                B = C{m}(:,2);
                            else
                                B = NaN(1,2);
                            end
                        else
                            A = C{m}(:,1:2);
                            if n > 0
                                B = C{m}(:,3:4);
                            else
                                B = NaN(1,2);
                            end
                        end
                        media = u_(n,layer{m}.r,layer{m}.theta,A,B,xi,eta,layer{m}.G,layer{m}.K,layer{m}.a_temp,layer{m}.b_temp,...
                                    layer{m}.P(2,:),layer{m}.dP(2,:),layer{m}.d2P(2,:), layer{m}.Z_xi, layer{m}.Z_eta, isSphere, layer{m}, omega_temp);
                        fieldNames = solidFieldNames;
                end
                for fieldName = fieldNames
                    if layer{m}.(['calc_' fieldName{1}])
                        [layer,hasCnvrgdTmp2] = updateSum(layer,m,media,fieldName{1},indices,hasCnvrgdTmp2,tiny,Eps);
                    end
                end
            end
            hasCnvrgdTmp(logical(hasCnvrgdTmp2),m) = 1;
        end
        hasCnvrgd(indices,:) = [hasCnvrgd(indices,2:end), prod(hasCnvrgdTmp,2)];
        indicesPrev = indices;
        indices = find(~prod(hasCnvrgd,2));
        [~,Zindices] = ismember(indices,indicesPrev);
        if length(indices) < length(indicesPrev)
            N_eps(setdiff(indicesPrev,indices)) = n;
        end
        if isempty(indices) % every element has converged
            break;
        end
        if displayIter
            fprintf('Completed calculation of term n = %d using %g seconds.\n', n, toc)
        end
        n = n + 1;
    catch ME
        flag = -~prod(hasCnvrgd,2);
        if strcmp(ME.identifier, 'e3Dss:infBessel')
            warning('e3Dss:infBessel','The summation ended prematurely at n = %d because a Bessel function evaluation was too large.', n)
        elseif strcmp(ME.identifier, 'e3Dss:singularK')
            warning('e3Dss:singularK','The summation ended prematurely at n = %d because the global matrix was singular to working precision.', n)
        else
            rethrow(ME)
        end
        break
    end
end
if strcmp(options.Display, 'final') || strcmp(options.Display, 'iter')
    if n-1 == N_max
        warning('e3Dss:N_max_reached','The summation did not converge using N_max = %d terms.', N_max)
    elseif ~any(flag)
        if singleModeSolution
            fprintf('Eps precision reached with a N_eps = 1 terms.\n')
        else
            fprintf('Eps precision reached with N_eps = %g terms.\n', max(N_eps))
        end
    end
end
% Scale functions with the amplitude of the incident wave
P_incIsVec = isa(options.P_inc,'function_handle') ;
if P_incIsVec
    P_inc = options.P_inc(omega);
else
    P_inc = options.P_inc;
end
for m = 1:M
    switch layer{m}.media
        case 'fluid'
            fieldNames = fluidFieldNames;
        case {'solid','viscoelastic'}
            fieldNames = solidFieldNames;
    end
    for fieldName = fieldNames
        if layer{m}.(['calc_' fieldName{1}])
            if P_incIsVec
                layer{m}.(fieldName{1}) = repmat(P_inc,1,size(layer{m}.(fieldName{1}),2)).*layer{m}.(fieldName{1});
            else
                layer{m}.(fieldName{1}) = repmat(P_inc,size(layer{m}.(fieldName{1}))).*layer{m}.(fieldName{1});
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Z = iterateZ(n,x,Z,Zindices,evalBessely)

if n == 0
    Z{1,2} = bessel_s(n,x,1);
    if evalBessely
        Z{2,2} = bessel_s(n,x,2);
    end
end
Z{1,1} = Z{1,2}(Zindices,:);
Z{1,2} = bessel_s(n+1,x,1);
if evalBessely
    Z{2,1} = Z{2,2}(Zindices,:);
    Z{2,2} = bessel_s(n+1,x,2);    
end

function [layer,hasCnvrgd] = updateSum(layer,m,media,fieldName,indices,hasCnvrgd,tiny,Eps)

layer{m}.(fieldName)(indices,:) = layer{m}.(fieldName)(indices,:) + media.(fieldName);
hasCnvrgd = hasCnvrgd.*prod(abs(media.(fieldName))./(abs(layer{m}.(fieldName)(indices,:))+tiny) < Eps, 2);

function fluid = p_(m,n,zeta,theta,C,k,P,dP,d2P,Z,layer,isSphere)
% Note that in the case of isSphere and zeta = 0: 
% --- dpdz =: dpdr and dpdx = dpdy = dpdt = 0
% --- nabla p =: d2pdr2, d2pdt2 := 0
% Also note that dpdt is scaled by csc(theta)

Q0 = P;
if layer.calc_p || layer.calc_dpdr || layer.calc_dpdt || layer.calc_d2pdr2 || layer.calc_d2pdr2
    j_n = Z{1,1};
    if layer.calc_dpdr
        dj_n = dbessel_s(n,zeta,1,Z);
    end
    if layer.calc_d2pdr2
        d2j_n = d2bessel_s(n,zeta,1,Z);
    end
    if m == 1
        y_n = Z{2,1};
        
        h_n = j_n + 1i*y_n;
        if layer.calc_dpdr
            dy_n = dbessel_s(n,zeta,2,Z);
            dh_n = dj_n + 1i*dy_n;
        end
        if layer.calc_d2pdr2
            d2y_n = d2bessel_s(n,zeta,2,Z);
            d2h_n = d2j_n + 1i*d2y_n;
        end
    elseif ~isSphere
        y_n = Z{2,1};
        if layer.calc_dpdr
            dy_n = dbessel_s(n,zeta,2,Z);
        end
        if layer.calc_d2pdr2
            d2y_n = d2bessel_s(n,zeta,2,Z);
        end
    end
end

if layer.calc_dpdt
    Q1 = Q_(1,theta,P,dP,d2P,true);
end
if layer.calc_d2pdt2
    Q2 = Q_(2,theta,P,dP,d2P);
end

if m == 1
    if layer.calc_p_0
        h_n_0 = repmat(1i^(-n-1)./k, 1, size(zeta,2));
        fluid.p_0 = C*Q0.*h_n_0;
    end
%     if layer.calc_dp_0dr
%         dh_n_0  = repmat(1i^(-n), size(zeta,1),size(zeta,2));
%     end
%     if layer.calc_d2p_0dr2
%         d2h_n_0 = repmat(1i^(-n+1)*k, 1, size(zeta,2));
%     end
    if layer.calc_p
        fluid.p = C*Q0.*h_n;
    end
    if layer.calc_dpdr
        fluid.dpdr = C.*k*Q0.*dh_n;
    end
    if layer.calc_dpdt % calculate a scaled dpdt (scaled by csc)
        fluid.dpdt = C*Q1.*h_n;
    end
    if layer.calc_d2pdr2
        fluid.d2pdr2 = C.*k.^2*Q0.*d2h_n;
    end
    if layer.calc_d2pdt2
        fluid.d2pdt2 = C*Q2.*h_n;
    end
elseif isSphere
    if layer.calc_p
        fluid.p = C*Q0.*j_n;
    end
    if layer.calc_dpdr
        fluid.dpdr = C.*k*Q0.*dj_n;
        indices = logical(zeta(1,:) < eps);
        if n == 1
            fluid.dpdr(:,indices) = repmat(k/3.*C,1,sum(indices));
        else
            fluid.dpdr(:,indices) = 0;
        end
    end
    if layer.calc_dpdt % calculate a scaled dpdt (scaled by csc)
        fluid.dpdt = C*Q1.*j_n;
        indices = logical(zeta(1,:) < eps);
        fluid.dpdt(:,indices) = 0;
    end
    if layer.calc_d2pdr2
        fluid.d2pdr2 = C.*k.^2*Q0.*d2j_n;
        indices = logical(zeta(1,:) < eps);
        if n == 0
            fluid.d2pdr2(:,indices) = repmat(-k.^2.*C,1,sum(indices));
        else
            fluid.d2pdr2(:,indices) = 0;
        end
    end
    if layer.calc_d2pdt2
        fluid.d2pdt2 = C*Q2.*j_n;
        indices = logical(zeta(1,:) < eps);
        fluid.d2pdt2(:,indices) = 0;
    end
else
    if layer.calc_p
        fluid.p = C(:,1)*Q0.*j_n + C(:,2)*Q0.*y_n;
    end
    if layer.calc_dpdr
        fluid.dpdr = C(:,1).*k*Q0.*dj_n + C(:,2).*k*Q0.*dy_n;
    end
    if layer.calc_dpdt % calculate a scaled dpdt (scaled by csc)
        fluid.dpdt = C(:,1)*Q1.*j_n + C(:,2)*Q1.*y_n;
    end
    if layer.calc_d2pdr2
        fluid.d2pdr2 = C(:,1).*k.^2*Q0.*d2j_n + C(:,2).*k.^2*Q0.*d2y_n;
    end
    if layer.calc_d2pdt2
        fluid.d2pdt2 = C(:,1)*Q2.*j_n + C(:,2)*Q2.*y_n;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function solid = u_(n,r,theta,A,B,xi,eta,G,K,a,b,P,dP,d2P,Zxi,Zeta,isSphere,layer,omega)
% Note that in the case of isSpherhe and r = 0: 
% --- u_z =: u_r and u_x = u_y = u_t = 0
% --- du_xdx =: du_rdr, du_ydy =: du_rdt, du_zdz =: du_tdr, 0 =: du_tdt
% --- sigma_11 =: sigma_rr, sigma_22 =: sigma_tt, sigma_33 =: sigma_pp, 0 =: sigma_rt
% Also note that u_t, du_tdr and du_rdt are scaled by csc(theta)

Q0 = Q_(0,theta,P,dP,d2P);
Q1s = Q_(1,theta,P,dP,d2P,true);
Q1 = sin(theta).*Q1s;
Q2 = Q_(2,theta,P,dP,d2P);

r2 = r.^2;

if layer.calc_u_r
    Q0r = Q0./r;
    u_r = A(:,1)*Q0r.*S_(1,1,n,xi,eta,Zxi);
    if n > 0
    	u_r = u_r + B(:,1)*Q0r.*T_(1,1,n,eta,Zeta);
    end
    if isSphere
        indices = logical(r < eps);
        if n == 1
            u_r(:,indices) = repmat((a.*A(:,1) - 2*b.*B(:,1))/3,1,sum(indices));
        else
            u_r(:,indices) = 0;
        end
    else
        u_r = u_r + A(:,2)*Q0r.*S_(1,2,n,xi,eta,Zxi);
        if n > 0
        	u_r = u_r + B(:,2)*Q0r.*T_(1,2,n,eta,Zeta);
        end
    end
    solid.u_r = u_r;
end

if layer.calc_u_t
    Q1sr = Q1s./r;
    u_t = A(:,1)*Q1sr.*S_(2,1,n,xi,eta,Zxi);
    if n > 0
    	u_t = u_t + B(:,1)*Q1sr.*T_(2,1,n,eta,Zeta);
    end
    if isSphere
        u_t(:,logical(r < eps)) = 0;
    else
        u_t = u_t + A(:,2)*Q1sr.*S_(2,2,n,xi,eta,Zxi);
        if n > 0
            u_t = u_t + B(:,2)*Q1sr.*T_(2,2,n,eta,Zeta);
        end
    end
    solid.u_t = u_t;
end

if layer.calc_du_rdr
    Q0r2 = Q0./r2;
    du_rdr = A(:,1)*Q0r2.*S_(3,1,n,xi,eta,Zxi);
    if n > 0
        du_rdr = du_rdr + B(:,1)*Q0r2.*T_(3,1,n,eta,Zeta);
    end
    if isSphere
        indices = logical(r < eps);
        if n == 0
            du_rdr(:,indices) = repmat(G/K*(4*a.^2-3*b.^2).*A(:,1)/9,1,sum(indices));
        elseif n == 2
            du_rdr(:,indices) = repmat(-(a.^2.*A(:,1)-3*b.^2.*B(:,1))/15,1,sum(indices));
        else
            du_rdr(:,indices) = 0;
        end
    else
        du_rdr = du_rdr + A(:,2)*Q0r2.*S_(3,2,n,xi,eta,Zxi);
        if n > 0
            du_rdr = du_rdr + B(:,2)*Q0r2.*T_(3,2,n,eta,Zeta);
        end
    end
    solid.du_rdr = du_rdr;
end

if layer.calc_du_rdt
    Q1sr = Q1s./r;
    du_rdt = A(:,1)*Q1sr.*S_(1,1,n,xi,eta,Zxi);
    if n > 0
        du_rdt = du_rdt + B(:,1)*Q1sr.*T_(1,1,n,eta,Zeta);
    end
    if isSphere
        indices = logical(r < eps);
        if n == 0
            du_rdt(:,indices) = repmat(G/K*(4*a.^2-3*b.^2).*A(:,1)/9,1,sum(indices));
        elseif n == 2
            du_rdt(:,indices) = repmat(-(a.^2.*A(:,1)-3*b.^2.*B(:,1))/15,1,sum(indices));
        else
            du_rdt(:,indices) = 0;
        end
    else
        du_rdt = du_rdt + A(:,2)*Q1sr.*S_(1,2,n,xi,eta,Zxi);
        if n > 0
        	du_rdt = du_rdt + B(:,2)*Q1sr.*T_(1,2,n,eta,Zeta);
        end
    end
    solid.du_rdt = du_rdt;
end

if layer.calc_du_tdr
    Q1sr2 = Q1s./r2;
    du_tdr = A(:,1)*Q1sr2.*S_(4,1,n,xi,eta,Zxi);
    if n > 0
        du_tdr = du_tdr + B(:,1)*Q1sr2.*T_(4,1,n,eta,Zeta);
    end
    if isSphere
        indices = logical(r < eps);
        if n == 0
            du_tdr(:,indices) = repmat(G/K*(4*a.^2-3*b.^2).*A(:,1)/9,1,sum(indices));
        elseif n == 2
            du_tdr(:,indices) = repmat(2*(a.^2.*A(:,1)-3*b.^2.*B(:,1))/15,1,sum(indices));
        else
            du_tdr(:,indices) = 0;
        end
    else
        du_tdr = du_tdr + A(:,2)*Q1sr2.*S_(4,2,n,xi,eta,Zxi);
        if n > 0
        	du_tdr = du_tdr + B(:,2)*Q1sr2.*T_(4,2,n,eta,Zeta);
        end
    end
    solid.du_tdr = du_tdr;
end

if layer.calc_du_tdt
    Q2r1 = Q2./r;
    du_tdt = A(:,1)*Q2r1.*S_(2,1,n,xi,eta,Zxi);
    if n > 0
        du_tdt = du_tdt + B(:,1)*Q2r1.*T_(2,1,n,eta,Zeta);
    end
    if isSphere
        du_tdt(:,logical(r < eps)) = 0;
    else
        du_tdt = du_tdt + A(:,2)*Q2r1.*S_(2,2,n,xi,eta,Zxi);
        if n > 0
        	du_tdt = du_tdt + B(:,2)*Q2r1.*T_(2,2,n,eta,Zeta);
        end
    end
    solid.du_tdt = du_tdt;
end

if layer.calc_sigma_rr
    Q0r2 = Q0./r2;
    sigma_rr = A(:,1)*Q0r2.*S_(5,1,n,xi,eta,Zxi);
    if n > 0
        sigma_rr = sigma_rr + B(:,1)*Q0r2.*T_(5,1,n,eta,Zeta);
    end
    if isSphere
        indices = logical(r < eps);
        if n == 0
            sigma_rr(:,indices) = repmat(5*(4*a.^2-3*b.^2).*A(:,1)/30,1,sum(indices));
        elseif n == 2
            sigma_rr(:,indices) = repmat((-2*a.^2.*A(:,1) + 6*b.^2.*B(:,1))/30,1,sum(indices));
        else
            sigma_rr(:,indices) = 0;
        end
    else
        sigma_rr = sigma_rr + A(:,2)*Q0r2.*S_(5,2,n,xi,eta,Zxi);
        if n > 0
        	sigma_rr = sigma_rr + B(:,2)*Q0r2.*T_(5,2,n,eta,Zeta);
        end
    end
	sigma_rr = 2*G*sigma_rr;
    
	solid.sigma_rr = sigma_rr;
end
if layer.calc_sigma_tt
    Q0r2 = Q0./r2;
    Q2r2 = Q2./r2;
    sigma_tt =   A(:,1)*Q0r2.*S_(6,1,n,xi,eta,Zxi)  ...
               + A(:,1)*Q2r2.*S_(2,1,n,xi,eta,Zxi);
    if n > 0
        sigma_tt = sigma_tt + B(:,1)*Q0r2.*T_(6,1,n,eta,Zeta) ...
                            + B(:,1)*Q2r2.*T_(2,1,n,eta,Zeta);
    end
    if isSphere
        indices = logical(r < eps);
        if n == 0
            sigma_tt(:,indices) = repmat(5*(4*a.^2-3*b.^2).*A(:,1)/30,1,sum(indices));
        elseif n == 2
            sigma_tt(:,indices) = repmat((-2*a.^2.*A(:,1) + 6*b.^2.*B(:,1))/30,1,sum(indices));
        else
            sigma_tt(:,indices) = 0;
        end
    else
        sigma_tt = sigma_tt + A(:,2)*Q0r2.*S_(6,2,n,xi,eta,Zxi) ...
                            + A(:,2)*Q2r2.*S_(2,2,n,xi,eta,Zxi);
        if n > 0
            sigma_tt = sigma_tt + B(:,2)*Q0r2.*T_(6,2,n,eta,Zeta) ...
                                + B(:,2)*Q2r2.*T_(2,2,n,eta,Zeta);
        end
    end
	sigma_tt = 2*G*sigma_tt;
    
	solid.sigma_tt = sigma_tt;
end
if layer.calc_sigma_pp
    Q0r2 = Q0./r2;
    cott_Q1r2 = cos(theta).*Q1s./r2;
    sigma_pp =   A(:,1)*Q0r2.*S_(6,1,n,xi,eta,Zxi) ...
               + A(:,1)*cott_Q1r2.*S_(2,1,n,xi,eta,Zxi);
    if n > 0
        sigma_pp = sigma_pp + B(:,1)*Q0r2.*T_(6,1,n,eta,Zeta) ...
                            + B(:,1)*cott_Q1r2.*T_(2,1,n,eta,Zeta);
    end
    if isSphere
        indices = logical(r < eps);
        if n == 0
            sigma_pp(:,indices) = repmat(5*(4*a.^2-3*b.^2).*A(:,1)/30,1,sum(indices));
        elseif n == 2
            sigma_pp(:,indices) = repmat((4*a.^2.*A(:,1) - 12*b.^2.*B(:,1))/30,1,sum(indices));
        else
            sigma_pp(:,indices) = 0;
        end
    else
        sigma_pp = sigma_pp + A(:,2)*Q0r2.*S_(6,2,n,xi,eta,Zxi) ...
                            + A(:,2)*cott_Q1r2.*S_(2,2,n,xi,eta,Zxi);
        if n > 0
            sigma_pp = sigma_pp + B(:,2)*Q0r2.*T_(6,2,n,eta,Zeta) ...
                                + B(:,2)*cott_Q1r2.*T_(2,2,n,eta,Zeta);
        end
    end
	sigma_pp = 2*G*sigma_pp;
    
	solid.sigma_pp = sigma_pp;
end
if layer.calc_sigma_rt
    Q1r2 = Q1./r2;
    sigma_rt = A(:,1)*Q1r2.*S_(7,1,n,xi,eta,Zxi);
    if n > 0
        sigma_rt = sigma_rt + B(:,1)*Q1r2.*T_(7,1,n,eta,Zeta);
    end
    if isSphere
        sigma_rt(:,logical(r < eps)) = 0;
    else
        sigma_rt = sigma_rt + A(:,2)*Q1r2.*S_(7,2,n,xi,eta,Zxi);
        if n > 0
            sigma_rt = sigma_rt + B(:,2)*Q1r2.*T_(7,2,n,eta,Zeta);
        end
    end
	sigma_rt = 2*G*sigma_rt;
    
	solid.sigma_rt = sigma_rt;
end
if layer.calc_errNav
    Q1sr2 = Q1s./r2;
    sigma_rt = A(:,1)*Q1sr2.*S_(7,1,n,xi,eta,Zxi);
    if n > 0
        sigma_rt = sigma_rt + B(:,1)*Q1sr2.*T_(7,1,n,eta,Zeta);
    end
    if ~isSphere
        sigma_rt = sigma_rt + A(:,2)*Q1sr2.*S_(7,2,n,xi,eta,Zxi);
        if n > 0
        	sigma_rt = sigma_rt + B(:,2)*Q1sr2.*T_(7,2,n,eta,Zeta);
        end
    end
	sigma_rt = 2*G*sigma_rt;
    
    r3 = r.^3;
    Q0r3 = Q0./r3;
    Q1r3 = Q1./r3;
    Q2r3 = Q2./r3;
    dsigma_rr_dr = A(:,1)*Q0r3.*S_(8,1,n,xi,eta,Zxi);
    if n > 0
    	dsigma_rr_dr = dsigma_rr_dr + B(:,1)*Q0r3.*T_(8,1,n,eta,Zeta);
    end
    if ~isSphere
        dsigma_rr_dr =  dsigma_rr_dr + A(:,2)*Q0r3.*S_(8,2,n,xi,eta,Zxi);
        if n > 0
        	dsigma_rr_dr =  dsigma_rr_dr + B(:,2)*Q0r3.*T_(8,2,n,eta,Zeta);
        end
    end
	dsigma_rr_dr = 2*G*dsigma_rr_dr;
      
    dsigma_rt_dr = A(:,1)*Q1r3.*S_(9,1,n,xi,eta,Zxi);
    if n > 0
        dsigma_rt_dr = dsigma_rt_dr + B(:,1)*Q1r3.*T_(9,1,n,eta,Zeta);
    end
    if ~isSphere
        dsigma_rt_dr =  dsigma_rt_dr + A(:,2)*Q1r3.*S_(9,2,n,xi,eta,Zxi);
        if n > 0
        	dsigma_rt_dr =  dsigma_rt_dr + B(:,2)*Q1r3.*T_(9,2,n,eta,Zeta);
        end
    end
	dsigma_rt_dr = 2*G*dsigma_rt_dr;

    dsigma_rt_dt = A(:,1)*Q2r3.*S_(7,1,n,xi,eta,Zxi);
    if n > 0
        dsigma_rt_dt = dsigma_rt_dt + B(:,1)*Q2r3.*T_(7,1,n,eta,Zeta);
    end
    if ~isSphere
        dsigma_rt_dt = dsigma_rt_dt + A(:,2)*Q2r3.*S_(7,2,n,xi,eta,Zxi);
        if n > 0
        	dsigma_rt_dt = dsigma_rt_dt + B(:,2)*Q2r3.*T_(7,2,n,eta,Zeta);
        end
    end
	dsigma_rt_dt = 2*G*dsigma_rt_dt;

    dsigma_diffr =   A(:,1)*Q1r3.*S_(6,1,n,xi,eta,Zxi) ...
                   + (-n^2-n+1)*A(:,1)*Q1r3.*S_(2,1,n,xi,eta,Zxi);
    if n > 0
        dsigma_diffr = dsigma_diffr +  B(:,1)*Q1r3.*T_(6,1,n,eta,Zeta) ...
                                    + (-n^2-n+1)*B(:,1)*Q1r3.*T_(2,1,n,eta,Zeta);
    end
    if ~isSphere
        dsigma_diffr = dsigma_diffr + A(:,2)*Q1r3.*S_(6,2,n,xi,eta,Zxi) ...
                                    + (-n^2-n+1)*A(:,2)*Q1r3.*S_(2,2,n,xi,eta,Zxi);
        if n > 0
            dsigma_diffr = dsigma_diffr + B(:,2)*Q1r3.*T_(6,2,n,eta,Zeta) ...
                                        + (-n^2-n+1)*B(:,2)*Q1r3.*T_(2,2,n,eta,Zeta);
        end
    end
	dsigma_diffr = 2*G*dsigma_diffr;
    
    R = repmat(r,size(xi,1),1);
    Theta = repmat(theta,size(xi,1),1);
    solid.navier1 = dsigma_rr_dr + dsigma_rt_dt + 1./R.*(2*sigma_rr-sigma_tt-sigma_pp+sigma_rt.*cos(Theta));
    solid.navier2 = dsigma_rt_dr + dsigma_diffr + 3./R.*sigma_rt.*sin(Theta);
    if isSphere
        indices = logical(r < eps);
        if n == 1
            rho = layer.rho;
            solid.navier1(:,indices) = repmat(-omega.^2*rho.*(a.*A(:,1) - 2*b.*B(:,1))/3,1,sum(indices));
        else
            solid.navier1(:,indices) = 0;
        end
        solid.navier2(:,indices) = 0;
    end
end
