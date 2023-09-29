function [layer,N_eps,flag,relTermMaxArr] = e3Dss(newLayers, newOptions)

% The function e3Dss (exact 3D scattering solutions) computes the solution
% to scattering problems on multilayered (elastic of fluid) spherical
% shells.
%
% The function handles the case in which both X and omega are vectors. 
%
% See README.md for details.
% 
% Author: Jon Vegard Venås
% E-mail: JonVegard.Venas@sintef.no
% Institute: SINTEF Digital
% Release: 2
% Release date: 21/03/2020

options = struct('d_vec',           [0;0;1],  ... 	 % Direction of the incident wave
                'omega',            2*pi*1e3, ...    % Angular frequency (can be an array of angular frequencies)
                'P_inc',            1,    ...     	 % Amplitude of incident wave
                'N_max',            Inf,  ...        % Upper limit for the number of terms in the series
                'Display',          'iter', ...     % Print options ('final', 'iter' or 'none')
                'BC',               'SHBC', ...      % Boundary condition at the innermost layer: 'IBC', 'SHBC', 'SSBC' or 'NNBC'
                'applyLoad',        'planeWave', ... % Incident wave type: I.e. 'planeWave', 'pointCharge', 'mechExcitation', 'surfExcitation', 'radialPulsation'
                'r_s',              NaN, ...         % Distance to load source (for loads 'pointCharge', 'mechExcitation' and 'surfExcitation')
                'p_inc_fromSeries', false, ...       % Calculate p_inc by series expansions (in terms of Bessel functions)
                'nu_a',             100, ...         % Value of nu at which scaled asymptotic expansions are used in bessel_c (set nu_a = -1 to turn off scaling)
                'z',                1, ...           % Impedance parameter for the case BC = 'IBC'
                'Eps',              eps,  ...        % Small parameter for series truncation
                'debug',            false, ...       % Print additional warning messages
                'saveRelTermMax',   false, ...       % Save the maximum of the relative terms added to the series
                'prec',             'double');       % Precision of the calculations
if nargin > 1
    options = updateOptions(options,newOptions);
end

fieldsToRemove = {'r','theta','phi','a_temp','a','b','b_temp','k','k_temp','G_temp','K_temp','Z_r_s','K','G','P','dP','d2P', ...
                  'Z_xi','Z_eta','Z_zeta','Z_zeta_i','Z_zeta_o','Z_xi_i','Z_eta_i','Z_xi_o','Z_eta_o','Z_xi_i_s','Z_xi_o_s','Z_eta_i_s','Z_eta_o_s','Z_zeta_i_s','Z_zeta_o_s','Z_r_s', ...
                  'calc_sigma_rr','calc_sigma_tt','calc_sigma_pp','calc_sigma_rt','calcStresses', ...
                  'dpdr','dpdt','d2pdr2','d2pdt2','dp_incdr','dp_incdt','u_r','u_t','sigma_rr','sigma_tt','sigma_pp','sigma_rt', ...
                  'calcCartesianDispDerivatives','calc_errors', ...
                  'calc_u_r','calc_u_t','calc_du_rdr','calc_du_rdt','calc_du_tdr','calc_du_tdt','calc_navier1','calc_navier2', ...
                  'calc_dpdr','calc_dpdt','calc_d2pdr2','calc_d2pdt2','calc_dp_incdr','calc_dp_incdt','calc_farFieldOnly'};

M = numel(newLayers);
fieldExist = cell(1,M);
for m = 1:M
    fieldExist{m} = false(size(fieldsToRemove));
    for i = 1:numel(fieldsToRemove)
        if isfield(newLayers{m},fieldsToRemove{i})
            fieldExist{m}(i) = true;
        end
    end
end
layer = getDefaultParameters(newLayers);

layer = updateOptions(layer,newLayers);


% Supress warning for ill-conditioned matrices in getCoeffs.m
if options.debug
    debugStr = 'on';
else
    debugStr = 'off';
end
poolobj = gcp('nocreate');
singMatrixWarning = 'off';
if isempty(poolobj)
    warning(singMatrixWarning, 'MATLAB:nearlySingularMatrix')
    warning(singMatrixWarning, 'MATLAB:singularMatrix')
    warning(singMatrixWarning, 'MATLAB:illConditionedMatrix')
    
    warning(debugStr, 'e3Dss:infBessel')
    warning(debugStr, 'e3Dss:divergeBessel')
    warning(debugStr, 'e3Dss:infWeight')
    warning(debugStr, 'e3Dss:singularK')
else
    pctRunOnAll(['warning(''' singMatrixWarning ''', ''MATLAB:singularMatrix'')'])
    pctRunOnAll(['warning(''' singMatrixWarning ''', ''MATLAB:nearlySingularMatrix'')'])
    pctRunOnAll(['warning(''' singMatrixWarning ''', ''MATLAB:illConditionedMatrix'')'])
    pctRunOnAll(['warning(''' debugStr ''', ''e3Dss:infBessel'')'])
    pctRunOnAll(['warning(''' debugStr ''', ''e3Dss:divergeBessel'')'])
    pctRunOnAll(['warning(''' debugStr ''', ''e3Dss:infWeight'')'])
    pctRunOnAll(['warning(''' debugStr ''', ''e3Dss:singularK'')'])
end
prec = options.prec;
PI = getC(prec,'pi');

omega = options.omega;
if isrow(omega)
    omega = omega.';
    options.omega = omega;
end
if isrow(options.z)
    options.z = options.z.';
end
if strcmp(layer{1}.media,'solid') % the outermost unbounded domain is solid
    error('This case is not implemented')
end
M = numel(layer);
noEvaluationPts = true;

switch options.applyLoad
    case {'planeWave','radialPulsation'}
        m_s = 1;
    case {'pointCharge','mechExcitation','surfExcitation','custom'}
        if ~isfield(options,'r_s')
            if strcmp(options.applyLoad, 'pointCharge')
                options.r_s = 2*layer{1}.R;
            elseif strcmp(options.applyLoad, 'surfExcitation') || strcmp(options.applyLoad, 'custom')
                options.r_s = layer{1}.R;
            end
        end
        if strcmp(options.applyLoad, 'surfExcitation') && ~isfield(options,'theta_s')
            options.theta_s = [0,PI/2];
        end
        if strcmp(options.applyLoad, 'custom') && ~isfield(options,'theta_s')
            options.theta_s = 1;
        end
        r_s = options.r_s;
        m_s = 1;
        for m = 1:M
            if r_s < layer{m}.R
                m_s = m_s + 1;
            else
                break
            end
        end
end
options.m_s = m_s;

% Update "calc_..." fields and calculate points in spherical coordinate
% system
for m = 1:M
    calc_err_pc = layer{m}.calc_err_pc;
    calc_err_dc = layer{m}.calc_err_dc;
    if size(layer{m}.lossFactor,2) == 1
        layer{m}.lossFactor(:,2) = layer{m}.lossFactor(:,1);
    end
    switch layer{m}.media
        case 'fluid'            
            calc_err_helmholtz = layer{m}.calc_err_helmholtz;
            calc_dp = or(layer{m}.calc_dp, calc_err_dc);

            calc_p_laplace = layer{m}.calc_p_laplace;

            layer{m}.calc_p = layer{m}.calc_p || calc_err_pc;
            layer{m}.calc_dpdr = any(calc_dp) || calc_p_laplace;
            layer{m}.calc_dp = calc_dp;
            layer{m}.calc_p_laplace = calc_p_laplace || calc_err_helmholtz;
            layer{m}.calc_dpdt = layer{m}.calc_dpdr;
            layer{m}.calc_d2pdr2 = calc_p_laplace || calc_err_helmholtz;
            layer{m}.calc_d2pdt2 = calc_p_laplace || calc_err_helmholtz;
            layer{m}.calc_errors = calc_err_dc || calc_err_pc || calc_err_helmholtz;
            layer{m}.calc_farFieldOnly = ~any([layer{m}.calc_p, layer{m}.calc_dpdr, layer{m}.calc_dpdt, layer{m}.calc_d2pdr2, layer{m}.calc_d2pdt2]);
            
            layer{m}.calc_p_inc = (layer{m}.calc_p_inc || layer{m}.calc_errors);
            layer{m}.calc_dp_inc = or(layer{m}.calc_dp_inc, layer{m}.calc_errors);
            layer{m}.calc_dp_incdr = and(any(layer{m}.calc_dp_inc), options.p_inc_fromSeries);
            layer{m}.calc_dp_incdt = and(any(layer{m}.calc_dp_inc), options.p_inc_fromSeries);
        case {'solid','viscoelastic'}
            calc_u = layer{m}.calc_u;
            calcCartesianDispDerivatives = any(layer{m}.calc_du(:));
            calc_errNav = any(layer{m}.calc_err_navier);
            layer{m}.calc_u = or(calc_u, calc_err_dc);
            if ~(isfield(layer{m},'calc_u_r') && layer{m}.calc_u_r)
                layer{m}.calc_u_r = calc_err_dc || calcCartesianDispDerivatives || any(calc_u) || calc_errNav;
            end
            if ~(isfield(layer{m},'calc_u_t') && layer{m}.calc_u_t)
                layer{m}.calc_u_t = calc_err_dc || calcCartesianDispDerivatives || any(calc_u) || calc_errNav;
            end
            layer{m}.calc_du_rdr = calcCartesianDispDerivatives;
            layer{m}.calc_du_rdt = calcCartesianDispDerivatives;
            layer{m}.calc_du_tdr = calcCartesianDispDerivatives;
            layer{m}.calc_du_tdt = calcCartesianDispDerivatives;
            calc_sigma = layer{m}.calc_sigma;
            calc_sigma_s = layer{m}.calc_sigma_s;
            if any(calc_sigma_s(2:end))
                error('Not yet implemented')
            end
            calcStresses = any(calc_sigma) || any(calc_sigma_s) || calc_errNav || calc_err_pc;
            for stressField = {'calc_sigma_rr','calc_sigma_tt','calc_sigma_pp','calc_sigma_rt'}
                if ~(isfield(layer{m},stressField{1}) && layer{m}.(stressField{1}))
                    layer{m}.(stressField{1}) = calcStresses;
                end
            end
            layer{m}.calcStresses = calcStresses;
            layer{m}.calc_navier1 = calc_errNav;
            layer{m}.calc_navier2 = calc_errNav;
            layer{m}.calc_sigma = or(calc_sigma, calcStresses);
            layer{m}.calcCartesianDispDerivatives = calcCartesianDispDerivatives;
            layer{m}.calc_errors = calc_err_dc || calc_err_pc || calc_errNav;
    end
    if isfield(layer{m},'X')
        %% Coordinate transformation
        % Due to symmetry, we can do a coordinate transformation such that we only
        % need to compute the solution for the special case k_vec = k*[0, 0, 1].
        noEvaluationPts = false;
        [r, theta, phi, A] = coordTransform(layer{m}.X, options.d_vec);
        layer{m}.r = r.';
        layer{m}.theta = theta.';
        layer{m}.phi = phi.';
    end
end
if isempty(omega)
    error('omega cannot be an empty array')
end
if noEvaluationPts
    error('No evaluation point (layer{m}.X) was provided')
end
options.SHBC = strcmp(options.BC,'SHBC');
options.SSBC = strcmp(options.BC,'SSBC');
options.IBC = strcmp(options.BC,'IBC');
if options.IBC && ~isfield(options,'z')
    options.z = 1;
end

if any(omega == 0) && (strcmp(options.applyLoad, 'mechExcitation') || ...
                      (strcmp(options.applyLoad, 'pointCharge') && options.r_s ~= 0) || ...
                      (strcmp(options.applyLoad, 'surfExcitation') && all(abs(options.theta_s - [0,PI]) < options.Eps)))
    warning('This case has no unique solution for omega=0, these values will be removed from the computation.')
    omega(omega == 0) = [];
    options.omega = omega;
end
if any([options.SHBC,options.SSBC,options.IBC]) && layer{end}.R == 0
    error('Boundary conditions may not be used at the origin')
end

%% Compute the solution with d_vec = [0, 0, 1]
[layer,N_eps,flag,relTermMaxArr] = e3Dss_0(layer, options);

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
                    dpdX{1} = zeros(nFreqs,n_X,prec);
                    dpdX{2} = zeros(nFreqs,n_X,prec);
                    dpdX{3} = zeros(nFreqs,n_X,prec);
                    for ii = 1:3
                        for jj = 1:3
                            dpdX{ii} = dpdX{ii} + A(ii,jj)*dpdX_m{jj};
                        end
                    end
                    for ii = 1:3
                        if layer{m}.calc_dp(ii)
                            layer{m}.dp{ii} = dpdX{ii};
                        end
                    end
                    
                    dpdt = layer{m}.dpdt; % use scaled version of dpdt for the laplace operator
                    if layer{m}.calc_p_laplace
                        d2pdr2 = layer{m}.d2pdr2;
                        d2pdt2 = layer{m}.d2pdt2;

                        temp = 2./R.*dpdr + 1./R.^2.*cos(Theta).*dpdt + 1./R.^2.*d2pdt2;                    
                        temp(indices) = 0;
                        layer{m}.p_laplace = d2pdr2 + temp;
                    end
                end
                P_incIsVec = isa(options.P_inc,'function_handle') ;
                if P_incIsVec
                    P_inc = options.P_inc(omega);
                else
                    P_inc = options.P_inc;
                end
                if m == m_s
                    switch options.applyLoad
                        case {'planeWave','pointCharge'}
                            if any(layer{m}.calc_dp_inc) || layer{m}.calc_p_inc
                                if options.p_inc_fromSeries
                                    if any(layer{m}.calc_dp_inc)
                                        dp_incdr = layer{m}.dp_incdr;
                                        dp_incdt = sin(Theta).*layer{m}.dp_incdt; % rescale dp_incdt
                                        dp_incdX_m = cell(3,1);

                                        dp_incdX_m{1} = dp_incdr.*sin(Theta).*cos(Phi) + dp_incdt.*cos(Theta).*cos(Phi)./R;
                                        dp_incdX_m{1}(indices) = 0;

                                        dp_incdX_m{2} = dp_incdr.*sin(Theta).*sin(Phi) + dp_incdt.*cos(Theta).*sin(Phi)./R;
                                        dp_incdX_m{2}(indices) = 0;

                                        dp_incdX_m{3} = dp_incdr.*cos(Theta) - dp_incdt.*sin(Theta)./R;
                                        dp_incdX_m{3}(indices) = dp_incdr(indices);
                                        layer{m}.dp_inc{ii} = cell(3,1);
                                        layer{m}.dp_inc{1} = zeros(nFreqs,n_X,prec);
                                        layer{m}.dp_inc{2} = zeros(nFreqs,n_X,prec);
                                        layer{m}.dp_inc{3} = zeros(nFreqs,n_X,prec);
                                        for ii = 1:3
                                            if layer{m}.calc_dp_inc(ii)
                                                for jj = 1:3
                                                    layer{m}.dp_inc{ii} = layer{m}.dp_inc{ii} + A(ii,jj)*dp_incdX_m{jj};
                                                end
                                            end
                                        end
                                    end
                                else
                                    k = omega/layer{m}.c_f;
                                    if strcmp(options.applyLoad,'planeWave')
                                        k_vec = k*options.d_vec.';
                                        if layer{m}.calc_p_inc
                                            layer{m}.p_inc = P_inc.*exp(1i*k_vec*layer{m}.X.');
                                        end
                                        for ii = 1:3
                                            if layer{m}.calc_dp_inc(ii)
                                                layer{m}.dp_inc{ii} = 1i*k_vec(:,ii).*layer{m}.p_inc;
                                            end
                                        end
                                    else
                                        x_s = r_s*options.d_vec.';
                                        Xxms = layer{m}.X-x_s;
                                        nXxms = norm2(Xxms).';
                                        p_inc = P_inc.*exp(1i*k.*nXxms)./nXxms;
                                        if layer{m}.calc_p_inc
                                            layer{m}.p_inc = p_inc;
                                        end
                                        for ii = 1:3
                                            if layer{m}.calc_dp_inc(ii)
                                                layer{m}.dp_inc{ii} = p_inc.*(1i*k - 1./nXxms).*Xxms(:,ii).'./nXxms;
                                            end
                                        end
                                    end
                                end
                            end
                        case 'radialPulsation' 
                            k = omega/layer{1}.c_f;
                            if layer{m}.calc_p_inc
                                layer{m}.p_inc = P_inc.*exp(1i*k.*(R-layer{1}.R))./R;
                            end
                            for ii = 1:3
                                if layer{m}.calc_dp_inc(ii)
                                    layer{m}.dp_inc{ii} = p_inc.*(1i*k - 1./R).*layer{m}.X(:,ii).'./R;
                                end
                            end
                        case 'surfExcitation'
                            if layer{m}.calc_p_inc
                                layer{m}.p_inc = zeros(nFreqs,n_X,prec);
                                layer{m}.p_inc(and(and(R == options.r_s,options.theta_s(1) < Theta), Theta > options.theta_s(2))) = 1;
                            end
                            for ii = 1:3
                                if layer{m}.calc_dp_inc(ii)
                                    layer{m}.dp_inc{ii} = zeros(nFreqs,n_X,prec);
                                end
                            end
                        otherwise
                            if any(layer{m}.calc_dp_inc) || layer{m}.calc_p_inc
                                error('Not valid')
                            end
                    end
                else
                    if layer{m}.calc_p_inc
                        layer{m}.p_inc = zeros(nFreqs,n_X,prec);
                    end
                    for ii = 1:3
                        if layer{m}.calc_dp_inc(ii)
                            layer{m}.dp_inc{ii} = zeros(nFreqs,n_X,prec);
                        end
                    end
                end
            case {'solid','viscoelastic'}
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
                    layer{m}.u = cell(3,1);
                    for ii = 1:3
                        if layer{m}.calc_u(ii)
                        	layer{m}.u{ii} = zeros(nFreqs,n_X,prec);
                            for jj = 1:3
                                layer{m}.u{ii} = layer{m}.u{ii} + A(ii,jj)*u_X_m{jj};
                            end
                        end
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
                            if layer{m}.calc_du(ii,jj)
                                layer{m}.du{ii,jj} = zeros(nFreqs,n_X,prec);
                                for ll = 1:3
                                    layer{m}.du{ii,jj} = layer{m}.du{ii,jj} + temp{ii, ll}*Ainv(ll,jj);
                                end
                            end
                        end
                    end
                    
                    if strcmp(layer{m}.media,'viscoelastic') && layer{m}.calc_p
                        Omega = repmat(omega,1,n_X);
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
                        c = repmat(4*1i*wv*muv/(3*rho) + ev.^2./kL.^2,1,n_X);

                        layer{m}.p = -1i*rho*c^2./Omega.*(layer{m}.du{1,1}+layer{m}.du{2,2}+layer{m}.du{3,3});
                    end
                end
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
                        if layer{m}.calc_sigma(vgtIdx)
                            layer{m}.sigma{vgtIdx} = zeros(nFreqs,n_X,prec);
                            for ii = 1:3
                                for jj = 1:3
                                    layer{m}.sigma{vgtIdx} = layer{m}.sigma{vgtIdx} + alpha(vgt(vgtIdx,1),ii)*alpha(vgt(vgtIdx,2),jj)*sigma_X_m{vgtinv(ii,jj)};
                                end
                            end
                        end
                    end
                    if layer{m}.calc_sigma_s(1)
                        layer{m}.sigma_s{1} = layer{m}.sigma_rr;
                    end
                end
        end
    end
end

% Calculate residual errors
SHBC = options.SHBC;
SSBC = options.SSBC;
if M > 1
    nextMedia = layer{2}.media;
else
    nextMedia = 'void';
end
for m = 1:M
    if layer{m}.calc_errors
        m2 = m + 1;
        Eps = options.Eps;
        R = layer{m}.R;
        isSphere = R == 0;
        n_X = size(layer{m}.X,1);
        X = layer{m}.X;
        
        % Calculate volumetric PDE residual errors
        switch layer{m}.media
            case 'fluid'
                if layer{m}.calc_err_helmholtz
                    k = omega/layer{m}.c_f;
                    K = repmat(k,1,n_X);

                    layer{m}.err_helmholtz = max(abs(layer{m}.p_laplace+K.^2.*layer{m}.p),[],2)./max(abs(K.^2.*layer{m}.p),[],2);
                end
            case {'solid','viscoelastic'}
                if any(layer{m}.calc_err_navier)
                    Theta = repmat(layer{m}.theta,nFreqs,1);
                    rho = layer{m}.rho;
                    u_r = layer{m}.u_r;
                    u_t = layer{m}.u_t.*sin(Theta); % rescale u_t
                    Omega = repmat(omega,1,n_X);
                    layer{m}.err_navier = cell(1,2);
                    if layer{m}.calc_err_navier(1)
                        layer{m}.err_navier{1} = max(abs(layer{m}.navier1 + rho*Omega.^2.*u_r),[],2)./max(abs(rho*Omega.^2.*u_r),[],2);
                    end
                    if layer{m}.calc_err_navier(2)
                        layer{m}.err_navier{2} = max(abs(layer{m}.navier2 + rho*Omega.^2.*u_t),[],2)./max(abs(rho*Omega.^2.*u_t),[],2);
                    end
                end
        end
        
        % Calculate surface BC residual errors
        P_inc = options.P_inc;
        if layer{m}.calc_err_pc && ~isSphere
            if m < M
                X2 = layer{m+1}.X;
                [indices, indices2] = findMatchingPoints(X,X2,Eps);
            else
                indices = abs(norm2(X) - R) < 10*Eps;
            end
            X_s = X(indices,:);
            switch layer{m}.media
                case 'fluid'
                    p_tot = layer{m}.p(:,indices);
                    if m == m_s
                        switch options.applyLoad
                            case 'planeWave'
                                p_tot = p_tot + layer{m}.p_inc(:,indices);
                            otherwise
                                error('Not implemented')
                        end
                    end
                    switch nextMedia
                        case {'solid','viscoelastic'}
                            sigma_cart = cell(6,1);
                            for i = 1:6
                                sigma_cart{i} = layer{m+1}.sigma{i}(:,indices2);
                            end

                            sigma_rr = zeros(size(sigma_cart{1}),prec);

                            phi_temp = atan2(X_s(:,2),X_s(:,1));
                            r_temp = sqrt(X_s(:,1).^2+X_s(:,2).^2+X_s(:,3).^2);
                            theta_temp = acos(X_s(:,3)./r_temp);
                            D = getStressTransformationMatrix(theta_temp,phi_temp,1);
                            for l = 1:6
                                D_kl = D(1, l, :);
                                D_kl = repmat(reshape(D_kl,1,[]),length(omega),1);
                                sigma_rr = sigma_rr + D_kl.*sigma_cart{l};
                            end
                            layer{m}.err_pc = max(abs(p_tot+sigma_rr),[],2)./max(abs(p_tot),[],2);
                        case 'fluid'
                            p_tot2 = layer{m2}.p(:,indices2);
                            if m2 == m_s
                                switch options.applyLoad
                                    case 'planeWave'
                                        p_tot2 = p_tot2 + layer{m2}.p_inc(:,indices);
                                    otherwise
                                        error('Not implemented')
                                end
                            end
                            layer{m}.err_pc = max(abs(p_tot-p_tot2),[],2)./max(abs(p_tot),[],2);
                    end
                case {'solid','viscoelastic'}
                    sigma_cart = cell(6,1);
                    for i = 1:6
                        sigma_cart{i} = layer{m}.sigma{i}(:,indices);
                    end

                    sigma_rr = zeros(size(sigma_cart{1}),prec);

                    phi_temp = atan2(X_s(:,2),X_s(:,1));
                    r_temp = sqrt(X_s(:,1).^2+X_s(:,2).^2+X_s(:,3).^2);
                    theta_temp = acos(X_s(:,3)./r_temp);
                    D = getStressTransformationMatrix(theta_temp,phi_temp,1);
                    for l = 1:6
                        D_kl = D(1, l, :);
                        D_kl = repmat(reshape(D_kl,1,[]),length(omega),1);
                        sigma_rr = sigma_rr + D_kl.*sigma_cart{l};
                    end
                    switch nextMedia
                        case 'void'
                            if SSBC
                                layer{m}.err_pc = max(abs(sigma_rr),[],2)/P_inc;
                            end
                        case 'fluid'
                            p_tot = layer{m+1}.p(:,indices2);
                            if m2 == m_s
                                k = omega/layer{m+1}.c_f;
                                k_vec = options.d_vec*k;
                                p_inc = P_inc*exp(1i*X_s*k_vec);
                                p_tot = p_tot + p_inc;
                            end
                            layer{m}.err_pc = max(abs(p_tot+sigma_rr),[],2)./max(abs(p_tot),[],2);
                        case {'solid','viscoelastic'}
                            sigma_cart2 = cell(6,1);
                            for i = 1:6
                                sigma_cart2{i} = layer{m2}.sigma{i}(:,indices2);
                            end

                            sigma_rr2 = zeros(size(sigma_cart2{1}),prec);

                            phi_temp = atan2(X_s(:,2),X_s(:,1));
                            r_temp = sqrt(X_s(:,1).^2+X_s(:,2).^2+X_s(:,3).^2);
                            theta_temp = acos(X_s(:,3)./r_temp);
                            D = getStressTransformationMatrix(theta_temp,phi_temp,1);
                            for l = 1:6
                                D_kl = D(1, l, :);
                                D_kl = repmat(reshape(D_kl,1,[]),length(omega),1);
                                sigma_rr2 = sigma_rr2 + D_kl.*sigma_cart2{l};
                            end
                            layer{m}.err_pc = max(abs(sigma_rr-sigma_rr2),[],2)./max(abs(sigma_rr),[],2);
                    end
            end
        end
        if layer{m}.calc_err_dc && ~isSphere
            if m < M
                X2 = layer{m+1}.X;
                [indices, indices2] = findMatchingPoints(X,X2,Eps);
            else
                indices = abs(norm2(X) - R) < 10*Eps;
            end
            switch layer{m}.media
                case 'fluid'
                    P_inc = options.P_inc;

                    rho = layer{m}.rho;

                    n_x = X_s(:,1)./norm2(X_s);
                    n_y = X_s(:,2)./norm2(X_s);
                    n_z = X_s(:,3)./norm2(X_s);
                    n_x = repmat(n_x,1,length(omega)).';
                    n_y = repmat(n_y,1,length(omega)).';
                    n_z = repmat(n_z,1,length(omega)).';

                    dp_s = cell(1,3);
                    for i = 1:3
                        dp_s{i} = layer{m}.dp{i}(:,indices);
                        if i == 1
                            Omega = repmat(omega,1,size(dp_s{i},2));
                        end
                        if m == m_s
                            dp_s{i} = dp_s{i} + layer{m}.dp_inc{i}(:,indices);
                        end
                        dp_s{i} = dp_s{i}./(rho*Omega.^2);
                    end
                    switch nextMedia
                        case 'void'
                            if SHBC
                                layer{m}.err_dc = max(abs(dp_s{1}.*n_x + dp_s{2}.*n_y + dp_s{3}.*n_z),[],2)/P_inc;
                            end
                        case 'fluid'
                            dp_s2 = cell(1,3);
                            rho2 = layer{m2}.rho;
                            for i = 1:3
                                dp_s2{i} = layer{m2}.dp{i}(:,indices2);
                                if m2 == m_s
                                    dp_s2{i} = dp_s2{i} + layer{m2}.dp_inc{i}(:,indices2);
                                end
                                dp_s2{i} = dp_s2{i}./(rho2*Omega.^2);
                            end
                            layer{m}.err_dc = max(abs((dp_s{1}-dp_s2{1}).*n_x + (dp_s{2}-dp_s2{2}).*n_y + (dp_s{3}-dp_s2{3}).*n_z),[],2)./max(abs(dp_s{1}.*n_x + dp_s{2}.*n_y + dp_s{3}.*n_z),[],1);                            
                        case {'solid','viscoelastic'}
                            u = cell(1,3);
                            for i = 1:3
                                u{i} = layer{m+1}.u{i}(:,indices2);
                            end
                            layer{m}.err_dc = max(abs(    (dp_s{1}-u{1}).*n_x ...
                                                        + (dp_s{2}-u{2}).*n_y ...
                                                        + (dp_s{3}-u{3}).*n_z),[],2)./max(abs(dp_s{1}.*n_x + dp_s{2}.*n_y + dp_s{3}.*n_z),[],1);
                    end
                case {'solid','viscoelastic'}
                    u = cell(1,3);
                    for i = 1:3
                        u{i} = layer{m}.u{i}(:,indices);
                    end
                    Omega = repmat(omega,1,size(u{1},2));
                    switch nextMedia
                        case 'fluid'
                            n_x = X_s(:,1)./norm2(X_s);
                            n_y = X_s(:,2)./norm2(X_s);
                            n_z = X_s(:,3)./norm2(X_s);
                            n_x = repmat(n_x,1,length(omega)).';
                            n_y = repmat(n_y,1,length(omega)).';
                            n_z = repmat(n_z,1,length(omega)).';
                            rho = layer{m+1}.rho;

                            dpdx = layer{m+1}.dp{1}(:,indices2);
                            dpdy = layer{m+1}.dp{2}(:,indices2);
                            dpdz = layer{m+1}.dp{3}(:,indices2);
                            layer{m}.err_dc = max(abs(    (dpdx-rho*Omega.^2.*u{1}).*n_x ...
                                                        + (dpdy-rho*Omega.^2.*u{2}).*n_y ...
                                                        + (dpdz-rho*Omega.^2.*u{3}).*n_z),[],2)./max(abs(dpdx.*n_x + dpdy.*n_y + dpdz.*n_z),[],2);
                        case {'solid','viscoelastic'}
                            u2 = cell(1,3);
                            for i = 1:3
                                u2{i} = layer{m2}.u{i}(:,indices2);
                            end
                            layer{m}.err_dc = max(abs(    (u{1}-u2{1}).*n_x ...
                                                        + (u{2}-u2{2}).*n_y ...
                                                        + (u{3}-u2{3}).*n_z),[],2)./max(abs(u{1}.*n_x + u{2}.*n_y + u{3}.*n_z),[],2);
                            
                    end
            end
        end
    end
    if M > 1 && ~(strcmp(nextMedia,'void') || strcmp(nextMedia,'origin'))
        if m+1 == M
            if layer{end}.R == 0
                nextMedia = 'origin';
            else
                nextMedia = 'void';
            end
        else
            nextMedia = layer{m+2}.media;
        end
    end
end

% Remove temporary fields
for m = 1:M
    for i = 1:numel(fieldsToRemove)
        if isfield(layer{m},fieldsToRemove{i}) && ~fieldExist{m}(i) && ~(isfield(layer{m},['calc_' fieldsToRemove{i}]) && layer{m}.(['calc_' fieldsToRemove{i}]))
            layer{m} = rmfield(layer{m},fieldsToRemove{i});
        end
    end
end

% Transpose back fields
for m = 1:M
    for field = fields(layer{m}).'
        fieldStr = field{1};
        if numel(fieldStr) > 5 && strcmp(fieldStr(1:5),'calc_')
            fieldName = fieldStr(6:end);
            d_f = numel(layer{m}.(fieldStr));
            for i = 1:d_f
                if layer{m}.(fieldStr)(i)
                    if isfield(layer{m},fieldName)
                        if d_f == 1 % Scalar fields are not represented as cell arrays
                            layer{m}.(fieldName) = layer{m}.(fieldName).';
                        else
                            layer{m}.(fieldName){i} = layer{m}.(fieldName){i}.';
                        end
                    end
                end
            end
        end
    end
end

function layer = getDefaultParameters(layer)

for m = 1:numel(layer)
    layer{m}.R         = 1;       % Inner radius of layer
    layer{m}.rho         = 1000;    % Mass density
    layer{m}.calc_err_dc = false;   % Calculate the errors for the displacement conditions
    layer{m}.calc_err_pc = false;   % Calculate the errors for the pressure conditions
    switch layer{m}.media
        case 'fluid'
            layer{m}.c_f          	    = 1500;       % Speed of sound
            layer{m}.calc_p_0           = false;      % Toggle calculation of the far field pattern
            layer{m}.calc_p       	    = false;      % Toggle calculation of the scattered pressure
            layer{m}.calc_dp      	    = false(1,3); % Toggle calculation of the three components of the gradient of the pressure
            layer{m}.calc_p_laplace	    = false;      % Toggle calculation of the Laplace operator of the scattered pressure fields
            layer{m}.calc_err_helmholtz	= false;      % Toggle calculation of the errors for the Helmholtz equation
            layer{m}.calc_p_inc         = false;      % Toggle calculation of the incident pressure
            layer{m}.lossFactor         = 0;          % Hysteretic loss factor (values around 0.001 for lightly damped materials, values around 0.01 for moderately damped materials and values around 0.1 for heavily damped materials)
            layer{m}.calc_dp_inc        = false(1,3); % Toggle calculation of the three components of the gradient of the incident pressure
        case {'solid','viscoelastic'}
            layer{m}.E                = 200e9;      % Youngs modulus for solid layers
            layer{m}.nu               = 0.3;        % Poisson ratio for solid layers
            layer{m}.calc_u           = false(1,3); % Toggle calculation of the three components of the displacement
            layer{m}.calc_du          = false(3,3); % Toggle calculation of the three cartesian derivatives of the three components of the displacement [du_xdx du_xdy du_xdz; 
                                                    %                                                                                                    du_ydx du_ydy du_ydz; 
                                                    %                                                                                                    du_zdx du_zdy du_zdz]
            layer{m}.calc_sigma       = false(1,6); % Toggle calculation of the six components of the stress field (cartesian coordinates) [sigma_xx sigma_yy sigma_zz sigma_yz sigma_xz sigma_xy]
            layer{m}.calc_sigma_s     = false(1,6); % Toggle calculation of the six components of the stress field (spherical coordinates) [sigma_rr sigma_tt sigma_pp sigma_tp sigma_rp sigma_rt]
            layer{m}.calc_err_navier  = false(1,2); % Toggle calculation of the errors for the Navier equation
            layer{m}.lossFactor       = [0,0];      % Hysteretic loss factor (first component for longitudinal velocities and second component for transverse velocities)
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [layer,N_eps,flag,relTermMaxArr] = e3Dss_0(layer, options)
% This function computes exact 3D scattering solutions when the axis of
% symmetry is the z-axis
%
%
% Note that dpdt, u_t, du_tdr, and du_rdt are scaled by csc(theta)
fluidFieldNames = {'p_0','p','dpdr','dpdt','d2pdr2','d2pdt2','p_inc','dp_incdr','dp_incdt'};
solidFieldNames = {'u_r','u_t','du_rdr','du_rdt','du_tdr','du_tdt','sigma_rr','sigma_tt','sigma_pp','sigma_rt','navier1','navier2'};

M = numel(layer);
m_s = options.m_s;
displayIter = strcmp(options.Display, 'iter');
prec = options.prec;
omega = options.omega;
nFreqs = length(omega);
Eps = options.Eps;
N_max = options.N_max;
nu_a = options.nu_a;
if any(omega == 0)
    computeForStaticCase = true;
else
    computeForStaticCase = false;
end
load(['miscellaneous/U_pol_' prec '.mat'],'U_pol','u_k','v_k')

R = inf;
for m = 1:M
    R_o = R;
    R = layer{m}.R;
    isOuterDomain = m == 1;
    isSphere = R == 0;
    if isfield(layer{m} ,'X')
        n_X = size(layer{m}.X,1);
    else
        n_X = 0;
    end
    layer{m}.P = zeros(2,n_X,prec); 
    layer{m}.dP = zeros(2,n_X,prec); 
    layer{m}.d2P = zeros(2,n_X,prec);
    besselIndices = getBesselFlags(layer{m},options,m,m_s,isOuterDomain,isSphere,true);
    switch layer{m}.media
        case 'fluid'
            % Modify sound speed to acount for the hysteretic loss factor
            c_f = layer{m}.c_f.*sqrt(1-1i*layer{m}.lossFactor(:,1));
            
            layer{m}.k = omega./c_f;
            
            for i = 1:numel(besselIndices)
                if besselIndices(i)
                    for j = 1:2
                        layer{m}.Z_zeta{i,j} = zeros(nFreqs-computeForStaticCase,n_X,prec);
                        if ~isSphere
                            layer{m}.Z_zeta_i{i,j} = zeros(nFreqs,1,prec);
                            layer{m}.Z_zeta_i_s{i,j} = zeros(nFreqs,1,prec);
                        end
                        if ~isinf(R_o)
                            layer{m}.Z_zeta_o{i,j} = zeros(nFreqs,1,prec);
                            layer{m}.Z_zeta_o_s{i,j} = zeros(nFreqs,1,prec);
                        end
                        if strcmp(options.applyLoad,'pointCharge')
                            layer{m}.Z_r_s{i,j} = zeros(nFreqs,1,prec);
                        end
                    end
                end
            end
            for fieldName = fluidFieldNames
                if layer{m}.(['calc_' fieldName{1}]) && ~(strcmp(fieldName{1},'p_inc') && ~options.p_inc_fromSeries)
                    layer{m}.(fieldName{1}) = zeros(nFreqs,n_X,prec);
                end
            end
        case {'solid','viscoelastic'}
            E = layer{m}.E;
            nu = layer{m}.nu;
            rho = layer{m}.rho;
            
            % Compute derived quantities            
            K = E./(3*(1-2*nu));
            G = E./(2*(1+nu));

            c_s_1 = sqrt((3*K+4*G)./(3*rho)); % longitudinal wave velocity 
            c_s_2 = sqrt(G./rho); % shear wave velocity
            
            % Add hysteris damping
            c_s_1 = c_s_1.*sqrt(1-1i*layer{m}.lossFactor(:,1));
            c_s_2 = c_s_2.*sqrt(1-1i*layer{m}.lossFactor(:,2));
            
            G = rho.*c_s_2.^2;
            K = rho.*c_s_1.^2 - 4*G/3;
            
            layer{m}.K = K;
            layer{m}.G = G;
            layer{m}.a = omega./c_s_1;
            layer{m}.b = omega./c_s_2;

            for i = 1:numel(besselIndices)
                if besselIndices(i)
                    for j = 1:2
                        layer{m}.Z_xi{i,j} = zeros(nFreqs-computeForStaticCase,n_X,prec);
                        layer{m}.Z_eta{i,j} = zeros(nFreqs-computeForStaticCase,n_X,prec);
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
            end
            for fieldName = solidFieldNames
                if layer{m}.(['calc_' fieldName{1}]) && ~(strcmp(fieldName{1},'p_inc') && ~options.p_inc_fromSeries)
                    layer{m}.(fieldName{1}) = zeros(nFreqs,n_X,prec);
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
hasDvrgd = zeros(nFreqs,1); % matrix of element that "has diverged"
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
relTermMaxFinal = zeros(nFreqs,1,prec); % Program terminated successfully unless error occurs (for each frequency)
if options.saveRelTermMax
    if isinf(N_max) % Track the relative term magnitude
        relTermMaxArr = zeros(nFreqs,1000,prec);
    else
        relTermMaxArr = zeros(nFreqs,N_max,prec);
    end
else
    relTermMaxArr = NaN(1,prec);
end
PI = getC(prec,'pi');
singleModeSolution = (strcmp(options.applyLoad,'pointCharge') && options.r_s == 0) || strcmp(options.applyLoad,'radialPulsation') ...
                     || (strcmp(options.applyLoad,'surfExcitation') && options.theta_s(1) == 0 && abs(options.theta_s(2) - PI) < Eps);

while n <= N_max && ~(singleModeSolution && n > 0)
    if displayIter
        tic
    end
    omega_temp = omega(indices);
    options.z_temp = getSubArray(options.z,indices);
    R = inf;
    for m = 1:M
        R_o = R;
        R = layer{m}.R;
        isSphere = R == 0;
        isOuterDomain = m == 1;
        besselIndices = getBesselFlags(layer{m},options,m,m_s,isOuterDomain,isSphere,true);
        switch layer{m}.media
            case 'fluid'
                evalPointCharge = strcmp(options.applyLoad,'pointCharge') && m == m_s;
                layer{m}.k_temp = layer{m}.k(indices,:);
                if ~isSphere
                    zeta_i = layer{m}.k_temp*R;
                    layer{m}.Z_zeta_i = iterate_Z(n,zeta_i,layer{m}.Z_zeta_i,Zindices,besselIndices,nu_a,U_pol,u_k,v_k,Eps);
                end
                if ~isOuterDomain
                    zeta_o = layer{m}.k_temp*R_o;
                    layer{m}.Z_zeta_o = iterate_Z(n,zeta_o,layer{m}.Z_zeta_o,Zindices,besselIndices,nu_a,U_pol,u_k,v_k,Eps);
                end
                if evalPointCharge
                    layer{m}.Z_r_s = iterate_Z(n,layer{m}.k_temp*options.r_s,layer{m}.Z_r_s,Zindices,besselIndices,nu_a,U_pol,u_k,v_k,Eps);
                end
            case {'solid','viscoelastic'}
                layer{m}.a_temp = layer{m}.a(indices,:);
                layer{m}.b_temp = layer{m}.b(indices,:);
                layer{m}.G_temp = getSubArray(layer{m}.G,indices);
                layer{m}.K_temp = getSubArray(layer{m}.K,indices);
                if ~isSphere
                    xi_i = layer{m}.a_temp*R;
                    eta_i = layer{m}.b_temp*R;
                    layer{m}.Z_xi_i  = iterate_Z(n,xi_i, layer{m}.Z_xi_i, Zindices,besselIndices,nu_a,U_pol,u_k,v_k,Eps);
                    layer{m}.Z_eta_i = iterate_Z(n,eta_i,layer{m}.Z_eta_i,Zindices,besselIndices,nu_a,U_pol,u_k,v_k,Eps);
                end
                if ~isinf(R_o)
                    xi_o = layer{m}.a_temp*R_o;
                    eta_o = layer{m}.b_temp*R_o;
                    layer{m}.Z_xi_o  = iterate_Z(n,xi_o,layer{m}.Z_xi_o, Zindices,besselIndices,nu_a,U_pol,u_k,v_k,Eps);
                    layer{m}.Z_eta_o = iterate_Z(n,eta_o,layer{m}.Z_eta_o,Zindices,besselIndices,nu_a,U_pol,u_k,v_k,Eps);
                end
        end
        % Compute exponent of scale for incident wave
        if m == m_s
            switch options.applyLoad
                case 'planeWave'
                    nu = n+0.5;
                    Rt_m = getIntermediateRadius(layer,m,isSphere,isOuterDomain);
                    zetat_m = Rt_m*layer{m}.k_temp;
                    exponentShift = exponent_(1,nu,zetat_m,nu_a);
                case 'pointCharge'
                    nu = n+0.5;
                    Rt_m = getIntermediateRadius(layer,m,isSphere,isOuterDomain);
                    zetat_m = Rt_m*layer{m}.k_temp;
                    r_s = options.r_s;
                    zetat_s = r_s*layer{m}.k_temp;
                    if Rt_m < r_s
                        exponentShift = exponent_(1,nu,zetat_m,nu_a) + exponent_(3,nu,zetat_s,nu_a);
                    else
                        exponentShift = exponent_(1,nu,zetat_s,nu_a) + exponent_(3,nu,zetat_m,nu_a);
                    end
                otherwise
                    exponentShift = 0;
            end
        end
    end
    C = getCoeffs(n, omega_temp, layer, options);
    hasCnvrgdTmp = zeros(length(indices),M); % temporary hasCnvrgd matrix
    hasDvrgdTmp = zeros(length(indices),M); % temporary hasDvrgd matrix
    relTermMax = -Inf(nFreqs,1,prec);
    for m = 1:M
        isSphere = layer{m}.R == 0;
        isOuterDomain = m == 1;
        besselIndices = getBesselFlags(layer{m},options,m,m_s,isOuterDomain,isSphere);
        hasCnvrgdTmp2 = ones(size(indices)); % temporary hasCnvrgd vector  
        hasDvrgdTmp2 = zeros(size(indices)); % temporary hasDvrgd vector    
        if isfield(layer{m},'r') && ~isempty(layer{m}.r)
            [layer{m}.P, layer{m}.dP, layer{m}.d2P] = legendreDerivs(n, cos(layer{m}.theta), layer{m}.P, layer{m}.dP, layer{m}.d2P);
            Rt_m = getIntermediateRadius(layer,m,isSphere,isOuterDomain);

            switch layer{m}.media
                case 'fluid'
                    zeta = layer{m}.k_temp*layer{m}.r;
                    if ~layer{m}.calc_farFieldOnly
                        layer{m}.Z_zeta = iterate_Z(n,zeta,layer{m}.Z_zeta,Zindices,besselIndices,nu_a,U_pol,u_k,v_k,Eps);
                    end
                    media = p_(m,n,zeta,C{m},Rt_m,layer,isSphere,isOuterDomain,options.applyLoad,besselIndices,nu_a,exponentShift,options);
                    fieldNames = fluidFieldNames;
                case {'solid','viscoelastic'}
                    xi = layer{m}.a_temp*layer{m}.r;
                    eta = layer{m}.b_temp*layer{m}.r;
                    layer{m}.Z_xi = iterate_Z(n,xi,layer{m}.Z_xi,Zindices,besselIndices,nu_a,U_pol,u_k,v_k,Eps);
                    layer{m}.Z_eta = iterate_Z(n,eta,layer{m}.Z_eta,Zindices,besselIndices,nu_a,U_pol,u_k,v_k,Eps);
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
                    media = u_(n,m,A,B,xi,eta,Rt_m,isSphere,layer,omega_temp,nu_a,exponentShift);
                    fieldNames = solidFieldNames;
            end
            
            % Update fields
            for fieldName = fieldNames
                if layer{m}.(['calc_' fieldName{1}]) && ~(strcmp(fieldName{1},'p_inc') && ~options.p_inc_fromSeries)
                    [layer,hasCnvrgdTmp2,hasDvrgdTmp2,relTermMax] = updateSum(layer,m,media,fieldName{1},indices,hasCnvrgdTmp2,hasDvrgdTmp2,relTermMax,tiny,Eps,prec);
                end
            end
        end
        hasCnvrgdTmp(logical(hasCnvrgdTmp2),m) = 1;
        hasDvrgdTmp(logical(hasDvrgdTmp2),m) = 1;
%         hasDvrgdTmp(any(C{m} == 0,2),m) = 1; % Assume that solution contains 1/Inf = 0 calculations (and flag it as divergent)
    end
    if options.saveRelTermMax
        relTermMaxArr(:,n+1) = relTermMax;
    end
    hasCnvrgd(indices,:) = [hasCnvrgd(indices,2:end), prod(hasCnvrgdTmp,2)];
    hasDvrgd(indices) = sum(hasDvrgdTmp,2);
    indicesPrev = indices;
    indices = find(and(~prod(hasCnvrgd,2),~hasDvrgd)); 
    relTermMaxFinal(indices) = relTermMax(indices);  % Update the final relTerm
    [~,Zindices] = ismember(indices,indicesPrev);
    if length(indices) < length(indicesPrev)
        N_eps(setdiff(indicesPrev,indices)) = n;
    end
    if isempty(indices) % every element has converged/diverged
        break;
    end
    if displayIter
        fprintf('n = %5d, Relative term: %10.5g, Elapsed time: %8.6g seconds.\n', n, double(max(relTermMax)), toc)
    end
    n = n + 1;
end
if singleModeSolution
    flag = -hasDvrgd;
else
    flag = -or(~hasCnvrgd(:,end),hasDvrgd);
end
if options.saveRelTermMax
    relTermMaxArr = relTermMaxArr(:,1:n-1);
end
if strcmp(options.Display, 'final') || strcmp(options.Display, 'iter')
    if n-1 == N_max
        warning('e3Dss:N_max_reached','The summation did not converge using N_max = %d terms (relTermMax = %f).', N_max, max(relTermMaxFinal))
    elseif ~any(flag)
        if singleModeSolution
            fprintf('Eps precision reached with a single term (N_eps = 1).\n')
        else
            fprintf('Eps precision reached with N_eps = %g terms.\n', max(N_eps))
        end
    else
        warning('e3Dss:diverged',['The summation ended prematurely at n = %d because some frequencies gave Inf/NaN computations (relTermMax = %g).\n' ...
                                  'Check the flag output variable to see what frequencies this applied to.'], n, max(relTermMaxFinal))
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
        if layer{m}.(['calc_' fieldName{1}]) && ~(strcmp(fieldName{1},'p_inc') && ~options.p_inc_fromSeries)
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
function x_temp = getSubArray(x,indices)
if numel(x) > 1
    x_temp = x(indices,:);
else
    x_temp = x;
end

function besselFlag = getBesselFlags(layer,options,m,m_s,isOuterDomain,isSphere,calcForCoeffs)
if nargin < 7
    calcForCoeffs = false;
end
besselFlag = false(1,3);
if isOuterDomain
    besselFlag(3) = true;
elseif isSphere
    besselFlag(1) = true;
else
    besselFlag(1:2) = true;
end
if strcmp(options.applyLoad,'planeWave') || strcmp(options.applyLoad,'pointCharge')
    calcForP_inc = (m == m_s && (calcForCoeffs || layer.calc_p_inc || layer.calc_dp_incdr || layer.calc_dp_incdt));
else
    calcForP_inc = false;
end
if strcmp(options.applyLoad,'planeWave')
    calcForP_inc = calcForP_inc*[1,0,0];
elseif strcmp(options.applyLoad,'pointCharge')
    calcForP_inc = calcForP_inc*[1,0,1];
end
besselFlag = or(besselFlag, calcForP_inc);

function [layer,hasCnvrgd,hasDvrgd,relTermMax] = updateSum(layer,m,media,fieldName,indices,hasCnvrgd,hasDvrgd,relTermMax,tiny,Eps,prec)
% Track any Inf/NaN in next term
hasDvrgdSub = or(hasDvrgd,or(any(isinf(media.(fieldName)),2),any(isnan(media.(fieldName)),2)));

% Add non-Inf/NaN terms to series
layer{m}.(fieldName)(indices(~hasDvrgdSub),:) = layer{m}.(fieldName)(indices(~hasDvrgdSub),:) + media.(fieldName)(~hasDvrgdSub,:);

% Compute relative error
relTerm = abs(media.(fieldName))./(abs(layer{m}.(fieldName)(indices,:))+tiny);

% Update converged/diverged tracking arrays
hasCnvrgd = hasCnvrgd.*prod(double(relTerm < Eps), 2);
hasDvrgd = hasDvrgd + hasDvrgdSub;

% Only track frequencies that have not converged
nonConvergedIdx = false(size(relTermMax,1),1);
nonConvergedIdx(indices) = true;
maxRelTerm = -Inf(size(relTermMax,1),1,prec);
maxRelTerm(nonConvergedIdx) = max(relTerm,[],2);
indices2 = logical(double(maxRelTerm > relTermMax));
relTermMax(indices2) = maxRelTerm(indices2);

function Rt_m = getIntermediateRadius(layer,m,isSphere,isOuterDomain)
if isOuterDomain
    Rt_m = layer{1}.R;
elseif isSphere
    Rt_m = layer{m-1}.R;
else
    Rt_m = (layer{m}.R+layer{m-1}.R)/2;
end

function fluid = p_(m,n,zeta,C,Rt_m,layer,isSphere,isOuterDomain,applyLoad,besselIndices,nu_a,exponentShift,options)
% Note that in the case of isSphere and zeta = 0: 
% --- dpdz =: dpdr and dpdx = dpdy = dpdt = 0, with same convention for p_inc
% --- nabla p =: d2pdr2, d2pdt2 := 0
% Also note that dpdt and dp_incdt are scaled by csc(theta)

fluid = [];
k = layer{m}.k_temp;
P = layer{m}.P(2,:);
dP = layer{m}.dP(2,:);
d2P = layer{m}.d2P(2,:);
Z = layer{m}.Z_zeta;
g = {g_(n,1,zeta,nu_a),g_(n,2,zeta,nu_a),g_(n,3,zeta,nu_a)};
r = layer{m}.r;
theta = layer{m}.theta;

Q0 = P;
if layer{m}.calc_dpdt
    Q1 = Q_(1,theta,P,dP,d2P,true);
end
if layer{m}.calc_d2pdt2
    Q2 = Q_(2,theta,P,dP,d2P);
end
zetat_m = Rt_m*k;
nu = n + 1/2;
if layer{m}.calc_p || layer{m}.calc_dpdr || layer{m}.calc_dpdt || layer{m}.calc_d2pdr2 || layer{m}.calc_d2pdr2 ...
        || layer{m}.calc_p_inc || layer{m}.calc_dp_incdr || layer{m}.calc_dp_incdt
    w = cell(1,3);
    wZt = cell(1,3);
    dwZt = cell(1,3);
    d2wZt = cell(1,3);
    for i = 1:numel(besselIndices)
        if besselIndices(i)
            wZt{i} = Z{i,1};
            if layer{m}.calc_dpdr
                dwZt{i} = dbessel_s(n,zeta,i,Z,false,g);
            end
            if layer{m}.calc_d2pdr2
                d2wZt{i} = d2bessel_s(n,zeta,i,Z,g);
            end
        end
    end
    if layer{m}.calc_p_inc && options.p_inc_fromSeries
        switch applyLoad
            case 'planeWave'
                s1 = exp(exponent_(1,nu,zeta,nu_a));
                fluid.p_inc = (2*n+1)*1i^n*Q0.*wZt{1}./s1;
            case 'pointCharge'
                s1 = exp(exponent_(1,nu,zeta,nu_a));
                s3 = exp(exponent_(3,nu,zeta,nu_a));
                Z_r_s = layer{m}.Z_r_s;
                r_s = options.r_s;
                s1_r_s = exp(exponent_(1,nu,k*r_s,nu_a));
                s3_r_s = exp(exponent_(3,nu,k*r_s,nu_a));
                fluid.p_inc = zeros(size(zeta),class(zeta));
                indices = repmat(r < r_s,numel(k),1);
                temp = (2*n+1)*1i*k*Q0.*wZt{1}.*Z_r_s{3,1}./s3_r_s./s1;
                fluid.p_inc(indices)  = temp(indices);
                temp = (2*n+1)*1i*k*Q0.*Z_r_s{1,1}.*wZt{3}./s1_r_s./s3;
                fluid.p_inc(~indices) = temp(~indices);
            otherwise
                fluid.p_inc = zeros(size(zeta),class(zeta));
        end
    end
    if layer{m}.calc_dp_incdr && options.p_inc_fromSeries
        switch applyLoad
            case 'planeWave'
                s1 = exp(exponent_(1,nu,zeta,nu_a));
                fluid.dp_incdr = (2*n+1)*1i^n*k*Q0.*dwZt{1}./s1;
            case 'pointCharge'
                s1 = exp(exponent_(1,nu,zeta,nu_a));
                s3 = exp(exponent_(3,nu,zeta,nu_a));
                Z_r_s = layer{m}.Z_r_s;
                r_s = options.r_s;
                s1_r_s = exp(exponent_(1,nu,k*r_s,nu_a));
                s3_r_s = exp(exponent_(3,nu,k*r_s,nu_a));
                fluid.dp_incdr = zeros(size(zeta),class(zeta));
                indices = repmat(r < r_s,numel(k),1);
                temp = (2*n+1)*1i*k*Q0.*dwZt{1}.*Z_r_s{3,1}./s3_r_s./s1;
                fluid.dp_incdr(indices)  = temp(indices);
                temp = (2*n+1)*1i*k*Q0.*Z_r_s{1,1}.*dwZt{3}./s1_r_s./s3;
                fluid.dp_incdr(~indices) = temp(~indices);
            otherwise
                fluid.dp_incdr = zeros(size(zeta),class(zeta));
        end
    end
    if layer{m}.calc_dp_incdt && options.p_inc_fromSeries
        switch applyLoad
            case 'planeWave'
                s1 = exp(exponent_(1,nu,zeta,nu_a));
                fluid.dp_incdt = (2*n+1)*1i^n*Q1.*wZt{1}./s1;
            case 'pointCharge'
                s1 = exp(exponent_(1,nu,zeta,nu_a));
                s3 = exp(exponent_(3,nu,zeta,nu_a));
                Z_r_s = layer{m}.Z_r_s;
                r_s = options.r_s;
                s1_r_s = exp(exponent_(1,nu,k*r_s,nu_a));
                s3_r_s = exp(exponent_(3,nu,k*r_s,nu_a));
                fluid.dp_incdt = zeros(size(zeta),class(zeta));
                indices = repmat(r < r_s,numel(k),1);
                temp = (2*n+1)*1i*k*Q1.*wZt{1}.*Z_r_s{3,1}./s3_r_s./s1;
                fluid.dp_incdt(indices)  = temp(indices);
                temp = (2*n+1)*1i*k*Q1.*Z_r_s{1,1}.*wZt{3}./s1_r_s./s3;
                fluid.dp_incdt(~indices) = temp(~indices);
            otherwise
                fluid.dp_incdt = zeros(size(zeta),class(zeta));
        end
    end
    for i = 1:numel(besselIndices)
        if besselIndices(i)
            w{i} = w_(n,i,zetat_m,zeta,nu_a,exponentShift);
            wZt{i} = w{i}.*wZt{i};
            if layer{m}.calc_dpdr
                dwZt{i} = w{i}.*dwZt{i};
            end
            if layer{m}.calc_d2pdr2
                d2wZt{i} = w{i}.*d2wZt{i};
            end
        end
    end
end

if isSphere
    indices = logical(zeta(1,:) < eps);
    sOrigin  = exp(exponent_(1,nu,zetat_m,nu_a)-exponentShift);
end
if isOuterDomain
    if layer{m}.calc_p_0
        cs = exp(exponent_(3,nu,zetat_m,nu_a)-exponentShift);
        if isinf(cs)
            warning('e3Dss:infWeight','A weight evaluation was too large')
        end
        h_n_0  = 1i^(-n-1)./k;
        fluid.p_0 = C*Q0.*h_n_0.*cs;
    end
%     if layer{m}.calc_dp_0dr
%         dh_n_0 = 1i^(-n);
%         fluid.dh_n_0  = C*Q0.*dh_n_0.*cs;
%     end
%     if layer{m}.calc_d2p_0dr2
%         d2h_n_0 = 1i^(-n+1)*k;
%         fluid.d2h_n_0 = C*Q0.*d2h_n_0.*cs;
%     end
    if layer{m}.calc_p
        fluid.p = C*Q0.*wZt{3};
    end
    if layer{m}.calc_dpdr
        fluid.dpdr = C.*k*Q0.*dwZt{3};
    end
    if layer{m}.calc_dpdt % calculate a scaled dpdt (scaled by csc)
        fluid.dpdt = C*Q1.*wZt{3};
    end
    if layer{m}.calc_d2pdr2
        fluid.d2pdr2 = C.*k.^2*Q0.*d2wZt{3};
    end
    if layer{m}.calc_d2pdt2
        fluid.d2pdt2 = C*Q2.*wZt{3};
    end
elseif isSphere
    if layer{m}.calc_p
        fluid.p = C*Q0.*wZt{1};
        if n == 0
            fluid.p(:,indices) = repmat(C,1,sum(indices)).*sOrigin;
        else
            fluid.p(:,indices) = 0;
        end
    end
    if layer{m}.calc_dpdr
        fluid.dpdr = C.*k*Q0.*dwZt{1};
        if n == 1
            fluid.dpdr(:,indices) = repmat(k/3.*C,1,sum(indices)).*sOrigin;
        else
            fluid.dpdr(:,indices) = 0;
        end
    end
    if layer{m}.calc_dpdt % calculate a scaled dpdt (scaled by csc)
        fluid.dpdt = C*Q1.*wZt{1};
        fluid.dpdt(:,indices) = 0;
    end
    if layer{m}.calc_d2pdr2
        fluid.d2pdr2 = C.*k.^2*Q0.*d2wZt{1};
        if n == 0
            fluid.d2pdr2(:,indices) = repmat(-k.^2.*C,1,sum(indices)).*sOrigin;
        else
            fluid.d2pdr2(:,indices) = 0;
        end
    end
    if layer{m}.calc_d2pdt2
        fluid.d2pdt2 = C*Q2.*wZt{1};
        fluid.d2pdt2(:,indices) = 0;
    end
else
    if layer{m}.calc_p
        fluid.p = C(:,1)*Q0.*wZt{1} + C(:,2)*Q0.*wZt{2};
    end
    if layer{m}.calc_dpdr
        fluid.dpdr = C(:,1).*k*Q0.*dwZt{1} + C(:,2).*k*Q0.*dwZt{2};
    end
    if layer{m}.calc_dpdt % calculate a scaled dpdt (scaled by csc)
        fluid.dpdt = C(:,1)*Q1.*wZt{1} + C(:,2)*Q1.*wZt{2};
    end
    if layer{m}.calc_d2pdr2
        fluid.d2pdr2 = C(:,1).*k.^2*Q0.*d2wZt{1} + C(:,2).*k.^2*Q0.*d2wZt{2};
    end
    if layer{m}.calc_d2pdt2
        fluid.d2pdt2 = C(:,1)*Q2.*wZt{1} + C(:,2)*Q2.*wZt{2};
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function solid = u_(n,m,A,B,xi,eta,Rt_m,isSphere,layer,omega,nu_a,exponentShift)
% Note that in the case of isSphere and r = 0: 
% --- u_z =: u_r and u_x = u_y = u_t = 0
% --- du_xdx =: du_rdr, du_ydy =: du_rdt, du_zdz =: du_tdr, 0 =: du_tdt
% --- sigma_11 =: sigma_rr, sigma_22 =: sigma_tt, sigma_33 =: sigma_pp, 0 =: sigma_rt
% Also note that u_t, du_tdr and du_rdt are scaled by csc(theta)

solid = [];
G = layer{m}.G_temp;
K = layer{m}.K_temp;
a = layer{m}.a_temp;
b = layer{m}.b_temp;
P = layer{m}.P(2,:);
dP = layer{m}.dP(2,:);
d2P = layer{m}.d2P(2,:);
Zxi = layer{m}.Z_xi;
Zeta = layer{m}.Z_eta;
r = layer{m}.r;
theta = layer{m}.theta;
gxi  = {g_(n,1,xi, nu_a), g_(n,2,xi, nu_a)};
geta = {g_(n,1,eta,nu_a), g_(n,2,eta,nu_a)};

Q0 = Q_(0,theta,P,dP,d2P);
Q1s = Q_(1,theta,P,dP,d2P,true);
Q1 = sin(theta).*Q1s;
Q2 = Q_(2,theta,P,dP,d2P);

r2 = r.^2;
xit_m = Rt_m*a;
etat_m = Rt_m*b;
nu = n + 1/2;
if isSphere
    indices = logical(r < eps);
    sOrigin_xi  = exp(exponent_(1,nu,xit_m, nu_a)-exponentShift);
    sOrigin_eta = exp(exponent_(1,nu,etat_m,nu_a)-exponentShift);
end
wjxi = w_(n,1,xit_m,xi,nu_a,exponentShift);
wjeta = w_(n,1,etat_m,eta,nu_a,exponentShift);
if ~isSphere
    wyxi = w_(n,2,xit_m,xi,nu_a,exponentShift);
    wyeta = w_(n,2,etat_m,eta,nu_a,exponentShift);
end
if layer{m}.calc_u_r
    Q0r = Q0./r;
    u_r = A(:,1)*Q0r.*wjxi.*S_(1,1,n,xi,eta,Zxi,gxi);
    if n > 0
    	u_r = u_r + B(:,1)*Q0r.*wjeta.*T_(1,1,n,eta,Zeta,geta);
    end
    if isSphere
        if n == 1
            u_r(:,indices) = repmat((a.*A(:,1).*sOrigin_xi - 2*b.*B(:,1).*sOrigin_eta)/3,1,sum(indices));
        else
            u_r(:,indices) = 0;
        end
    else
        u_r = u_r + A(:,2)*Q0r.*wyxi.*S_(1,2,n,xi,eta,Zxi,gxi);
        if n > 0
        	u_r = u_r + B(:,2)*Q0r.*wyeta.*T_(1,2,n,eta,Zeta,geta);
        end
    end
    solid.u_r = u_r;
end
if layer{m}.calc_u_t
    Q1sr = Q1s./r;
    u_t = A(:,1)*Q1sr.*wjxi.*S_(2,1,n,xi,eta,Zxi,gxi);
    if n > 0
    	u_t = u_t + B(:,1)*Q1sr.*wjeta.*T_(2,1,n,eta,Zeta,geta);
    end
    if isSphere
        u_t(:,logical(r < eps)) = 0;
    else
        u_t = u_t + A(:,2)*Q1sr.*wyxi.*S_(2,2,n,xi,eta,Zxi,gxi);
        if n > 0
            u_t = u_t + B(:,2)*Q1sr.*wyeta.*T_(2,2,n,eta,Zeta,geta);
        end
    end
    solid.u_t = u_t;
end

if layer{m}.calc_du_rdr
    Q0r2 = Q0./r2;
    du_rdr = A(:,1)*Q0r2.*wjxi.*S_(3,1,n,xi,eta,Zxi,gxi);
    if n > 0
        du_rdr = du_rdr + B(:,1)*Q0r2.*wjeta.*T_(3,1,n,eta,Zeta,geta);
    end
    if isSphere
        if n == 0
            du_rdr(:,indices) = repmat(G./K.*(4*a.^2-3*b.^2).*A(:,1).*sOrigin_xi/9,1,sum(indices));
        elseif n == 2
            du_rdr(:,indices) = repmat(-(a.^2.*A(:,1).*sOrigin_xi - 3*b.^2.*B(:,1).*sOrigin_eta)/15,1,sum(indices));
        else
            du_rdr(:,indices) = 0;
        end
    else
        du_rdr = du_rdr + A(:,2)*Q0r2.*wyxi.*S_(3,2,n,xi,eta,Zxi,gxi);
        if n > 0
            du_rdr = du_rdr + B(:,2)*Q0r2.*wyeta.*T_(3,2,n,eta,Zeta,geta);
        end
    end
    solid.du_rdr = du_rdr;
end

if layer{m}.calc_du_rdt
    Q1sr = Q1s./r;
    du_rdt = A(:,1)*Q1sr.*wjxi.*S_(1,1,n,xi,eta,Zxi,gxi);
    if n > 0
        du_rdt = du_rdt + B(:,1)*Q1sr.*wjeta.*T_(1,1,n,eta,Zeta,geta);
    end
    if isSphere
        if n == 0
            du_rdt(:,indices) = repmat(G./K.*(4*a.^2-3*b.^2).*A(:,1).*sOrigin_xi/9,1,sum(indices));
        elseif n == 2
            du_rdt(:,indices) = repmat(-(a.^2.*A(:,1).*sOrigin_xi - 3*b.^2.*B(:,1).*sOrigin_eta)/15,1,sum(indices));
        else
            du_rdt(:,indices) = 0;
        end
    else
        du_rdt = du_rdt + A(:,2)*Q1sr.*wyxi.*S_(1,2,n,xi,eta,Zxi,gxi);
        if n > 0
        	du_rdt = du_rdt + B(:,2)*Q1sr.*wyeta.*T_(1,2,n,eta,Zeta,geta);
        end
    end
    solid.du_rdt = du_rdt;
end

if layer{m}.calc_du_tdr
    Q1sr2 = Q1s./r2;
    du_tdr = A(:,1)*Q1sr2.*wjxi.*S_(4,1,n,xi,eta,Zxi,gxi);
    if n > 0
        du_tdr = du_tdr + B(:,1)*Q1sr2.*wjeta.*T_(4,1,n,eta,Zeta,geta);
    end
    if isSphere
        if n == 0
            du_tdr(:,indices) = repmat(G./K.*(4*a.^2-3*b.^2).*A(:,1).*sOrigin_xi/9,1,sum(indices));
        elseif n == 2
            du_tdr(:,indices) = repmat(2*(a.^2.*A(:,1).*sOrigin_xi - 3*b.^2.*B(:,1).*sOrigin_eta)/15,1,sum(indices));
        else
            du_tdr(:,indices) = 0;
        end
    else
        du_tdr = du_tdr + A(:,2)*Q1sr2.*wyxi.*S_(4,2,n,xi,eta,Zxi,gxi);
        if n > 0
        	du_tdr = du_tdr + B(:,2)*Q1sr2.*wyeta.*T_(4,2,n,eta,Zeta,geta);
        end
    end
    solid.du_tdr = du_tdr;
end

if layer{m}.calc_du_tdt
    Q2r1 = Q2./r;
    du_tdt = A(:,1)*Q2r1.*wjxi.*S_(2,1,n,xi,eta,Zxi,gxi);
    if n > 0
        du_tdt = du_tdt + B(:,1)*Q2r1.*wjeta.*T_(2,1,n,eta,Zeta,geta);
    end
    if isSphere
        du_tdt(:,logical(r < eps)) = 0;
    else
        du_tdt = du_tdt + A(:,2)*Q2r1.*wyxi.*S_(2,2,n,xi,eta,Zxi,gxi);
        if n > 0
        	du_tdt = du_tdt + B(:,2)*Q2r1.*wyeta.*T_(2,2,n,eta,Zeta,geta);
        end
    end
    solid.du_tdt = du_tdt;
end

if layer{m}.calc_sigma_rr
    Q0r2 = Q0./r2;
    sigma_rr = A(:,1)*Q0r2.*wjxi.*S_(5,1,n,xi,eta,Zxi,gxi);
    if n > 0
        sigma_rr = sigma_rr + B(:,1)*Q0r2.*wjeta.*T_(5,1,n,eta,Zeta,geta);
    end
    if isSphere
        if n == 0
            sigma_rr(:,indices) = repmat(5*(4*a.^2-3*b.^2).*A(:,1).*sOrigin_xi/30,1,sum(indices));
        elseif n == 2
            sigma_rr(:,indices) = repmat((-2*a.^2.*A(:,1).*sOrigin_xi + 6*b.^2.*B(:,1).*sOrigin_eta)/30,1,sum(indices));
        else
            sigma_rr(:,indices) = 0;
        end
    else
        sigma_rr = sigma_rr + A(:,2)*Q0r2.*wyxi.*S_(5,2,n,xi,eta,Zxi,gxi);
        if n > 0
        	sigma_rr = sigma_rr + B(:,2)*Q0r2.*wyeta.*T_(5,2,n,eta,Zeta,geta);
        end
    end
	sigma_rr = 2*G.*sigma_rr;
    
	solid.sigma_rr = sigma_rr;
end
if layer{m}.calc_sigma_tt
    Q0r2 = Q0./r2;
    Q2r2 = Q2./r2;
    sigma_tt =   A(:,1)*Q0r2.*wjxi.*S_(6,1,n,xi,eta,Zxi,gxi)  ...
               + A(:,1)*Q2r2.*wjxi.*S_(2,1,n,xi,eta,Zxi,gxi);
    if n > 0
        sigma_tt = sigma_tt + B(:,1)*Q0r2.*wjeta.*T_(6,1,n,eta,Zeta,geta) ...
                            + B(:,1)*Q2r2.*wjeta.*T_(2,1,n,eta,Zeta,geta);
    end
    if isSphere
        if n == 0
            sigma_tt(:,indices) = repmat(5*(4*a.^2-3*b.^2).*A(:,1).*sOrigin_xi/30,1,sum(indices));
        elseif n == 2
            sigma_tt(:,indices) = repmat((-2*a.^2.*A(:,1).*sOrigin_xi + 6*b.^2.*B(:,1).*sOrigin_eta)/30,1,sum(indices));
        else
            sigma_tt(:,indices) = 0;
        end
    else
        sigma_tt = sigma_tt + A(:,2)*Q0r2.*wyxi.*S_(6,2,n,xi,eta,Zxi,gxi) ...
                            + A(:,2)*Q2r2.*wyxi.*S_(2,2,n,xi,eta,Zxi,gxi);
        if n > 0
            sigma_tt = sigma_tt + B(:,2)*Q0r2.*wyeta.*T_(6,2,n,eta,Zeta,geta) ...
                                + B(:,2)*Q2r2.*wyeta.*T_(2,2,n,eta,Zeta,geta);
        end
    end
	sigma_tt = 2*G.*sigma_tt;
    
	solid.sigma_tt = sigma_tt;
end
if layer{m}.calc_sigma_pp
    Q0r2 = Q0./r2;
    cott_Q1r2 = cos(theta).*Q1s./r2;
    sigma_pp =   A(:,1)*Q0r2.*wjxi.*S_(6,1,n,xi,eta,Zxi,gxi) ...
               + A(:,1)*cott_Q1r2.*wjxi.*S_(2,1,n,xi,eta,Zxi,gxi);
    if n > 0
        sigma_pp = sigma_pp + B(:,1)*Q0r2.*wjeta.*T_(6,1,n,eta,Zeta,geta) ...
                            + B(:,1)*cott_Q1r2.*wjeta.*T_(2,1,n,eta,Zeta,geta);
    end
    if isSphere
        if n == 0
            sigma_pp(:,indices) = repmat(5*(4*a.^2-3*b.^2).*A(:,1).*sOrigin_xi/30,1,sum(indices));
        elseif n == 2
            sigma_pp(:,indices) = repmat((4*a.^2.*A(:,1).*sOrigin_xi - 12*b.^2.*B(:,1).*sOrigin_eta)/30,1,sum(indices));
        else
            sigma_pp(:,indices) = 0;
        end
    else
        sigma_pp = sigma_pp + A(:,2)*Q0r2.*wyxi.*S_(6,2,n,xi,eta,Zxi,gxi) ...
                            + A(:,2)*cott_Q1r2.*wyxi.*S_(2,2,n,xi,eta,Zxi,gxi);
        if n > 0
            sigma_pp = sigma_pp + B(:,2)*Q0r2.*wyeta.*T_(6,2,n,eta,Zeta,geta) ...
                                + B(:,2)*cott_Q1r2.*wyeta.*T_(2,2,n,eta,Zeta,geta);
        end
    end
	sigma_pp = 2*G.*sigma_pp;
    
	solid.sigma_pp = sigma_pp;
end
if layer{m}.calc_sigma_rt
    Q1r2 = Q1./r2;
    sigma_rt = A(:,1)*Q1r2.*wjxi.*S_(7,1,n,xi,eta,Zxi,gxi);
    if n > 0
        sigma_rt = sigma_rt + B(:,1)*Q1r2.*wjeta.*T_(7,1,n,eta,Zeta,geta);
    end
    if isSphere
        sigma_rt(:,logical(r < eps)) = 0;
    else
        sigma_rt = sigma_rt + A(:,2)*Q1r2.*wyxi.*S_(7,2,n,xi,eta,Zxi,gxi);
        if n > 0
            sigma_rt = sigma_rt + B(:,2)*Q1r2.*wyeta.*T_(7,2,n,eta,Zeta,geta);
        end
    end
	sigma_rt = 2*G.*sigma_rt;
    
	solid.sigma_rt = sigma_rt;
end
if layer{m}.calc_navier1 || layer{m}.calc_navier2
    Q1sr2 = Q1s./r2;
    sigma_rt = A(:,1)*Q1sr2.*wjxi.*S_(7,1,n,xi,eta,Zxi,gxi);
    if n > 0
        sigma_rt = sigma_rt + B(:,1)*Q1sr2.*wjeta.*T_(7,1,n,eta,Zeta,geta);
    end
    if ~isSphere
        sigma_rt = sigma_rt + A(:,2)*Q1sr2.*wyxi.*S_(7,2,n,xi,eta,Zxi,gxi);
        if n > 0
        	sigma_rt = sigma_rt + B(:,2)*Q1sr2.*wyeta.*T_(7,2,n,eta,Zeta,geta);
        end
    end
	sigma_rt = 2*G.*sigma_rt;
    
    r3 = r.^3;
    Q0r3 = Q0./r3;
    Q1r3 = Q1./r3;
    Q2r3 = Q2./r3;
    dsigma_rr_dr = A(:,1)*Q0r3.*wjxi.*S_(8,1,n,xi,eta,Zxi,gxi);
    if n > 0
    	dsigma_rr_dr = dsigma_rr_dr + B(:,1)*Q0r3.*wjeta.*T_(8,1,n,eta,Zeta,geta);
    end
    if ~isSphere
        dsigma_rr_dr =  dsigma_rr_dr + A(:,2)*Q0r3.*wyxi.*S_(8,2,n,xi,eta,Zxi,gxi);
        if n > 0
        	dsigma_rr_dr =  dsigma_rr_dr + B(:,2)*Q0r3.*wyeta.*T_(8,2,n,eta,Zeta,geta);
        end
    end
	dsigma_rr_dr = 2*G.*dsigma_rr_dr;
      
    dsigma_rt_dr = A(:,1)*Q1r3.*wjxi.*S_(9,1,n,xi,eta,Zxi,gxi);
    if n > 0
        dsigma_rt_dr = dsigma_rt_dr + B(:,1)*Q1r3.*wjeta.*T_(9,1,n,eta,Zeta,geta);
    end
    if ~isSphere
        dsigma_rt_dr =  dsigma_rt_dr + A(:,2)*Q1r3.*wyxi.*S_(9,2,n,xi,eta,Zxi,gxi);
        if n > 0
        	dsigma_rt_dr =  dsigma_rt_dr + B(:,2)*Q1r3.*wyeta.*T_(9,2,n,eta,Zeta,geta);
        end
    end
	dsigma_rt_dr = 2*G.*dsigma_rt_dr;

    dsigma_rt_dt = A(:,1)*Q2r3.*wjxi.*S_(7,1,n,xi,eta,Zxi,gxi);
    if n > 0
        dsigma_rt_dt = dsigma_rt_dt + B(:,1)*Q2r3.*wjeta.*T_(7,1,n,eta,Zeta,geta);
    end
    if ~isSphere
        dsigma_rt_dt = dsigma_rt_dt + A(:,2)*Q2r3.*wyxi.*S_(7,2,n,xi,eta,Zxi,gxi);
        if n > 0
        	dsigma_rt_dt = dsigma_rt_dt + B(:,2)*Q2r3.*wyeta.*T_(7,2,n,eta,Zeta,geta);
        end
    end
	dsigma_rt_dt = 2*G.*dsigma_rt_dt;

    dsigma_diffr =   A(:,1)*Q1r3.*wjxi.*S_(6,1,n,xi,eta,Zxi,gxi) ...
                   + (-n^2-n+1)*A(:,1)*Q1r3.*wjxi.*S_(2,1,n,xi,eta,Zxi,gxi);
    if n > 0
        dsigma_diffr = dsigma_diffr +  B(:,1)*Q1r3.*wjeta.*T_(6,1,n,eta,Zeta,geta) ...
                                    + (-n^2-n+1)*B(:,1)*Q1r3.*wjeta.*T_(2,1,n,eta,Zeta,geta);
    end
    if ~isSphere
        dsigma_diffr = dsigma_diffr + A(:,2)*Q1r3.*wyxi.*S_(6,2,n,xi,eta,Zxi,gxi) ...
                                    + (-n^2-n+1)*A(:,2)*Q1r3.*wyxi.*S_(2,2,n,xi,eta,Zxi,gxi);
        if n > 0
            dsigma_diffr = dsigma_diffr + B(:,2)*Q1r3.*wyeta.*T_(6,2,n,eta,Zeta,geta) ...
                                        + (-n^2-n+1)*B(:,2)*Q1r3.*wyeta.*T_(2,2,n,eta,Zeta,geta);
        end
    end
	dsigma_diffr = 2*G.*dsigma_diffr;
    
    R = repmat(r,size(xi,1),1);
    Theta = repmat(theta,size(xi,1),1);
    solid.navier1 = dsigma_rr_dr + dsigma_rt_dt + 1./R.*(2*sigma_rr-sigma_tt-sigma_pp+sigma_rt.*cos(Theta));
    solid.navier2 = dsigma_rt_dr + dsigma_diffr + 3./R.*sigma_rt.*sin(Theta);
    if isSphere
        if n == 1
            rho = layer{m}.rho;
            solid.navier1(:,indices) = repmat(-omega.^2.*rho.*(a.*A(:,1).*sOrigin_xi - 2*b.*B(:,1).*sOrigin_eta)/3,1,sum(indices));
        else
            solid.navier1(:,indices) = 0;
        end
        solid.navier2(:,indices) = 0;
    end
end
