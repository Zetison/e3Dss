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
                'Display',          'final', ...     % Print options ('final', 'iter' or 'none')
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

fieldsToRemove = {'r','theta','phi','a_temp','a','b','b_temp','G_temp','K_temp','Z_r_s','K','G','P','dP','d2P', ...
                  'Z_xi','Z_eta','Z_xi_i','Z_eta_i','Z_xi_o','Z_eta_o','Z_xi_i_s','Z_xi_o_s','Z_eta_i_s','Z_eta_o_s','Z_r_s', ...
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
    calc_err_helmholtz = layer{m}.calc_err_helmholtz;
    calc_dp = or(layer{m}.calc_dp, calc_err_dc);

    calc_p_laplace = layer{m}.calc_p_laplace;

    layer{m}.calc_p = layer{m}.calc_p;
    layer{m}.calc_dpdr = any(calc_dp);
    layer{m}.calc_dp = calc_dp;
    layer{m}.calc_p_laplace = calc_p_laplace || calc_err_helmholtz;
    layer{m}.calc_dpdt = layer{m}.calc_dpdr;
    layer{m}.calc_d2pdr2 = calc_p_laplace || calc_err_helmholtz;
    layer{m}.calc_d2pdt2 = calc_p_laplace || calc_err_helmholtz;
    layer{m}.calc_errors = calc_err_dc || calc_err_pc || calc_err_helmholtz;
    layer{m}.calc_farFieldOnly = ~any([layer{m}.calc_p, layer{m}.calc_dpdr, layer{m}.calc_dpdt, layer{m}.calc_d2pdr2, layer{m}.calc_d2pdt2]);

    layer{m}.calc_p_inc = (layer{m}.calc_p_inc || layer{m}.calc_errors) && m == m_s;
    layer{m}.calc_dp_inc = and(or(layer{m}.calc_dp_inc, layer{m}.calc_errors), m == m_s);
    layer{m}.calc_dp_incdr = and(any(layer{m}.calc_dp_inc), options.p_inc_fromSeries);
    layer{m}.calc_dp_incdt = and(any(layer{m}.calc_dp_inc), options.p_inc_fromSeries);

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
if noEvaluationPts
    error('No evaluation point (layer{m}.X) was provided')
end
options.SHBC = strcmp(options.BC,'SHBC');
options.SSBC = strcmp(options.BC,'SSBC');
options.IBC = strcmp(options.BC,'IBC');
if options.IBC && ~isfield(options,'z')
    options.z = ones(1,prec);
end

if any(omega == 0) && ( strcmp(options.applyLoad, 'mechExcitation') || ...
                       (strcmp(options.applyLoad, 'pointCharge') && options.r_s ~= 0) || ...
                       (strcmp(options.applyLoad, 'surfExcitation') && abs(options.theta_s - PI) > options.Eps))
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
        theta = repmat(layer{m}.theta,nFreqs,1);
        phi = repmat(layer{m}.phi,nFreqs,1);
        r = repmat(layer{m}.r,nFreqs,1);
        
        indices = logical(r < options.Eps);
        n_X = size(layer{m}.X,1);

        if any(layer{m}.calc_dp) || layer{m}.calc_p_laplace
            dpdr = layer{m}.dpdr;
            dpdt = sin(theta).*layer{m}.dpdt; % rescale dpdt
            dpdX_m = cell(3,1);

            dpdX_m{1} = dpdr.*sin(theta).*cos(phi) + dpdt.*cos(theta).*cos(phi)./r;
            dpdX_m{1}(indices) = 0;

            dpdX_m{2} = dpdr.*sin(theta).*sin(phi) + dpdt.*cos(theta).*sin(phi)./r;
            dpdX_m{2}(indices) = 0;

            dpdX_m{3} = dpdr.*cos(theta) - dpdt.*sin(theta)./r;
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

                temp = 2./r.*dpdr + 1./r.^2.*cos(theta).*dpdt + 1./r.^2.*d2pdt2;                    
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
        switch options.applyLoad
            case {'planeWave','pointCharge'}
                if any(layer{m}.calc_dp_inc) && m == m_s
                    if options.p_inc_fromSeries
                        dp_incdr = layer{m}.dp_incdr;
                        dp_incdt = sin(theta).*layer{m}.dp_incdt; % rescale dp_incdt
                        dp_incdX_m = cell(3,1);

                        dp_incdX_m{1} = dp_incdr.*sin(theta).*cos(phi) + dp_incdt.*cos(theta).*cos(phi)./r;
                        dp_incdX_m{1}(indices) = 0;

                        dp_incdX_m{2} = dp_incdr.*sin(theta).*sin(phi) + dp_incdt.*cos(theta).*sin(phi)./r;
                        dp_incdX_m{2}(indices) = 0;

                        dp_incdX_m{3} = dp_incdr.*cos(theta) - dp_incdt.*sin(theta)./r;
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
                    else
                        a = layer{1}.a;
                        if strcmp(options.applyLoad,'planeWave')
                            a_vec = a*options.d_vec.';
                            layer{1}.p_inc = P_inc.*exp(1i*a_vec*layer{1}.X.');
                            for ii = 1:3
                                if layer{m}.calc_dp_inc(ii)
                                    layer{m}.dp_inc{ii} = 1i*a_vec(:,ii).*layer{1}.p_inc;
                                end
                            end
                        else
                            x_s = r_s*options.d_vec.';
                            Xxms = layer{m}.X-x_s;
                            nXxms = norm2(Xxms).';
                            p_inc = P_inc.*exp(1i*a.*nXxms)./nXxms;
                            if layer{m}.calc_p_inc
                                layer{m}.p_inc = p_inc;
                            end
                            for ii = 1:3
                                if layer{m}.calc_dp_inc(ii)
                                    layer{m}.dp_inc{ii} = p_inc.*(1i*a - 1./nXxms).*Xxms(:,ii).'./nXxms;
                                end
                            end
                        end
                    end
                end
            case 'radialPulsation' 
                if layer{m}.calc_p_inc
                    a = layer{1}.a;
                    layer{m}.p_inc = P_inc.*exp(1i*a.*(r-layer{1}.R))./r;
                end
            otherwise
                if any(layer{m}.calc_dp_inc) || layer{m}.calc_p_inc
                    error('Not valid')
                end
        end
        if any(layer{m}.calc_u)
            u_X_m = cell(3,1);
            u_r = layer{m}.u_r;
            u_t = layer{m}.u_t.*sin(theta); % rescale u_t
            u_X_m{1} = u_r.*sin(theta).*cos(phi) + u_t.*cos(theta).*cos(phi);
            u_X_m{2} = u_r.*sin(theta).*sin(phi) + u_t.*cos(theta).*sin(phi);
            u_X_m{3} = u_r.*cos(theta) - u_t.*sin(theta);
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
            du_m{1,2} = layer{m}.du_rdt.*sin(theta); % rescale du_rdt
            du_m{1,3} = zeros(nFreqs,n_X,prec);
            du_m{2,1} = layer{m}.du_tdr.*sin(theta); % rescale du_tdr
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
        if layer{m}.calc_err_helmholtz
            a = layer{m}.a;
            aa = repmat(a,1,n_X);

            layer{m}.err_helmholtz = max(abs(layer{m}.Phi_laplace+aa.^2.*layer{m}.Phi),[],2)./max(abs(aa.^2.*layer{m}.Phi),[],2);
        end
        if any(layer{m}.calc_err_navier)
            theta = repmat(layer{m}.theta,nFreqs,1);
            rho = layer{m}.rho;
            u_r = layer{m}.u_r;
            u_t = layer{m}.u_t.*sin(theta); % rescale u_t
            Omega = repmat(omega,1,n_X);
            layer{m}.err_navier = cell(1,2);
            if layer{m}.calc_err_navier(1)
                layer{m}.err_navier{1} = max(abs(layer{m}.navier1 + rho*Omega.^2.*u_r),[],2)./max(abs(rho*Omega.^2.*u_r),[],2);
            end
            if layer{m}.calc_err_navier(2)
                layer{m}.err_navier{2} = max(abs(layer{m}.navier2 + rho*Omega.^2.*u_t),[],2)./max(abs(rho*Omega.^2.*u_t),[],2);
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
            if m == m_s
                a = layer{m+1}.a;
                a_vec = options.d_vec*a;
                p_inc = P_inc*exp(1i*X_s*a_vec);
                sigma_rr = sigma_rr + p_inc;
            end
            switch nextMedia
                case 'void'
                    if SSBC
                        layer{m}.err_pc = max(abs(sigma_rr),[],2)/P_inc;
                    end
                case {'solid','fluid'}
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
                    if m2 == m_s
                        a = layer{m+1}.a;
                        a_vec = options.d_vec*a;
                        p_inc = P_inc*exp(1i*X_s*a_vec);
                        sigma_rr2 = sigma_rr2 + p_inc;
                    end
                    layer{m}.err_pc = max(abs(sigma_rr-sigma_rr2),[],2)./max(abs(sigma_rr),[],2);
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
                case 'void'
                    if SHBC
                        layer{m}.err_dc = max(abs(u{1}.*n_x + u{2}.*n_y + u{3}.*n_z),[],2)./P_inc;
                    end
                case {'solid','fluid'}
                    u = cell(1,3);
                    for i = 1:3
                        u{i} = layer{m}.u{i}(:,indices);
                    end
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
    layer{m}.R               = 1;          % Inner radius of layer
    layer{m}.rho             = 1000;       % Mass density
    switch layer{m}.media
        case 'fluid'
            layer{m}.c  	= 1500;       % Speed of sound
            layer{m}.mu 	= 0;          % Shear coefficient of viscosity
            layer{m}.mu_b	= 0;          % Bulk coefficient of viscosity
        case 'solid'
            layer{m}.E   	= 200e9;      % Youngs modulus for solid layers
            layer{m}.nu  	= 0.3;        % Poisson ratio for solid layers
    end
    layer{m}.lossFactor    	 = [0,0];      % Hysteretic loss factor (values around 0.001 for lightly damped materials, values around 0.01 for moderately damped materials and values around 0.1 for heavily damped materials, first component for longitudinal velocities and second component for transverse velocities)
    layer{m}.sqrtHysterisis  = false;      % If true, the square root is used when including the lossfactor (i.e. c -> c*sqrt(1-1i*lossFactor)) as opposed to the linear version (i.e. c -> c*(1-1i*lossFactor)) 
    
    layer{m}.calc_err_dc     = false;      % Calculate the errors for the displacement conditions
    layer{m}.calc_err_pc  	 = false;      % Calculate the errors for the pressure conditions    
    layer{m}.calc_p_0    	 = false;      % Toggle calculation of the far field pattern
    layer{m}.calc_p       	 = false;      % Toggle calculation of the scattered pressure
    layer{m}.calc_dp     	 = false(1,3); % Toggle calculation of the three components of the gradient of the pressure
    layer{m}.calc_p_laplace	 = false;      % Toggle calculation of the Laplace operator of the scattered pressure fields
    layer{m}.calc_p_inc   	 = false;      % Toggle calculation of the incident pressure
    layer{m}.calc_dp_inc   	 = false(1,3); % Toggle calculation of the three components of the gradient of the incident pressure
    layer{m}.calc_u       	 = false(1,3); % Toggle calculation of the three components of the displacement
    layer{m}.calc_du       	 = false(3,3); % Toggle calculation of the three cartesian derivatives of the three components of the displacement [du_xdx du_xdy du_xdz; 
                                           %                                                                                                    du_ydx du_ydy du_ydz; 
                                           %                                                                                                    du_zdx du_zdy du_zdz]
    layer{m}.calc_sigma  	 = false(1,6); % Toggle calculation of the six components of the stress field (cartesian coordinates) [sigma_xx sigma_yy sigma_zz sigma_yz sigma_xz sigma_xy]
    layer{m}.calc_sigma_s 	 = false(1,6); % Toggle calculation of the six components of the stress field (spherical coordinates) [sigma_rr sigma_tt sigma_pp sigma_tp sigma_rp sigma_rt]
    layer{m}.calc_err_navier = false(1,2); % Toggle calculation of the errors for the Navier equation
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [layer,N_eps,flag,relTermMaxArr] = e3Dss_0(layer, options)
% This function computes exact 3D scattering solutions when the axis of
% symmetry is the z-axis
%
%
% Note that dpdt, u_t, du_tdr, and du_rdt are scaled by csc(theta)
fluidFieldNames = {};
fieldNames = {'p_0','p','dpdr','dpdt','d2pdr2','d2pdt2','p_inc','dp_incdr','dp_incdt','u_r','u_t','du_rdr','du_rdt','du_tdr','du_tdt','sigma_rr','sigma_tt','sigma_pp','sigma_rt','navier1','navier2'};

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
    rho = layer{m}.rho;
    switch layer{m}.media
        case 'fluid'
            mu = layer{m}.mu;
            mu_b = layer{m}.mu_b;
            c = layer{m}.c;
            
            % Compute derived quantities            
            K = rho.*c.^2 - 1i*omega.*mu_b;
            G = -1i*omega.*mu;
            
        case 'solid'
            E = layer{m}.E;
            nu = layer{m}.nu;
            
            % Compute derived quantities            
            K = E./(3*(1-2*nu));
            G = E./(2*(1+nu));
    end
    c_s_1 = sqrt((3*K+4*G)./(3*rho)); % longitudinal wave velocity 
    c_s_2 = sqrt(G./rho); % shear wave velocity

    % Add hysteris damping
    if layer{m}.sqrtHysterisis
        c_s_1 = c_s_1.*sqrt(1-1i*layer{m}.lossFactor(:,1));
        c_s_2 = c_s_2.*sqrt(1-1i*layer{m}.lossFactor(:,2));
    else
        c_s_1 = c_s_1.*(1-1i*layer{m}.lossFactor(:,1));
        c_s_2 = c_s_2.*(1-1i*layer{m}.lossFactor(:,2));
    end

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
                if strcmp(options.applyLoad,'pointCharge')
                    layer{m}.Z_r_s{i,j} = zeros(nFreqs,1,prec);
                end
            end
        end
    end
    for fieldName = fieldNames
        if layer{m}.(['calc_' fieldName{1}]) && ~(strcmp(fieldName{1},'p_inc') && ~options.p_inc_fromSeries)
            layer{m}.(fieldName{1}) = zeros(nFreqs,n_X,prec);
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
            case 'solid'
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
singleModeSolution = (strcmp(options.applyLoad,'pointCharge') && options.r_s == 0) || strcmp(options.applyLoad,'radialPulsation') ...
                     || (strcmp(options.applyLoad,'surfExcitation') && options.theta_s(1) == 0 && abs(options.theta_s(2) - pi) < Eps);

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
            case 'solid'
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
            media = getTerm_n(layer,m,n,xi,eta,A,B,Rt_m,omega_temp,isSphere,isOuterDomain,applyLoad,besselIndices,nu_a,exponentShift,options);
            
            % Update fields
            for fieldName = fieldNames
                if layer{m}.(['calc_' fieldName{1}]) && ~(strcmp(fieldName{1},'p_inc') && ~options.p_inc_fromSeries)
                    [layer,hasCnvrgdTmp2,hasDvrgdTmp2,relTermMax] = updateSum(layer,m,media,fieldName{1},indices,hasCnvrgdTmp2,hasDvrgdTmp2,relTermMax,tiny,Eps,prec,n);
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
function Z = iterate_Z(n,x,Z,Zindices,besselIndices,nu_a,U_pol,u_k,v_k,Eps)
if n == 0
    for i = 1:numel(besselIndices)
        if besselIndices(i)
            Z{i,2} = bessel_s(n,x,i,nu_a,U_pol,u_k,v_k,Eps);
        end
    end
end
for i = 1:numel(besselIndices)
    if besselIndices(i)
        Z{i,1} = Z{i,2}(Zindices,:);
        Z{i,2} = bessel_s(n+1,x,i,nu_a,U_pol,u_k,v_k,Eps);
    end
end

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

function [layer,hasCnvrgd,hasDvrgd,relTermMax] = updateSum(layer,m,media,fieldName,indices,hasCnvrgd,hasDvrgd,relTermMax,tiny,Eps,prec,n)
% Track any Inf/NaN in next term
hasDvrgdSub = or(hasDvrgd,or(any(isinf(media.(fieldName)),2),any(isnan(media.(fieldName)),2)));

% Add non-Inf/NaN terms to series
layer{m}.(fieldName)(indices(~hasDvrgdSub),:) = layer{m}.(fieldName)(indices(~hasDvrgdSub),:) + media.(fieldName)(~hasDvrgdSub,:);

% Compute relative error
relTerm = abs(media.(fieldName))./(abs(layer{m}.(fieldName)(indices,:))+tiny);

% Update converged/diverged tracking arrays
hasCnvrgd = hasCnvrgd.*prod(relTerm < Eps, 2);
hasDvrgd = hasDvrgd + hasDvrgdSub;

% Only track frequencies that have not converged
nonConvergedIdx = false(size(relTermMax,1),1);
nonConvergedIdx(indices) = true;
maxRelTerm = -Inf(size(relTermMax,1),1,prec);
maxRelTerm(nonConvergedIdx) = max(relTerm,[],2);
indices2 = maxRelTerm > relTermMax;
relTermMax(indices2) = maxRelTerm(indices2);

function Rt_m = getIntermediateRadius(layer,m,isSphere,isOuterDomain)
if isOuterDomain
    Rt_m = layer{1}.R;
elseif isSphere
    Rt_m = layer{m-1}.R;
else
    Rt_m = (layer{m}.R+layer{m-1}.R)/2;
end
