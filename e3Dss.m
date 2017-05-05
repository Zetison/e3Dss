function data = e3Dss(vv, newOptions)

%% Dirichlet and Robin conditions has not been implemented.
options = struct('d_vec',   [0;0;1], ...
                 'omega',   2*pi*1e3, ...
                 'R_i',     [], ...
                 'R_o',     1, ...
                 'P_inc',   1, ...
                 'E',       [], ...
                 'nu',      [], ...
                 'rho_s',   [], ...
                 'rho_f',   [], ...
                 'c_f',     1500, ...
                 'Eps',     eps, ...
                 'usePlaneWave',        1, ...
                 'usePointChargeWave',  0, ...
                 'r_s',             NaN, ...
                 'calc_p',          1, ...
                 'calc_dpdx',       0, ...
                 'calc_dpdy',       0, ...
                 'calc_dpdz',       0, ...
                 'calc_u_x',        0, ...
                 'calc_u_y',        0, ...
                 'calc_u_z',        0, ...
                 'calc_u_r',        0, ...
                 'calc_u_t',        0, ...
                 'calc_sigma_xx',   0, ...
                 'calc_sigma_yy',   0, ...
                 'calc_sigma_zz',   0, ...
                 'calc_sigma_yz',   0, ...
                 'calc_sigma_xz',   0, ...
                 'calc_sigma_xy',   0, ...
                 'calc_sigma_rr',   0, ...
                 'calc_sigma_tt',   0, ...
                 'calc_sigma_pp',   0, ...
                 'calc_sigma_rt',   0, ...
                 'calc_navier',     0, ...
                 'calc_p_laplace',  0, ...
                 'calc_farField',   0, ...
                 'calc_errorsDisplacementCondition', 0, ...
                 'calc_errorsPressureCondition', 	 0, ...
                 'calc_errorsHelmholtz',             0, ...
                 'calc_errorsNavier',                0, ...
                 'useSymbolicPrecision',             0, ...
                 'N_max', inf);
if nargin > 1
    newOptionFields = fieldnames(newOptions);
    for j = 1:numel(newOptionFields)
        options.(newOptionFields{j}) = newOptions.(newOptionFields{j});
    end
end

if ~iscell(vv)
    vv = {vv};
end
useSymbolicPrecision = options.useSymbolicPrecision;
if useSymbolicPrecision
    options.d_vec = vpa(options.d_vec);
    options.omega = vpa(options.omega);
    options.R_i = vpa(options.R_i);
    options.R_o = vpa(options.R_o);
    options.P_inc = vpa(options.P_inc);
    options.E = vpa(options.E);
    options.nu = vpa(options.nu);
    options.rho_s = vpa(options.rho_s);
    options.rho_f = vpa(options.rho_f);
    options.c_f = vpa(options.c_f);
    options.Eps = vpa(options.Eps);
    for i = 1:length(vv)
        vv{i} = vpa(vv{i});
    end
end
calc_errorsPressureCondition = options.calc_errorsDisplacementCondition;
calc_errorsDisplacementCondition = options.calc_errorsPressureCondition;
calc_errorsHelmholtz = options.calc_errorsHelmholtz;
calc_errorsNavier = options.calc_errorsNavier;
omega = options.omega;
E = options.E;
rho_f = options.rho_f;
R_o = options.R_o;
R_i = options.R_i;
d_vec = options.d_vec;
calc_p = options.calc_p || calc_errorsPressureCondition;
calc_dpdx = options.calc_dpdx || calc_errorsDisplacementCondition;
calc_dpdy = options.calc_dpdy || calc_errorsDisplacementCondition;
calc_dpdz = options.calc_dpdz || calc_errorsDisplacementCondition;
calc_u_x = options.calc_u_x;
calc_u_y = options.calc_u_y;
calc_u_z = options.calc_u_z;
calc_u_r = options.calc_u_r || calc_errorsDisplacementCondition;
calc_u_t = options.calc_u_t || calc_errorsDisplacementCondition;
calc_sigma_xx = options.calc_sigma_xx;
calc_sigma_yy = options.calc_sigma_yy;
calc_sigma_zz = options.calc_sigma_zz;
calc_sigma_yz = options.calc_sigma_yz;
calc_sigma_xz = options.calc_sigma_xz;
calc_sigma_xy = options.calc_sigma_xy;
calc_p_laplace = options.calc_p_laplace;
calc_sigma_rr = options.calc_sigma_rr;
calc_sigma_tt = options.calc_sigma_tt;
calc_sigma_pp = options.calc_sigma_pp;
calc_sigma_rt = options.calc_sigma_rt;
calc_navier = options.calc_navier || calc_errorsNavier;
calc_errors = calc_errorsDisplacementCondition || calc_errorsPressureCondition || calc_errorsHelmholtz ||calc_errorsNavier;

options.calc_dpdr = calc_dpdx || calc_dpdy || calc_dpdz || calc_p_laplace;
options.calc_dpdt = options.calc_dpdr;
options.calc_d2pdr2 = calc_p_laplace || calc_errorsHelmholtz;
options.calc_d2pdt2 = calc_p_laplace || calc_errorsHelmholtz;
options.calc_u_r = calc_u_r || calc_u_x || calc_u_y || calc_u_z || calc_navier;
options.calc_u_t = calc_u_t || calc_u_x || calc_u_y || calc_u_z || calc_navier;

calcStresses =  calc_sigma_xx || calc_sigma_yy || calc_sigma_zz || calc_sigma_yz || calc_sigma_xz ...
                             || calc_sigma_xy || calc_navier || calc_errorsPressureCondition;
options.calc_sigma_rr = calcStresses;
options.calc_sigma_tt = calcStresses;
options.calc_sigma_pp = calcStresses;
options.calc_sigma_rt = calcStresses;
options.calc_navier = calc_navier;

if isrow(omega)
    options.omega = omega';
end
if options.usePointChargeWave 
    if isnan(options.r_s)
        options.r_s = 2*R_o(1);
    elseif abs(options.r_s) < R_o(1)
        error('d must satisfy |d| > R_o(1)')
    end
    warning(['It is not implemented an efficient routine to evaluate Equation (D.7). ' ...
             'The built in MATLAB routine is used for this purpose.'])
end
if options.usePointChargeWave
    options.usePlaneWave = false;
end
if ~(options.usePointChargeWave || options.usePlaneWave)
    error('Only plane waves and point charge waves are implemented for the incident wave.')
end
if length(E) < length(R_o)
    SHBC = true;
    ESBC = false;
    SSBC = false;  
else
    SHBC = false;
    if length(R_o) > length(R_i) % Inner domain is a solid sphere
        ESBC = true;
        SSBC = false;  
    else
        ESBC = false;
        if length(R_o) == length(rho_f)
            SSBC = true;
        else
            SSBC = false;        
        end
    end
end
options.ESBC = ESBC;
options.SHBC = SHBC;
options.SSBC = SSBC;

if omega(1) == 0 && options.usePointChargeWave
    error('This case has a non unique solution')
end
    

%% Coordinate transformation
% Due to symmetry, we can do a coordinate transformation such that we only
% need to compute the solution for the special case k_vec = k*[0, 0, 1].
r = cell(size(vv));
theta = cell(size(vv));
phi = cell(size(vv));
for j = 1:length(vv)
    if ~isempty(vv{j})
        [r{j}, theta{j}, phi{j}, A] = coordTransform(vv{j}, d_vec);
        r{j} = r{j}.';
        theta{j} = theta{j}.';
        phi{j} = phi{j}.';    
    end
end

%% Compute the solution with d_vec = [0, 0, 1]
data_0 = e3Dss_0(r, theta, options);

%% Coordinate transformation (Transform back to original Cartesian coordinate system)

% Allocate memory
nFreqs = length(omega);

M = length(options.R_o);
if SHBC || ESBC || SSBC
    data(M).p = [];
else
    data(M+1).p = [];
end
m = 1;
for j = 1:length(vv)
    if ~isempty(r{j})
        if mod(j,2)
            if calc_p
                data(m).p = zeros(size(vv{j},1),nFreqs);
                if useSymbolicPrecision
                    data(m).p = vpa(data(m).p);
                end
            end
            if calc_dpdx
                data(m).dpdx = zeros(size(vv{j},1),nFreqs);
                if useSymbolicPrecision
                    data(m).dpdx = vpa(data(m).dpdx);
                end
            end
            if calc_dpdy
                data(m).dpdy = zeros(size(vv{j},1),nFreqs);
                if useSymbolicPrecision
                    data(m).dpdy = vpa(data(m).dpdy);
                end
            end
            if calc_dpdz
                data(m).dpdz = zeros(size(vv{j},1),nFreqs);
                if useSymbolicPrecision
                    data(m).dpdz = vpa(data(m).dpdz);
                end
            end
            if calc_p_laplace
                data(m).p_laplace = zeros(size(vv{j},1),nFreqs);
                if useSymbolicPrecision
                    data(m).p_laplace = vpa(data(m).p_laplace);
                end
            end
        else
            if calc_u_x
                data(m).u_x = zeros(size(vv{j},1),nFreqs);
                if useSymbolicPrecision
                    data(m).u_x = vpa(data(m).u_x);
                end
            end
            if calc_u_y
                data(m).u_y = zeros(size(vv{j},1),nFreqs);
                if useSymbolicPrecision
                    data(m).u_y = vpa(data(m).u_x);
                end
            end
            if calc_u_z
                data(m).u_z = zeros(size(vv{j},1),nFreqs);
                if useSymbolicPrecision
                    data(m).u_z = vpa(data(m).u_x);
                end
            end
        end
    end
    if ~mod(j,2)
        m = m + 1;
    end
end

m = 1;
for j = 1:length(vv)
    if calc_errors
        if mod(j,2)
            data(m).v_fluid = vv{j};
        else
            data(m).v_solid = vv{j};
        end
    end
    if ~isempty(r{j})
        Theta = repmat(theta{j},nFreqs,1);
        Phi = repmat(phi{j},nFreqs,1);
        R = repmat(r{j},nFreqs,1);
        indices = logical(R < eps);
        if mod(j,2)
            if calc_dpdx || calc_dpdy || calc_dpdz || calc_p_laplace
                dpdr = data_0(m).dpdr;
                dpdt = sin(Theta).*data_0(m).dpdt; % rescale dpdt
                dpdX_m = cell(3,1);

                dpdX_m{1} = dpdr.*sin(Theta).*cos(Phi) + dpdt.*cos(Theta).*cos(Phi)./R;
                dpdX_m{1}(indices) = 0;
                
                dpdX_m{2} = dpdr.*sin(Theta).*sin(Phi) + dpdt.*cos(Theta).*sin(Phi)./R;
                dpdX_m{2}(indices) = 0;

                dpdX_m{3} = dpdr.*cos(Theta) - dpdt.*sin(Theta)./R;
                dpdX_m{3}(indices) = dpdr(indices);

                dpdX = cell(3,1);
                dpdX{1} = zeros(size(vv{j},1),nFreqs);
                dpdX{2} = zeros(size(vv{j},1),nFreqs);
                dpdX{3} = zeros(size(vv{j},1),nFreqs);
                for ii = 1:3
                    for jj = 1:3
                        dpdX{ii} = dpdX{ii} + A(ii,jj)*dpdX_m{jj}.';
                    end
                end
                if calc_dpdx
                    data(m).dpdx = dpdX{1};
                end
                if calc_dpdy
                    data(m).dpdy = dpdX{2};
                end
                if calc_dpdz
                    data(m).dpdz = dpdX{3};
                end
                dpdt = data_0(m).dpdt; % use scaled version of dpdt for the laplace operator
                if calc_p_laplace
                    d2pdr2 = data_0(m).d2pdr2;
                    d2pdt2 = data_0(m).d2pdt2;

                    temp = 2./R.*dpdr + 1./R.^2.*cos(Theta).*dpdt + 1./R.^2.*d2pdt2;                    
                    temp(indices) = 0;
                    data(m).p_laplace = d2pdr2 + temp;
                    data(m).p_laplace = data(m).p_laplace.';
                end
            end
            if calc_p
                data(m).p = data_0(m).p.';
            end
        else
            if calcStresses
                % Transform the stresses in the spherical coordinate system to the 
                % Cartesian coordinate system
                sigma_m = cell(6,1);
                sigma_m{1} = data_0(m).sigma_rr;
                sigma_m{2} = data_0(m).sigma_tt;
                sigma_m{3} = data_0(m).sigma_pp;
                sigma_m{4} = zeros(nFreqs,size(vv{j},1));
                sigma_m{5} = zeros(nFreqs,size(vv{j},1));
                sigma_m{6} = data_0(m).sigma_rt;
                if useSymbolicPrecision
                    sigma_m{4} = vpa(sigma_m{4});
                    sigma_m{5} = vpa(sigma_m{5});
                end
                D = getStressTransformationMatrix(theta{j},phi{j},2);
                sigma_X_m = cell(6,1);
                for ii = 1:6
                    sigma_X_m{ii} = zeros(nFreqs,size(vv{j},1));
                    if useSymbolicPrecision
                        sigma_X_m{ii} = vpa(sigma_X_m{ii});
                    end
                end
                for ii = 1:6
                    for jj = 1:6
                        D_kl = D(ii, jj, :);
                        D_kl = D_kl(:);
                        sigma_X_m{ii} = sigma_X_m{ii} + repmat(D_kl.',nFreqs,1).*sigma_m{jj};
                    end
                end
                if ESBC && m == M
                    sigma_X_m{1}(indices) = sigma_m{1}(indices);
                    sigma_X_m{2}(indices) = sigma_m{2}(indices);
                    sigma_X_m{3}(indices) = sigma_m{3}(indices);
                    sigma_X_m{4}(indices) = 0;
                    sigma_X_m{5}(indices) = 0;
                    sigma_X_m{6}(indices) = 0;
                end

                alpha = zeros(3);
                if useSymbolicPrecision
                    alpha = vpa(alpha);
                end
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
                sigma_X = cell(6,1);
                for ii = 1:6
                    sigma_X{ii} = zeros(size(vv{j},1),nFreqs);
                    if useSymbolicPrecision
                        sigma_X{ii} = vpa(sigma_X{ii});
                    end
                end
                for vgtIdx = 1:6
                    for ii = 1:3
                        for jj = 1:3
                            sigma_X{vgtIdx} = sigma_X{vgtIdx} + alpha(vgt(vgtIdx,1),ii)*alpha(vgt(vgtIdx,2),jj)*sigma_X_m{vgtinv(ii,jj)}.';
                        end
                    end
                end
                if calc_sigma_rr
                    data(m).sigma_rr = data_0(m).sigma_rr.';
                end
                if calc_sigma_tt
                    data(m).sigma_tt = data_0(m).sigma_tt.';
                end
                if calc_sigma_pp
                    data(m).sigma_pp = data_0(m).sigma_pp.';
                end
                if calc_sigma_rt
                    data(m).sigma_rt = data_0(m).sigma_rt.';
                end
                if calc_sigma_xx
                    data(m).sigma_xx = sigma_X{1};
                end
                if calc_sigma_yy
                    data(m).sigma_yy = sigma_X{2};
                end
                if calc_sigma_zz
                    data(m).sigma_zz = sigma_X{3};
                end
                if calc_sigma_yz
                    data(m).sigma_yz = sigma_X{4};
                end
                if calc_sigma_xz
                    data(m).sigma_xz = sigma_X{5};
                end
                if calc_sigma_xy
                    data(m).sigma_xy = sigma_X{6};
                end
            end
            if calc_u_x || calc_u_y || calc_u_z
                u_X_m = cell(3,1);
                u_r = data_0(m).u_r;
                u_t = data_0(m).u_t;
                u_X_m{1} = u_r.*sin(Theta).*cos(Phi) + u_t.*cos(Theta).*cos(Phi);
                u_X_m{2} = u_r.*sin(Theta).*sin(Phi) + u_t.*cos(Theta).*sin(Phi);
                u_X_m{3} = u_r.*cos(Theta) - u_t.*sin(Theta);
                if ESBC && m == M
                    u_X_m{1}(indices) = 0;
                    u_X_m{2}(indices) = 0;
                    u_X_m{3}(indices) = u_r(indices);
                end
                u_X = cell(3,1);
                u_X{1} = zeros(size(vv{j},1),nFreqs);
                u_X{2} = zeros(size(vv{j},1),nFreqs);
                u_X{3} = zeros(size(vv{j},1),nFreqs);
                if useSymbolicPrecision
                    u_X{1} = vpa(u_X{1});
                    u_X{2} = vpa(u_X{2});
                    u_X{3} = vpa(u_X{3});
                end
                for ii = 1:3
                    for jj = 1:3
                        u_X{ii} = u_X{ii} + A(ii,jj)*u_X_m{jj}.';
                    end
                end
                if calc_u_x
                    data(m).u_x = u_X{1};
                end
                if calc_u_y
                    data(m).u_y = u_X{2};
                end
                if calc_u_z
                    data(m).u_z = u_X{3};
                end
            end
            if calc_u_r || calc_navier
                data(m).u_r = data_0(m).u_r.';
            end
            if calc_u_t || calc_navier
                data(m).u_t = data_0(m).u_t.';
            end
            if calc_navier
                data(m).navier1 = data_0(m).navier1.';
                data(m).navier2 = data_0(m).navier2.';
            end
        end
    end
    if ~mod(j,2)
        m = m + 1;
    end
end
if calc_errors
    for m = 1:M+1
        if m ~= M+1
            for investigate = {'innerSurface','outerSurface'}
                if m ~= M+1 && ~((SHBC || ESBC) && m == M && strcmp(investigate{1},'innerSurface'))
                    iS = strcmp(investigate{1},'innerSurface');
                    Eps = options.Eps;
                    if SSBC && m == M && strcmp(investigate{1},'innerSurface')
                        v_f = vv{2*m};
                    else
                        v_f = vv{2*(m+iS)-1};
                        c_f = options.c_f(m+iS);
                    end
                    if m == M && SHBC
                        indices_f = abs(norm2(vv{2*m-1}) - R_o(m)) < 10*Eps;
                    elseif SSBC && m == M && strcmp(investigate{1},'innerSurface')
                        indices_f = abs(norm2(vv{2*m}) - R_i(m)) < 10*Eps;    
                        indices_s = indices_f;
                    else
                        [indices_f, indices_s] = findMatchingPoints(v_f,vv{2*m},Eps);
                    end
                    v_f = v_f(indices_f,:);
                    P_inc = options.P_inc;


                    if calc_errorsDisplacementCondition && ~(SSBC && m == M && strcmp(investigate{1},'innerSurface'))
                        rho_f = options.rho_f(m+iS);
                        
                        n_x = v_f(:,1)./norm2(v_f);
                        n_y = v_f(:,2)./norm2(v_f);
                        n_z = v_f(:,3)./norm2(v_f);
                        n_x = repmat(n_x,1,length(omega));
                        n_y = repmat(n_y,1,length(omega));
                        n_z = repmat(n_z,1,length(omega));
                        
                        dpdx = data(m+iS).dpdx(indices_f,:);
                        dpdy = data(m+iS).dpdy(indices_f,:);
                        dpdz = data(m+iS).dpdz(indices_f,:);
                        if m == 1 && ~iS
                            k = omega/c_f(1);
                            k_vec = d_vec*k;
                            p_inc = P_inc*exp(1i*v_f*k_vec);
                            dpdx = dpdx + 1i*p_inc.*repmat(k_vec(1,:),size(p_inc,1),1);
                            dpdy = dpdy + 1i*p_inc.*repmat(k_vec(2,:),size(p_inc,1),1);
                            dpdz = dpdz + 1i*p_inc.*repmat(k_vec(3,:),size(p_inc,1),1);
                        end
                        if SHBC && m == M && strcmp(investigate{1},'outerSurface')
                            err_dc = max(abs(dpdx.*n_x + dpdy.*n_y + dpdz.*n_z),[],1)/P_inc;
                        else
                            u_x = data(m).u_x(indices_s,:);
                            u_y = data(m).u_y(indices_s,:);
                            u_z = data(m).u_z(indices_s,:);
                            Omega = repmat(omega,size(u_x,1),1);
                            err_dc = max(abs(     (dpdx-rho_f*Omega.^2.*u_x).*n_x ...
                                                + (dpdy-rho_f*Omega.^2.*u_y).*n_y ...
                                                + (dpdz-rho_f*Omega.^2.*u_z).*n_z),[],1)./max(abs(dpdx.*n_x + dpdy.*n_y + dpdz.*n_z),[],1);
                        end
                        if ~isfield(data(m),'err_dc') || isempty(data(m).err_dc)
                            data(m).err_dc = err_dc;
                        else
                            data(m).err_dc = max([data(m).err_dc; err_dc],[],1);
                        end
                    end

                    if calc_errorsPressureCondition
                        if ~(SHBC && m == M)
                            if ~(SSBC && m == M && strcmp(investigate{1},'innerSurface'))
                                p_tot = data(m+iS).p(indices_f,:);
                                if m == 1 && ~iS
                                    k = omega/c_f(1);
                                    k_vec = d_vec*k;
                                    p_inc = P_inc*exp(1i*v_f*k_vec);
                                    p_tot = p_tot + p_inc;
                                end
                            end
                            sigma_cart = cell(6,1);
                            sigma_cart{1} = data(m).sigma_xx(indices_s,:);
                            sigma_cart{2} = data(m).sigma_yy(indices_s,:);
                            sigma_cart{3} = data(m).sigma_zz(indices_s,:);
                            sigma_cart{4} = data(m).sigma_yz(indices_s,:);
                            sigma_cart{5} = data(m).sigma_xz(indices_s,:);
                            sigma_cart{6} = data(m).sigma_xy(indices_s,:);

                            sigma_rr = zeros(size(sigma_cart{1}));
                            if useSymbolicPrecision
                                sigma_rr = vpa(sigma_rr);
                            end
                            phi = atan2(v_f(:,2),v_f(:,1));
                            r = sqrt(v_f(:,1).^2+v_f(:,2).^2+v_f(:,3).^2);
                            theta = acos(v_f(:,3)./r);
                            D = getStressTransformationMatrix(theta,phi,1);
                            for l = 1:6
                                D_kl = D(1, l, :);
                                D_kl = repmat(D_kl(:),1,length(omega));
                                sigma_rr = sigma_rr + D_kl.*sigma_cart{l};
                            end
                            if SSBC && m == M && strcmp(investigate{1},'innerSurface')
                                err_pc = max(abs(sigma_rr),[],1)/P_inc;
                            else
                                err_pc = max(abs(p_tot+sigma_rr),[],1)./max(abs(p_tot),[],1);
                            end
                            if ~isfield(data(m), 'err_pc') || isempty(data(m).err_pc)
                                data(m).err_pc = err_pc;
                            else
                                data(m).err_pc = max([data(m).err_pc; err_pc],[],1);
                            end
                        end
                    end
                end
            end
            if calc_errorsNavier && ~(SHBC && m == M)
                rho_s = options.rho_s(m);
                Omega = repmat(omega,size(data(m).u_r,1),1);
                data(m).err_navier1 = max(abs(data(m).navier1 + rho_s*Omega.^2.*data(m).u_r),[],1)./max(abs(rho_s*Omega.^2.*data(m).u_r),[],1);
                data(m).err_navier2 = max(abs(data(m).navier2 + rho_s*Omega.^2.*data(m).u_t),[],1)./max(abs(rho_s*Omega.^2.*data(m).u_t),[],1);
            end
        end
        if calc_errorsHelmholtz && ~(m == M+1 && (SSBC || SHBC || ESBC))
            p_laplace = data(m).p_laplace;
            k = omega/options.c_f(m);
            K = repmat(k,size(data(m).p,1),1);

            data(m).err_helmholtz = max(abs(p_laplace+K.^2.*data(m).p),[],1)./max(abs(K.^2.*data(m).p),[],1);
        end
    end
end
data(1).flag = data_0(1).flag;
data(1).N_eps = data_0(1).N_eps;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = e3Dss_0(r, theta, options)
% This function computes exact 3D scattering solutions when the axis of
% symmetry is the z-axis
%
%
% Note that dpdt are scaled by csc(theta)

omega = options.omega;
c_f = options.c_f;
E = options.E;
R_o = options.R_o;

k = omega*(1./c_f);

SSBC = options.SSBC;
ESBC = options.ESBC;
SHBC = options.SHBC;
useSymbolicPrecision = options.useSymbolicPrecision;

Eps = options.Eps;
nFreqs = length(omega);
M = length(R_o);
N_max = options.N_max;

% Compute derived quantities
if ~(SHBC && M == 1)
    nu = options.nu;
    rho_s = options.rho_s;
    K = E./(3*(1-2*nu));
    G = E./(2*(1+nu));
    options.G = G;

    c_s_1 = sqrt((3*K+4*G)./(3*rho_s)); % longitudinal wave velocity 
    c_s_2 = sqrt(G./rho_s); % shear wave velocity
    
    a = omega*(1./c_s_1);
    b = omega*(1./c_s_2);
end


if any(omega == 0)
    computeForStaticCase = true;
else
    computeForStaticCase = false;
end
% Allocate memory
P = cell(length(r),1);
dP = cell(length(r),1);
d2P = cell(length(r),1);
if ~(SHBC && M == 1)
    Z(M).xi = cell(2,2); % Spherical Bessel functions evaluated at a*r
    Z(M).eta = cell(2,2); % Spherical Bessel functions evaluated at b*r
end
if SHBC || ESBC || SSBC
    data(M).p = [];
    Z(M).zeta = cell(2,2); % Spherical Bessel functions evaluated at k*r
else
    data(M+1).p = [];
    Z(M+1).zeta = cell(2,2); % Spherical Bessel functions evaluated at k*r
end

m = 1;
for j = 1:length(r)
    if ~isempty(r{j})
        P{j} = zeros(2,length(theta{j})); 
        dP{j} = zeros(2,length(theta{j})); 
        d2P{j} = zeros(2,length(theta{j}));
        if useSymbolicPrecision
            P{j} = vpa(P{j});
            dP{j} = vpa(dP{j});
            d2P{j} = vpa(d2P{j});
        end
        if mod(j,2)
            Z(m).zeta{1,1} = zeros(nFreqs-computeForStaticCase,length(theta{j}));
            Z(m).zeta{1,2} = zeros(nFreqs-computeForStaticCase,length(theta{j}));
            Z(m).zeta{2,1} = zeros(nFreqs-computeForStaticCase,length(theta{j}));
            Z(m).zeta{2,2} = zeros(nFreqs-computeForStaticCase,length(theta{j}));
            if useSymbolicPrecision
                Z(m).zeta{1,1} = vpa(Z(m).zeta{1,1});
                Z(m).zeta{1,2} = vpa(Z(m).zeta{1,2});
                Z(m).zeta{2,1} = vpa(Z(m).zeta{2,1});
                Z(m).zeta{2,2} = vpa(Z(m).zeta{2,2});
            end
            if options.calc_p
                data(m).p = zeros(nFreqs,length(r{j}));
                if useSymbolicPrecision
                    data(m).p = vpa(data(m).p);
                end
            end
            if options.calc_dpdr
                data(m).dpdr = zeros(nFreqs,length(r{j}));
                if useSymbolicPrecision
                    data(m).dpdr = vpa(data(m).dpdr);
                end
            end
            if options.calc_dpdt
                data(m).dpdt = zeros(nFreqs,length(r{j}));
                if useSymbolicPrecision
                    data(m).dpdt = vpa(data(m).dpdt);
                end
            end
            if options.calc_d2pdr2
                data(m).d2pdr2 = zeros(nFreqs,length(r{j}));
                if useSymbolicPrecision
                    data(m).d2pdr2 = vpa(data(m).d2pdr2);
                end
            end
            if options.calc_d2pdt2
                data(m).d2pdt2 = zeros(nFreqs,length(r{j}));
                if useSymbolicPrecision
                    data(m).d2pdt2 = vpa(data(m).d2pdt2);
                end
            end
        else
            Z(m).xi{1,1} = zeros(nFreqs-computeForStaticCase,length(theta{j}));
            Z(m).xi{1,2} = zeros(nFreqs-computeForStaticCase,length(theta{j}));
            Z(m).xi{2,1} = zeros(nFreqs-computeForStaticCase,length(theta{j}));
            Z(m).xi{2,2} = zeros(nFreqs-computeForStaticCase,length(theta{j}));

            Z(m).eta{1,1} = zeros(nFreqs-computeForStaticCase,length(theta{j}));
            Z(m).eta{1,2} = zeros(nFreqs-computeForStaticCase,length(theta{j}));
            Z(m).eta{2,1} = zeros(nFreqs-computeForStaticCase,length(theta{j}));
            Z(m).eta{2,2} = zeros(nFreqs-computeForStaticCase,length(theta{j}));
            if useSymbolicPrecision                
                Z(m).xi{1,1} = vpa(Z(m).xi{1,1});
                Z(m).xi{1,2} = vpa(Z(m).xi{1,2});
                Z(m).xi{2,1} = vpa(Z(m).xi{2,1});
                Z(m).xi{2,2} = vpa(Z(m).xi{2,2});
                
                Z(m).eta{1,1} = vpa(Z(m).eta{1,1});
                Z(m).eta{1,2} = vpa(Z(m).eta{1,2});
                Z(m).eta{2,1} = vpa(Z(m).eta{2,1});
                Z(m).eta{2,2} = vpa(Z(m).eta{2,2});
            end
            if options.calc_u_r
                data(m).u_r = zeros(nFreqs,length(r{j}));
                if useSymbolicPrecision
                    data(m).u_r = vpa(data(m).u_r);
                end
            end
            if options.calc_u_t
                data(m).u_t = zeros(nFreqs,length(r{j}));
                if useSymbolicPrecision
                    data(m).u_t = vpa(data(m).u_t);
                end
            end
            if options.calc_sigma_rr
                data(m).sigma_rr = zeros(nFreqs,length(r{j}));
                if useSymbolicPrecision
                    data(m).sigma_rr = vpa(data(m).sigma_rr);
                end
            end
            if options.calc_sigma_tt
                data(m).sigma_tt = zeros(nFreqs,length(r{j}));
                if useSymbolicPrecision
                    data(m).sigma_tt = vpa(data(m).sigma_tt);
                end
            end
            if options.calc_sigma_pp
                data(m).sigma_pp = zeros(nFreqs,length(r{j}));
                if useSymbolicPrecision
                    data(m).sigma_pp = vpa(data(m).sigma_pp);
                end
            end
            if options.calc_sigma_rt
                data(m).sigma_rt = zeros(nFreqs,length(r{j}));
                if useSymbolicPrecision
                    data(m).sigma_rt = vpa(data(m).sigma_rt);
                end
            end
            if options.calc_navier
                data(m).navier1 = zeros(nFreqs,length(r{j}));
                data(m).navier2 = zeros(nFreqs,length(r{j}));
                if useSymbolicPrecision
                    data(m).navier1 = vpa(data(m).navier1);
                    data(m).navier2 = vpa(data(m).navier2);
                end
            end
        end
    end
    if ~mod(j,2)
        m = m + 1;
    end
end

indices = (1:nFreqs)';
Zindices = indices;
m = 1;
if computeForStaticCase
    if isa(options.P_inc,'function_handle')
        P_inc = options.P_inc(0);
    else
        P_inc = options.P_inc;
    end
    staticIdx = find(omega == 0);
    for j = 1:length(r)
        if ~isempty(r{j})
            if mod(j,2)
                if m > 1
                    data(m).p(staticIdx,:) = P_inc*ones(1,length(r{j}));
                end
            else
                A = -P_inc/(3*K);
                if options.calc_u_r
                    data(m).u_r(staticIdx,:) = A*repmat(r{j},nFreqs,1);
                end
                if options.calc_sigma_rr
                    data(m).sigma_rr(staticIdx,:) = A*ones(1,length(r{j}));
                end
                if options.calc_sigma_tt
                    data(m).sigma_tt(staticIdx,:) = A*ones(1,length(r{j}));
                end
                if options.calc_sigma_pp
                    data(m).sigma_pp(staticIdx,:) = A*ones(1,length(r{j}));
                end
            end
        end
        if ~mod(j,2)
            m = m + 1;
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
tiny = 1e-200; % To avoid dividing by zero.
countUpWards = 1;
if countUpWards
    n = 0;
else
    n = round(3*max(r{1})*max(k(:,1)));
end
data(1).N_eps = NaN(nFreqs,1);
flag = zeros(size(hasCnvrgd,1),1); % Program terminated successfully unless error occurs (for each frequency)
while n >= 0 && n <= N_max
    try % and hope that no spherical Bessel functions are evaluated to be too large
        omega_temp = omega(indices);
        k_temp = k(indices,:);
        if SHBC && M == 1
            a_temp = [];
            b_temp = [];
        else
            a_temp = a(indices,:);
            b_temp = b(indices,:);
        end
        
        CC = getCoeffs(n, omega_temp, a_temp, b_temp, k_temp, options);

        m = 1;
        hasCnvrgdTmp = zeros(length(indices),length(r)); % temporary hasCnvrgd matrix
        for j = 1:length(r)    
            hasCnvrgdTmp2 = ones(size(indices)); % temporary hasCnvrgd vector    
            if ~isempty(r{j})
                if countUpWards
                    [P{j}, dP{j}, d2P{j}] = legendreDerivs(n, cos(theta{j}), P{j}, dP{j}, d2P{j});
                else
                    P{j}(2,:) = legendre(n, cos(theta{j}));
                    dP{j}(2,:) = legendreDeriv(n, cos(theta{j}));
                    d2P{j}(2,:) = legendreDeriv2(n, cos(theta{j}));
                end
                if mod(j,2)
                    if m == 1
                        C = CC(:,1);
                    elseif m == M + 1
                        C = CC(:,end);
                    else
                        C = CC(:,6*(m-1):6*(m-1)+1);
                    end
                    zeta = k_temp(:,m)*r{j};
                    if ~options.calc_farField
                        if n == 0
                            Z(m).zeta{1,2} = bessel_s(n,zeta,1);
                            if m < M+1
                                Z(m).zeta{2,2} = bessel_s(n,zeta,2);
                            end
                        end
                        Z(m).zeta{1,1} = Z(m).zeta{1,2}(Zindices,:);
                        Z(m).zeta{1,2} = bessel_s(n+1,zeta,1);
                        if m < M+1
                            Z(m).zeta{2,1} = Z(m).zeta{2,2}(Zindices,:);
                            Z(m).zeta{2,2} = bessel_s(n+1,zeta,2);    
                        end
                    end
                    fluid = p_(m,n,zeta,theta{j},M,C,k_temp(:,m),P{j}(2,:), dP{j}(2,:), d2P{j}(2,:), Z(m).zeta, options);
                    if options.calc_p
                        data(m).p(indices,:) = data(m).p(indices,:) + fluid.p;
                        hasCnvrgdTmp2 = hasCnvrgdTmp2.*prod(abs(fluid.p)./(abs(data(m).p(indices,:))+tiny) < Eps,2);
                    end
                    if options.calc_dpdr
                        data(m).dpdr(indices,:) = data(m).dpdr(indices,:) + fluid.dpdr;
                        hasCnvrgdTmp2 = hasCnvrgdTmp2.*prod(abs(fluid.dpdr)./(abs(data(m).dpdr(indices,:))+tiny) < Eps,2);
                    end
                    if options.calc_dpdt
                        data(m).dpdt(indices,:) = data(m).dpdt(indices,:) + fluid.dpdt;
                        hasCnvrgdTmp2 = hasCnvrgdTmp2.*prod(abs(fluid.dpdt)./(abs(data(m).dpdt(indices,:))+tiny) < Eps,2);
                    end
                    if options.calc_d2pdr2
                        data(m).d2pdr2(indices,:) = data(m).d2pdr2(indices,:) + fluid.d2pdr2;
                        hasCnvrgdTmp2 = hasCnvrgdTmp2.*prod(abs(fluid.d2pdr2)./(abs(data(m).d2pdr2(indices,:))+tiny) < Eps,2);
                    end
                    if options.calc_d2pdt2
                        data(m).d2pdt2(indices,:) = data(m).d2pdt2(indices,:) + fluid.d2pdt2;
                        hasCnvrgdTmp2 = hasCnvrgdTmp2.*prod(abs(fluid.d2pdt2)./(abs(data(m).d2pdt2(indices,:))+tiny) < Eps,2);
                    end
                else
                    C_indices = 6*(m-1)+2:6*m-1;
                    
                    xi = a_temp(:,m)*r{j};
                    eta = b_temp(:,m)*r{j};
                    if n == 0
                        Z(m).xi{1,2} = bessel_s(n,xi,1);
                        Z(m).eta{1,2} = bessel_s(n,eta,1);
                        if ~(ESBC && m == M)
                            Z(m).xi{2,2} = bessel_s(n,xi,2);
                            Z(m).eta{2,2} = bessel_s(n,eta,2);
                        end
                    end     
                    Z(m).xi{1,1} = Z(m).xi{1,2}(Zindices,:);
                    Z(m).eta{1,1} = Z(m).eta{1,2}(Zindices,:);
                    Z(m).xi{1,2} = bessel_s(n+1,xi,1);
                    Z(m).eta{1,2} = bessel_s(n+1,eta,1);
                    if ~(ESBC && m == M)
                        Z(m).xi{2,1} = Z(m).xi{2,2}(Zindices,:);
                        Z(m).eta{2,1} = Z(m).eta{2,2}(Zindices,:);
                        Z(m).xi{2,2} = bessel_s(n+1,xi,2); 
                        Z(m).eta{2,2} = bessel_s(n+1,eta,2);
                    end
                    if ESBC && m == M
                        A = CC(:,C_indices(1));
                        B = CC(:,C_indices(3));
                    else
                        A = CC(:,C_indices(1:2));
                        B = CC(:,C_indices(3:4));
                    end
                    options.indices = indices;
                    solid = u_(m,n,r{j},theta{j},M,A,B,xi,eta,G(m),a_temp(:,m),b_temp(:,m),...
                                P{j}(2,:),dP{j}(2,:),d2P{j}(2,:), Z(m).xi, Z(m).eta, ESBC, options);

                    if options.calc_u_r
                        data(m).u_r(indices,:) = data(m).u_r(indices,:) + solid.u_r;
                        hasCnvrgdTmp2 = hasCnvrgdTmp2.*prod(abs(solid.u_r)./(abs(data(m).u_r(indices,:))+tiny) < Eps,2);
                    end
                    if options.calc_u_t
                        data(m).u_t(indices,:) = data(m).u_t(indices,:) + solid.u_t;
                        hasCnvrgdTmp2 = hasCnvrgdTmp2.*prod(abs(solid.u_t)./(abs(data(m).u_t(indices,:))+tiny) < Eps,2);
                    end
                    if options.calc_sigma_rr
                        data(m).sigma_rr(indices,:) = data(m).sigma_rr(indices,:) + solid.sigma_rr;
                        hasCnvrgdTmp2 = hasCnvrgdTmp2.*prod(abs(solid.sigma_rr)./(abs(data(m).sigma_rr(indices,:))+tiny) < Eps,2);
                    end
                    if options.calc_sigma_tt
                        data(m).sigma_tt(indices,:) = data(m).sigma_tt(indices,:) + solid.sigma_tt;
                        hasCnvrgdTmp2 = hasCnvrgdTmp2.*prod(abs(solid.sigma_tt)./(abs(data(m).sigma_tt(indices,:))+tiny) < Eps,2);
                    end
                    if options.calc_sigma_pp
                        data(m).sigma_pp(indices,:) = data(m).sigma_pp(indices,:) + solid.sigma_pp;
                        hasCnvrgdTmp2 = hasCnvrgdTmp2.*prod(abs(solid.sigma_pp)./(abs(data(m).sigma_pp(indices,:))+tiny) < Eps,2);
                    end
                    if options.calc_sigma_rt
                        data(m).sigma_rt(indices,:) = data(m).sigma_rt(indices,:) + solid.sigma_rt;
                        hasCnvrgdTmp2 = hasCnvrgdTmp2.*prod(abs(solid.sigma_rt)./(abs(data(m).sigma_rt(indices,:))+tiny) < Eps,2);
                    end
                    if options.calc_navier
                        data(m).navier1(indices,:) = data(m).navier1(indices,:) + solid.navier1;
                        hasCnvrgdTmp2 = hasCnvrgdTmp2.*prod(abs(solid.navier1)./(abs(data(m).navier1(indices,:))+tiny) < Eps,2);
                    end
                    if options.calc_navier
                        data(m).navier2(indices,:) = data(m).navier2(indices,:) + solid.navier2;
                        hasCnvrgdTmp2 = hasCnvrgdTmp2.*prod(abs(solid.navier2)./(abs(data(m).navier2(indices,:))+tiny) < Eps,2);
                    end
                end
            end
            hasCnvrgdTmp(logical(hasCnvrgdTmp2),j) = 1;
            
            if ~mod(j,2)
                m = m + 1;
            end
        end
        if countUpWards 
            hasCnvrgd(indices,:) = [hasCnvrgd(indices,2:end), prod(hasCnvrgdTmp,2)];
            indicesPrev = indices;
            indices = find(~prod(hasCnvrgd,2));
            [~,Zindices] = ismember(indices,indicesPrev);
            if length(indices) < length(indicesPrev)
                data(1).N_eps(setdiff(indicesPrev,indices)) = n;
            end
            if isempty(indices) % every element has converged
                break;
            end
            n = n + 1;
        else
            n = n - 1;
        end
    catch        
        flag = -~prod(hasCnvrgd,2);
        warning(['The summation ended prematurely at n = ' num2str(n) ...
                 ' because a Bessel function evaluation was too large or that the' ...
                 ' global matrix was singular to working precision.'])
        break
    end
end
data(1).flag = flag;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fluid = p_(m,n,zeta,theta,M,C,k,P,dP,d2P,Z,options)
% Note that in the case of m = M+1 and zeta = 0: 
% --- dpdz =: dpdr and dpdx = dpdy = dpdt = 0
% --- nabla p =: d2pdr2, d2pdt2 := 0

Q0 = P;
if options.calc_farField
    h_n   = repmat(1i^(-n-1)./k, 1, size(zeta,2));
    if options.calc_dpdr
        dh_n  = repmat(1i^(-n), size(zeta,1),size(zeta,2));
    end
    if options.calc_d2pdr2
        d2h_n = repmat(1i^(-n+1)*k, 1, size(zeta,2));
    end
else
    j_n = Z{1,1};
    if options.calc_dpdr
        dj_n = dbessel_s(n,zeta,1,Z);
    end
    if options.calc_d2pdr2
        d2j_n = d2bessel_s(n,zeta,1,Z);
    end
    if m == 1
        y_n = Z{2,1};
        
        h_n = j_n + 1i*y_n;
        if options.calc_dpdr
            dy_n = dbessel_s(n,zeta,2,Z);
            dh_n = dj_n + 1i*dy_n;
        end
        if options.calc_d2pdr2
            d2y_n = d2bessel_s(n,zeta,2,Z);
            d2h_n = d2j_n + 1i*d2y_n;
        end
    elseif m < M+1
        y_n = Z{2,1};
        if options.calc_dpdr
            dy_n = dbessel_s(n,zeta,2,Z);
        end
        if options.calc_d2pdr2
            d2y_n = d2bessel_s(n,zeta,2,Z);
        end
    end
end

if options.calc_dpdt
    Q1 = Q_(1,theta,P,dP,d2P,true);
end
if options.calc_d2pdt2
    Q2 = Q_(2,theta,P,dP,d2P);
end

if m == 1
    if options.calc_p
        fluid.p = C*Q0.*h_n;
    end
    if options.calc_dpdr
        fluid.dpdr = C.*k*Q0.*dh_n;
    end
    if options.calc_dpdt % calculate a scaled dpdt (scaled by csc)
        fluid.dpdt = C*Q1.*h_n;
    end
    if options.calc_d2pdr2
        fluid.d2pdr2 = C.*k.^2*Q0.*d2h_n;
    end
    if options.calc_d2pdt2
        fluid.d2pdt2 = C*Q2.*h_n;
    end
elseif m == M+1
    if options.calc_p
        fluid.p = C*Q0.*j_n;
    end
    if options.calc_dpdr
        fluid.dpdr = C.*k*Q0.*dj_n;
        indices = logical(zeta(1,:) < eps);
        if n == 1
            fluid.dpdr(:,indices) = repmat(k/3.*C,1,sum(indices));
        else
            fluid.dpdr(:,indices) = 0;
        end
    end
    if options.calc_dpdt % calculate a scaled dpdt (scaled by csc)
        fluid.dpdt = C*Q1.*j_n;
        indices = logical(zeta(1,:) < eps);
        fluid.dpdt(:,indices) = 0;
    end
    if options.calc_d2pdr2
        fluid.d2pdr2 = C.*k.^2*Q0.*d2j_n;
        indices = logical(zeta(1,:) < eps);
        if n == 0
            fluid.d2pdr2(:,indices) = repmat(-k.^2.*C,1,sum(indices));
        else
            fluid.d2pdr2(:,indices) = 0;
        end
    end
    if options.calc_d2pdt2
        fluid.d2pdt2 = C*Q2.*j_n;
        indices = logical(zeta(1,:) < eps);
        fluid.d2pdt2(:,indices) = 0;
    end
else
    if options.calc_p
        fluid.p = C(:,1)*Q0.*j_n + C(:,2)*Q0.*y_n;
    end
    if options.calc_dpdr
        fluid.dpdr = C(:,1).*k*Q0.*dj_n + C(:,2).*k*Q0.*dy_n;
    end
    if options.calc_dpdt % calculate a scaled dpdt (scaled by csc)
        fluid.dpdt = C(:,1)*Q1.*j_n + C(:,2)*Q1.*y_n;
    end
    if options.calc_d2pdr2
        fluid.d2pdr2 = C(:,1).*k.^2*Q0.*d2j_n + C(:,2).*k.^2*Q0.*d2y_n;
    end
    if options.calc_d2pdt2
        fluid.d2pdt2 = C(:,1)*Q2.*j_n + C(:,2)*Q2.*y_n;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function solid = u_(m,n,r,theta,M,A,B,xi,eta,G,a,b,P,dP,d2P,Zxi,Zeta,ESBC,options)
% Note that in the case of "ESBC" and r = 0: 
% --- u_z =: u_r and u_x = u_y = u_t = 0
% --- sigma_11 =: sigma_rr, sigma_22 =: sigma_tt, sigma_33 =: sigma_pp
% --- 0 =: sigma_rt

Q0 = Q_(0,theta,P,dP,d2P);
Q1s = Q_(1,theta,P,dP,d2P,true);
Q1 = sin(theta).*Q1s;
Q2 = Q_(2,theta,P,dP,d2P);

r2 = r.^2;

if options.calc_u_r
    Q0r = Q0./r;
    u_r =   A(:,1)*Q0r.*S_(1,1,n,xi,eta,Zxi) ...
          + B(:,1)*Q0r.*T_(1,1,n,eta,Zeta);
    if ~(ESBC && m == M)
        u_r = u_r + A(:,2)*Q0r.*S_(1,2,n,xi,eta,Zxi) ...
                  + B(:,2)*Q0r.*T_(1,2,n,eta,Zeta);
    end
    if ESBC
        indices = logical(r < eps);
        if n == 1
            u_r(:,indices) = repmat((a.*A(:,1) - 2*b.*B(:,1))/3,1,sum(indices));
        else
            u_r(:,indices) = 0;
        end
    end
    solid.u_r = u_r;
end

if options.calc_u_t
    Q1r = Q1./r;
    u_t =   A(:,1)*Q1r.*S_(2,1,n,xi,eta,Zxi) ...
          + B(:,1)*Q1r.*T_(2,1,n,eta,Zeta);
    if ~(ESBC && m == M)
        u_t = u_t + A(:,2)*Q1r.*S_(2,2,n,xi,eta,Zxi) ...
                  + B(:,2)*Q1r.*T_(2,2,n,eta,Zeta);
    end
    if ESBC
        u_t(:,logical(r < eps)) = 0;
    end
    solid.u_t = u_t;
end
if options.calc_sigma_rr
    Q0r2 = Q0./r2;
    sigma_rr =   A(:,1)*Q0r2.*S_(5,1,n,xi,eta,Zxi) ...
               + B(:,1)*Q0r2.*T_(5,1,n,eta,Zeta);
    if ~(ESBC && m == M)
        sigma_rr = sigma_rr + A(:,2)*Q0r2.*S_(5,2,n,xi,eta,Zxi) ...
                            + B(:,2)*Q0r2.*T_(5,2,n,eta,Zeta);
    end
	sigma_rr = 2*G*sigma_rr;
    if ESBC
        indices = logical(r < eps);
        if n == 0
            sigma_rr(:,indices) = G/15*repmat(5*(4*a.^2-3*b.^2).*A(:,1),1,sum(indices));
        elseif n == 2
            sigma_rr(:,indices) = G/15*repmat(-2*a.^2.*A(:,1) + 6*b.^2.*B(:,1),1,sum(indices));
        else
            sigma_rr(:,indices) = 0;
        end
    end
	solid.sigma_rr = sigma_rr;
end
if options.calc_sigma_tt
    Q0r2 = Q0./r2;
    Q2r2 = Q2./r2;
    sigma_tt =   A(:,1)*Q0r2.*S_(6,1,n,xi,eta,Zxi)  ...
               + B(:,1)*Q0r2.*T_(6,1,n,eta,Zeta) ...
               + A(:,1)*Q2r2.*S_(2,1,n,xi,eta,Zxi) ...
               + B(:,1)*Q2r2.*T_(2,1,n,eta,Zeta);
    if ~(ESBC && m == M)
        sigma_tt = sigma_tt + A(:,2)*Q0r2.*S_(6,2,n,xi,eta,Zxi) ...
                            + B(:,2)*Q0r2.*T_(6,2,n,eta,Zeta) ...
                            + A(:,2)*Q2r2.*S_(2,2,n,xi,eta,Zxi) ...
                            + B(:,2)*Q2r2.*T_(2,2,n,eta,Zeta);
    end
	sigma_tt = 2*G*sigma_tt;
    if ESBC
        indices = logical(r < eps);
        if n == 0
            sigma_tt(:,indices) = G/15*repmat(5*(4*a.^2-3*b.^2).*A(:,1),1,sum(indices));
        elseif n == 2
            sigma_tt(:,indices) = G/15*repmat(-2*a.^2.*A(:,1) + 6*b.^2.*B(:,1),1,sum(indices));
        else
            sigma_tt(:,indices) = 0;
        end
    end
	solid.sigma_tt = sigma_tt;
end
if options.calc_sigma_pp
    Q0r2 = Q0./r2;
    cott_Q1r2 = cos(theta).*Q1s./r2;
    sigma_pp =   A(:,1)*Q0r2.*S_(6,1,n,xi,eta,Zxi) ...
               + B(:,1)*Q0r2.*T_(6,1,n,eta,Zeta) ...
               + A(:,1)*cott_Q1r2.*S_(2,1,n,xi,eta,Zxi) ...
               + B(:,1)*cott_Q1r2.*T_(2,1,n,eta,Zeta);
    if ~(ESBC && m == M)
        sigma_pp = sigma_pp + A(:,2)*Q0r2.*S_(6,2,n,xi,eta,Zxi) ...
                            + B(:,2)*Q0r2.*T_(6,2,n,eta,Zeta) ...
                            + A(:,2)*cott_Q1r2.*S_(2,2,n,xi,eta,Zxi) ...
                            + B(:,2)*cott_Q1r2.*T_(2,2,n,eta,Zeta);
    end
	sigma_pp = 2*G*sigma_pp;
    if ESBC
        indices = logical(r < eps);
        if n == 0
            sigma_pp(:,indices) = G/15*repmat(5*(4*a.^2-3*b.^2).*A(:,1),1,sum(indices));
        elseif n == 2
            sigma_pp(:,indices) = G/15*repmat(4*a.^2.*A(:,1) - 12*b.^2.*B(:,1),1,sum(indices));
        else
            sigma_pp(:,indices) = 0;
        end
    end
	solid.sigma_pp = sigma_pp;
end
if options.calc_sigma_rt
    Q1r2 = Q1./r2;
    sigma_rt =   A(:,1)*Q1r2.*S_(7,1,n,xi,eta,Zxi) ...
               + B(:,1)*Q1r2.*T_(7,1,n,eta,Zeta);
    if ~(ESBC && m == M)
        sigma_rt = sigma_rt + A(:,2)*Q1r2.*S_(7,2,n,xi,eta,Zxi) ...
                            + B(:,2)*Q1r2.*T_(7,2,n,eta,Zeta);
    end
	sigma_rt = 2*G*sigma_rt;
    if ESBC
        sigma_rt(:,logical(r < eps)) = 0;
    end
	solid.sigma_rt = sigma_rt;
end
if options.calc_navier
    Q1sr2 = Q1s./r2;
    sigma_rt =   A(:,1)*Q1sr2.*S_(7,1,n,xi,eta,Zxi) ...
               + B(:,1)*Q1sr2.*T_(7,1,n,eta,Zeta);
    if ~(ESBC && m == M)
        sigma_rt = sigma_rt + A(:,2)*Q1sr2.*S_(7,2,n,xi,eta,Zxi) ...
                            + B(:,2)*Q1sr2.*T_(7,2,n,eta,Zeta);
    end
	sigma_rt = 2*G*sigma_rt;
    
    r3 = r.^3;
    Q0r3 = Q0./r3;
    Q1r3 = Q1./r3;
    Q2r3 = Q2./r3;
    dsigma_rr_dr =   A(:,1)*Q0r3.*S_(8,1,n,xi,eta,Zxi) ...
                   + B(:,1)*Q0r3.*T_(8,1,n,eta,Zeta);
    if ~(ESBC && m == M)
        dsigma_rr_dr =  dsigma_rr_dr + A(:,2)*Q0r3.*S_(8,2,n,xi,eta,Zxi) ...
                                     + B(:,2)*Q0r3.*T_(8,2,n,eta,Zeta);
    end
	dsigma_rr_dr = 2*G*dsigma_rr_dr;
      
    dsigma_rt_dr =   A(:,1)*Q1r3.*S_(9,1,n,xi,eta,Zxi) ...
                   + B(:,1)*Q1r3.*T_(9,1,n,eta,Zeta);
    if ~(ESBC && m == M)
        dsigma_rt_dr =  dsigma_rt_dr + A(:,2)*Q1r3.*S_(9,2,n,xi,eta,Zxi) ...
                                     + B(:,2)*Q1r3.*T_(9,2,n,eta,Zeta);
    end
	dsigma_rt_dr = 2*G*dsigma_rt_dr;

    dsigma_rt_dt =   A(:,1)*Q2r3.*S_(7,1,n,xi,eta,Zxi) ...
                   + B(:,1)*Q2r3.*T_(7,1,n,eta,Zeta);
    if ~(ESBC && m == M)
        dsigma_rt_dt = dsigma_rt_dt + A(:,2)*Q2r3.*S_(7,2,n,xi,eta,Zxi) ...
                                    + B(:,2)*Q2r3.*T_(7,2,n,eta,Zeta);
    end
	dsigma_rt_dt = 2*G*dsigma_rt_dt;

    dsigma_diffr =   A(:,1)*Q1r3.*S_(6,1,n,xi,eta,Zxi) ...
                   + B(:,1)*Q1r3.*T_(6,1,n,eta,Zeta) ...
                   + (-n^2-n+1)*A(:,1)*Q1r3.*S_(2,1,n,xi,eta,Zxi) ...
                   + (-n^2-n+1)*B(:,1)*Q1r3.*T_(2,1,n,eta,Zeta);
    if ~(ESBC && m == M)
        dsigma_diffr = dsigma_diffr + A(:,2)*Q1r3.*S_(6,2,n,xi,eta,Zxi) ...
                                    + B(:,2)*Q1r3.*T_(6,2,n,eta,Zeta) ...
                                    + (-n^2-n+1)*A(:,2)*Q1r3.*S_(2,2,n,xi,eta,Zxi) ...
                                    + (-n^2-n+1)*B(:,2)*Q1r3.*T_(2,2,n,eta,Zeta);
    end
	dsigma_diffr = 2*G*dsigma_diffr;
    
    R = repmat(r,size(xi,1),1);
    Theta = repmat(theta,size(xi,1),1);
    solid.navier1 = dsigma_rr_dr + dsigma_rt_dt + 1./R.*(2*sigma_rr-sigma_tt-sigma_pp+sigma_rt.*cos(Theta));
    solid.navier2 = dsigma_rt_dr + dsigma_diffr + 3./R.*sigma_rt.*sin(Theta);
    if ESBC && m == M
        indices = logical(r < eps);
        if n == 1
            omega = options.omega(options.indices);
            rho_s = options.rho_s(end);
            solid.navier1(:,indices) = repmat(-omega.^2*rho_s.*(a.*A(:,1) - 2*b.*B(:,1))/3,1,sum(indices));
        else
            solid.navier1(:,indices) = 0;
        end
        solid.navier2(:,indices) = 0;
    end
end
