function C = getCoeffs(n, omega, layer, options)

prec = options.prec;
SHBC = options.SHBC;
SSBC = options.SSBC;
IBC = options.IBC;
M = numel(layer);
m_s = options.m_s;

%% Calculate submatrices
H1 = cell(M,1);
D1 = cell(M,1);
singleLayerSupport = strcmp(options.applyLoad,'mechExcitation') || strcmp(options.applyLoad,'surfExcitation') || strcmp(options.applyLoad,'custom');

dofs = zeros(M,1);
if M > 1
    nextMedia = layer{2}.media;
else
    nextMedia = 'void';
end
for m = 1:M
    R_i = layer{m}.R_i;
    isSphere = R_i == 0;
    isOuterDomain = m == 1;
    supportAtR_i = ~(singleLayerSupport && options.r_s ~= R_i) && ~isSphere;
    if IBC
        zs = options.z/1i*reshape(omega,1,1,numel(omega));
    end
    if m == m_s && supportAtR_i
        rho = layer{m}.rho;
        if singleLayerSupport
            k = NaN(size(layer{1}.k_temp)); % Not used
            Z_zeta = NaN;                   % Not used
        else
            k = layer{m}.k_temp;
            Z_zeta = layer{m}.Z_zeta_i;
        end
        if strcmp(options.applyLoad,'pointCharge')
            Z_r_s = layer{m}.Z_r_s;
        else
            Z_r_s = NaN;
        end
        if SSBC && M == m_s
            D1{m} = D2_(n,k,R_i,Z_zeta,Z_r_s,options);
        elseif SHBC && M == m_s
            D1{m} = D1_(n,k,R_i,Z_zeta,Z_r_s,rho,omega,options);
        elseif IBC && M == m_s
            D1{m} = -D2_(n,k,R_i,Z_zeta,Z_r_s,options) + zs.*D1_(n,k,R_i,Z_zeta,Z_r_s,rho,omega,options); % Ayres1987ars equation (38)
        else
            D1{m} = cat(1,D1_(n,k,R_i,Z_zeta,Z_r_s,rho,omega,options),...    
                          D2_(n,k,R_i,Z_zeta,Z_r_s,options));
        end
    elseif m+1 == m_s && supportAtR_i
        rho = layer{m+1}.rho;
        if strcmp(options.applyLoad,'pointCharge')
            Z_r_s = layer{m+1}.Z_r_s;
        else
            Z_r_s = NaN;
        end
        if singleLayerSupport
            k = NaN(size(layer{1}.k_temp)); % Not used
            Z_zeta = NaN;                   % Not used
        else
            k = layer{m+1}.k_temp;
            Z_zeta = layer{m+1}.Z_zeta_o;
        end
        if strcmp(layer{m}.media,'solid') || strcmp(layer{m}.media,'viscoelastic')
            void = zeros(n > 0,1,numel(omega),class(R_i));
            D1{m} = cat(1,void,...
                          D2_(n,k,R_i,Z_zeta,Z_r_s,options),...    
                          D1_(n,k,R_i,Z_zeta,Z_r_s,rho,omega,options));
        else
            D1{m} = cat(1,D1_(n,k,R_i,Z_zeta,Z_r_s,rho,omega,options),...    
                          D2_(n,k,R_i,Z_zeta,Z_r_s,options));
        end
    end
    switch layer{m}.media
        case 'fluid'
            k = layer{m}.k_temp;
            rho = layer{m}.rho;
            if ~strcmp(nextMedia,'origin')
                Z_zeta = layer{m}.Z_zeta_i;
            end
            switch nextMedia
                case 'void'
                    if SHBC
                        H1{m} = dp_dr_s_(n,rho,k,omega,R_i,Z_zeta,isSphere,isOuterDomain);
                    elseif SSBC
                        H1{m} = p_(Z_zeta,isSphere,isOuterDomain);
                    elseif IBC
                        H1{m} = p_(Z_zeta,isSphere,isOuterDomain) - zs.*dp_dr_s_(n,rho,k,omega,R_i,Z_zeta,isSphere,isOuterDomain);
                    end
                case 'fluid'
                        k2 = layer{m+1}.k_temp;
                        rho2 = layer{m+1}.rho;
                        Z_zeta2 = layer{m+1}.Z_zeta_o;
                        isSphere2 = layer{m+1}.R_i == 0;
                        dp_dr_s = dp_dr_s_(n,rho, k, omega,R_i,Z_zeta, isSphere,isOuterDomain);
                        dp_dr_s2 = dp_dr_s_(n,rho2,k2,omega,R_i,Z_zeta2,isSphere2,0);
                        p = p_(Z_zeta, isSphere,isOuterDomain);
                        p2 = p_(Z_zeta2,isSphere2,0);
                        
                        H1{m} = cat(1,cat(2, dp_dr_s, -dp_dr_s2),...    
                                      cat(2, p, -p2));
                case {'solid','viscoelastic'}
                        isSphere2 = layer{m+1}.R_i == 0;
                        G = layer{m+1}.G;
                        a = layer{m+1}.a_temp;
                        b = layer{m+1}.b_temp;
                        Z_xi  = layer{m+1}.Z_xi_o;
                        Z_eta = layer{m+1}.Z_eta_o;
                        dp_dr_s = dp_dr_s_(n,rho,k,omega,R_i,Z_zeta,isSphere,isOuterDomain);
                        p = p_(Z_zeta, isSphere,isOuterDomain);
                        u_r = u_r_(n,a,b,R_i,Z_xi,Z_eta,isSphere2);
                        sigma_rr = sigma_rr_(n,a,b,R_i,Z_xi,Z_eta,isSphere2,G);
                        sigma_rt = sigma_rt_(n,a,b,R_i,Z_xi,Z_eta,isSphere2,G);
                        void = zeros(size(sigma_rt,1),size(p,2),size(sigma_rt,3),class(a));
                        
                        H1{m} = cat(1,cat(2, dp_dr_s,u_r),...
                                      cat(2, p, sigma_rr),... 
                                      cat(2, void,sigma_rt));
            end
            if isSphere || isOuterDomain
                dofs(m) = 1;
            else
                dofs(m) = 2;
            end
        case {'solid','viscoelastic'}
            G = layer{m}.G;
            a = layer{m}.a_temp;
            b = layer{m}.b_temp;

            if ~strcmp(nextMedia,'origin')
                Z_xi  = layer{m}.Z_xi_i;
                Z_eta = layer{m}.Z_eta_i;
            end
            switch nextMedia
                case 'void'
                    if SHBC
                        u_r = u_r_(n,a,b,R_i,Z_xi,Z_eta,isSphere);
                        u_t = u_t_(n,a,b,R_i,Z_xi,Z_eta,isSphere);
                        
                        H1{m} = cat(1,u_r,u_t);
                    elseif SSBC
                        sigma_rr = sigma_rr_(n,a,b,R_i,Z_xi,Z_eta,isSphere,G);
                        sigma_rt = sigma_rt_(n,a,b,R_i,Z_xi,Z_eta,isSphere,G);
                        
                        H1{m} = cat(1,sigma_rr,sigma_rt);
                    elseif IBC
                        error('Impedance boundary conditions can only be implemented on a fluid media')
                    end
                case 'fluid'
                        k = layer{m+1}.k_temp;
                        rho = layer{m+1}.rho;
                        isSphere2 = layer{m+1}.R_i == 0;
                        Z_zeta = layer{m+1}.Z_zeta_o;
                        dp_dr_s = dp_dr_s_(n,rho,k,omega,R_i,Z_zeta,isSphere2,0);
                        p = p_(Z_zeta, isSphere2,0);
                        u_r = u_r_(n,a,b,R_i,Z_xi,Z_eta,isSphere);
                        sigma_rr = sigma_rr_(n,a,b,R_i,Z_xi,Z_eta,isSphere,G);
                        sigma_rt = sigma_rt_(n,a,b,R_i,Z_xi,Z_eta,isSphere,G);
                        void = zeros(size(sigma_rt,1),size(p,2),numel(a),class(a));
                        
                        H1{m} = cat(1,cat(2, sigma_rt,void),...
                                      cat(2, sigma_rr, p),...
                                      cat(2, u_r,dp_dr_s));
                case {'solid','viscoelastic'}
                        G2 = layer{m+1}.G;
                        a2 = layer{m+1}.a_temp;
                        b2 = layer{m+1}.b_temp;
                        isSphere2 = layer{m+1}.R_i == 0;
                        Z_xi2  = layer{m+1}.Z_xi_o;
                        Z_eta2 = layer{m+1}.Z_eta_o;
                        u_r = u_r_(n,a, b, R_i,Z_xi, Z_eta, isSphere);
                        u_r2 = u_r_(n,a2,b2,R_i,Z_xi2,Z_eta2,isSphere2);
                        u_t = u_t_(n,a, b, R_i,Z_xi, Z_eta, isSphere);
                        u_t2 = u_t_(n,a2,b2,R_i,Z_xi2,Z_eta2,isSphere2);
                        sigma_rr = sigma_rr_(n,a, b, R_i,Z_xi, Z_eta, isSphere,G);
                        sigma_rr2 = sigma_rr_(n,a2,b2,R_i,Z_xi2,Z_eta2,isSphere2,G2);
                        sigma_rt = sigma_rt_(n,a, b, R_i,Z_xi, Z_eta, isSphere,G);
                        sigma_rt2 = sigma_rt_(n,a2,b2,R_i,Z_xi2,Z_eta2,isSphere2,G2);
                        
                        H1{m} = cat(1,cat(2, u_r,-u_r2),...
                                      cat(2, sigma_rr,-sigma_rr2),...
                                      cat(2, sigma_rt,-sigma_rt2),...
                                      cat(2, u_t,-u_t2));
            end
            if isSphere || isOuterDomain
                if n == 0
                    dofs(m) = 1;
                else
                    dofs(m) = 2;
                end
            else
                if n == 0
                    dofs(m) = 2;
                else
                    dofs(m) = 4;
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
%% Calculate coefficients CC for each frequency
C = cell(M,1);
if M == 1
    if SSBC || SHBC || IBC
        C{1} = D1{1}(:)./H1{1}(:);
    else
        error('Not implemented')
    end
    return
end
systemSize = sum(dofs);
CC = zeros(length(omega),systemSize,prec);

% The matrix H and vector D should be allocated outside the for-loop if parfor is not used
% H = zeros(systemSize,prec); % global matrix
% D = zeros(systemSize,1,prec); % righ hand side
for j = 1:length(omega)
% parfor j = 1:length(omega)
    H = zeros(systemSize,prec); % global matrix
    D = zeros(systemSize,1,prec); % righ hand side
    I = 1;
    J = 1;
    for m = 1:M
        if ~isempty(H1{m})
            H_temp = H1{m}(:,:,j);
            ii = size(H_temp,1);
            jj = size(H_temp,2);
            H(I:I+ii-1,J:J+jj-1) = H_temp;
            if ~isempty(D1{m})
                ii2 = size(D1{m},1);
                D(I:I+ii2-1) = D1{m}(:,:,j);
            end
            I = I + ii;
            J = J + dofs(m);
        end
    end
    Pinv = diag(1./max(abs(H)));
    if any(isinf(Pinv(:)))
        error('e3Dss:singularK','K was singular')
    end
    H2 = H*Pinv;
    Pinv2 = diag(1./max(abs(H2),[],2));
    H2 = Pinv2*H2;
    CC(j,:) = diag(Pinv).*(H2\(Pinv2*D));
    
    % Uncomment the following to get the spy matrix in the paper
%     if n == 300
%         fileName = 'results/spy_H';
%         figure(1)
%         spy2(H)
%         cond(full(H))
%         extraAxisOptions = {...
%             'axis on top=true', ...
%             'at={(0,0)}', ...
%             'xtick={1,2,...,18}', ...
%             'ytick={1,2,...,18}', ...
%             'xlabel=$j$', ...
%             'ylabel=$i$', ...
%             'colorbar style={ylabel={$H_{ij,300}$}, ytick={-300,-200,...,200}, yticklabels={$10^{-300}$, $10^{-200}$, $10^{-100}$, $10^{0}$, $10^{100}$, $10^{200}$}}'};
%         
% %         matlab2tikz([fileName '_1.tex'], 'height', '3.2094in', ...
% %             'extraAxisOptions', extraAxisOptions, 'relativeDataPath', '../../../results/S135') % , 'imagesAsPng', false
%         figure(2)
%         spy2(Pinv2*H*Pinv)
%         cond(full(Pinv2*H*Pinv))
%         extraAxisOptions{7} = 'colorbar style={ylabel={$\tilde{H}_{ij,300}$}, ytick={-10,-8,...,0}, yticklabels={$10^{-10}$, $10^{-8}$, $10^{-6}$, $10^{-4}$, $10^{-2}$, $10^0$}}';
% %         matlab2tikz([fileName '_2.tex'], 'height', '3.2094in', ...
% %             'extraAxisOptions', extraAxisOptions, 'relativeDataPath', '../../../results/S135')
%         keyboard
%     end
end
J = 1;
for m = 1:M
    C{m} = CC(:,J:J+dofs(m)-1);
    J = J + dofs(m);
end

function D1 = D1_(n,k,R,Z,Z_r_s,rho,omega,options)
zeta = k*R;
switch options.applyLoad
    case 'planeWave'
        D1 = (2*n+1)*1i^n/R*(n*Z{1,1} - zeta.*Z{1,2});
    case 'pointCharge'
        r_s = options.r_s;
        if R < r_s
            D1 = (2*n+1)*1i/R*k.*(n*Z{1,1} - zeta.*Z{1,2}).*(Z_r_s{1,1}+1i*Z_r_s{2,1});
        else
            D1 = (2*n+1)*1i/R*k.*Z_r_s{1,1}.*(n*(Z{1,1}+1i*Z{2,1}) - zeta.*(Z{1,2}+1i*Z{2,2}));
        end
    case 'mechExcitation'
        D1 = zeros(size(k));
    case 'surfExcitation'
        D1 = zeros(size(k));
    case 'radialPulsation'
        if n == 0
            D1 = -(1/R+1i*k);
        else
            D1 = zeros(size(k));
        end
    case 'custom'
        D1 = zeros(size(k));
end
D1 = 1./(rho*omega.^2).*D1;
D1 = reshape(D1,1,1,numel(D1));

function D2 = D2_(n,k,R,Z,Z_r_s,options)
switch options.applyLoad
    case 'planeWave'  
        D2 = (2*n+1)*1i^n.*Z{1,1};
    case 'pointCharge'
        r_s = options.r_s;
        if R < r_s
            D2 = (2*n+1)*1i*k.*Z{1,1}.*(Z_r_s{1,1}+1i*Z_r_s{2,1});
        else
            D2 = (2*n+1)*1i*k.*Z_r_s{1,1}.*(Z{1,1}+1i*Z{2,1});
        end
    case 'mechExcitation'
        D2 = (2*n+1)/(4*pi*R^2)*ones(size(k));
    case 'surfExcitation'
        theta_s = options.theta_s;
        D2 = (  legendre_(n-1,cos(theta_s(2))) - legendre_(n+1,cos(theta_s(2))) ...
              -(legendre_(n-1,cos(theta_s(1))) - legendre_(n+1,cos(theta_s(1)))))/2*ones(size(k));
    case 'radialPulsation' 
        if n == 0
            D2 = ones(size(k));
        else
            D2 = zeros(size(k));
        end
    case 'custom' 
        theta_s = options.theta_s;
        Psi = @(theta) exp(-theta.^2/theta_s(1));
        integrand = @(theta) Psi(theta).*legendre_(n,cos(theta)).*sin(theta);
        D2 = (2*n+1)/2*integral(integrand,0,pi)*ones(numel(k),1,class(R));
end
D2 = reshape(-D2,1,1,numel(D2));
                    
function H = p_(Z,isSphere,isOuterDomain)

if isSphere
    H = zeros(1,1,length(Z{1,1}),class(Z{1,1}));
    H(1,1,:) = Z{1,1};
elseif isOuterDomain
    H = zeros(1,1,length(Z{1,1}),class(Z{1,1}));
    H(1,1,:) = Z{1,1}+1i*Z{2,1};   
else
    H = zeros(1,2,length(Z{1,1}),class(Z{1,1}));

    H(1,1,:) = Z{1,1};
    H(1,2,:) = Z{2,1};
end

function H = dp_dr_s_(n,rho,k,omega,R,Z,isSphere,isOuterDomain)

zeta = k*R;
if isSphere
    H = zeros(1,1,length(k),class(R));
    zeta_dj_n = n*Z{1,1} - zeta.*Z{1,2}; % = zeta*dj_n
    temp = -1./(rho*omega.^2);
    H(1,1,:) = temp.*zeta_dj_n;
elseif isOuterDomain
    H = zeros(1,1,length(k),class(R));
    zeta_dj_n = n*(Z{1,1}+1i*Z{2,1}) - zeta.*(Z{1,2}+1i*Z{2,2}); % = zeta*dj_n
    temp = -1./(rho*omega.^2);
    H(1,1,:) = temp.*zeta_dj_n;
else
    H = zeros(1,2,length(k),class(R));

    zeta_dj_n = n*Z{1,1} - zeta.*Z{1,2}; % = zeta*dj_n
    zeta_dy_n = n*Z{2,1} - zeta.*Z{2,2}; % = zeta*dy_n
    temp = -1./(rho*omega.^2);
    H(1,1,:) = temp.*zeta_dj_n;
    H(1,2,:) = temp.*zeta_dy_n;
end
H = H/R;

function H = u_r_(n,a,b,R,Z_xi,Z_eta,isSphere)

xi = a*R;
eta = b*R;
if isSphere
    if n == 0
        H = zeros(1,1,length(a),class(R));

        H(1,1,:) = S_(1, 1, n, xi, eta, Z_xi);
    else
        H = zeros(1,2,length(a),class(R));

        H(1,1,:) = S_(1, 1, n, xi, eta, Z_xi);
        H(1,2,:) = T_(1, 1, n, eta, Z_eta);
    end
else
    if n == 0
        H = zeros(1,2,length(a),class(R));

        H(1,1,:) = S_(1, 1, n, xi, eta, Z_xi);
        H(1,2,:) = S_(1, 2, n, xi, eta, Z_xi);
    else
        H = zeros(1,4,length(a),class(R));

        H(1,1,:) = S_(1, 1, n, xi, eta, Z_xi);
        H(1,2,:) = S_(1, 2, n, xi, eta, Z_xi);
        H(1,3,:) = T_(1, 1, n, eta, Z_eta);
        H(1,4,:) = T_(1, 2, n, eta, Z_eta);
    end
end
H = H/R;

function H = u_t_(n,a,b,R,Z_xi,Z_eta,isSphere)

if n == 0
    if isSphere
        H = zeros(0,1,length(a),class(R));
    else
        H = zeros(0,2,length(a),class(R));
    end
else
    xi = a*R;
    eta = b*R;
    if isSphere
        H = zeros(1,2,length(a),class(R));

        H(1,1,:) = S_(2, 1, n, xi, eta, Z_xi);
        H(1,2,:) = T_(2, 1, n, eta, Z_eta);
    else
        H = zeros(1,4,length(a),class(R));

        H(1,1,:) = S_(2, 1, n, xi, eta, Z_xi);
        H(1,2,:) = S_(2, 2, n, xi, eta, Z_xi);
        H(1,3,:) = T_(2, 1, n, eta, Z_eta);
        H(1,4,:) = T_(2, 2, n, eta, Z_eta);
    end
    H = H/R;
end

function H = sigma_rr_(n,a,b,R,Z_xi,Z_eta,isSphere,G)

xi = a*R;
eta = b*R;
if isSphere
    if n == 0
        H = zeros(1,1,length(a),class(R));

        H(1,1,:) = S_(5, 1, n, xi, eta, Z_xi);
    else
        H = zeros(1,2,length(a),class(R));

        H(1,1,:) = S_(5, 1, n, xi, eta, Z_xi);
        H(1,2,:) = T_(5, 1, n, eta, Z_eta);
    end
else
    if n == 0
        H = zeros(1,2,length(a),class(R));

        H(1,1,:) = S_(5, 1, n, xi, eta, Z_xi);
        H(1,2,:) = S_(5, 2, n, xi, eta, Z_xi);
    else
        H = zeros(1,4,length(a),class(R));

        H(1,1,:) = S_(5, 1, n, xi, eta, Z_xi);
        H(1,2,:) = S_(5, 2, n, xi, eta, Z_xi);

        H(1,3,:) = T_(5, 1, n, eta, Z_eta);
        H(1,4,:) = T_(5, 2, n, eta, Z_eta);
    end
end
H = 2*G/R^2*H;

function H = sigma_rt_(n,a,b,R,Z_xi,Z_eta,isSphere,G)

if n == 0
    if isSphere
        H = zeros(0,1,length(a),class(R));
    else
        H = zeros(0,2,length(a),class(R));
    end
else
    xi = a*R;
    eta = b*R;
    if isSphere
        H = zeros(1,2,length(a),class(R));

        H(1,1,:) = S_(7, 1, n, xi, eta, Z_xi);
        H(1,2,:) = T_(7, 1, n, eta, Z_eta);
    else
        H = zeros(1,4,length(a),class(R));
        H(1,1,:) = S_(7, 1, n, xi, eta, Z_xi);
        H(1,2,:) = S_(7, 2, n, xi, eta, Z_xi);
        H(1,3,:) = T_(7, 1, n, eta, Z_eta);
        H(1,4,:) = T_(7, 2, n, eta, Z_eta);
    end
    H = 2*G/R^2*H;
end








