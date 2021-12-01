function C = getCoeffs(n, omega, layer, options)

prec = options.prec;
SHBC = options.SHBC;
SSBC = options.SSBC;
IBC = options.IBC;
M = numel(layer);
m_s = options.m_s;
nu_a = options.nu_a;

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
    if m == 1
        Rt_m = layer{1}.R_i;
    elseif isSphere
        Rt_m = layer{m-1}.R_i;
    else
        Rt_m = (layer{m}.R_i+layer{m-1}.R_i)/2;
    end
    m2 = m + 1;
    if m2 <= M
        isSphere2 = layer{m2}.R_i == 0;
        if isSphere2
            Rt_m2 = layer{M-1}.R_i;
        else
            Rt_m2 = (layer{m2}.R_i+layer{m2-1}.R_i)/2;
        end
    end
    supportAtR_i = ~(singleLayerSupport && options.r_s ~= R_i) && ~isSphere;
    if strcmp(layer{m}.media,'fluid')
        k = layer{m}.k_temp;
        if ~strcmp(nextMedia,'origin')
            Z_zeta = layer{m}.Z_zeta_i;
            zeta_i = k.*R_i;
            gzeta_i  = {g_(n,1,zeta_i,nu_a),g_(n,2,zeta_i,nu_a),g_(n,3,zeta_i,nu_a)};
        end
        if IBC
            zs = options.z_temp./(1i.*k.*layer{m}.rho.*layer{m}.c_f);
        end
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
            D1{m} = D1_(n,k,R_i,Z_zeta,Z_r_s,rho,omega,options,gzeta_i);
        elseif IBC && M == m_s
            D1{m} = D2_(n,k,R_i,Z_zeta,Z_r_s,options) - reshape(zs,1,1,[]).*D1_(n,k,R_i,Z_zeta,Z_r_s,rho,omega,options,gzeta_i); % Ayres1987ars equation (38)
        else
            D1{m} = cat(1,D1_(n,k,R_i,Z_zeta,Z_r_s,rho,omega,options,gzeta_i),...    
                          D2_(n,k,R_i,Z_zeta,Z_r_s,options));
        end
    elseif m2 == m_s && supportAtR_i
        rho2 = layer{m2}.rho;
        if strcmp(options.applyLoad,'pointCharge')
            Z_r_s2 = layer{m2}.Z_r_s;
        else
            Z_r_s2 = NaN;
        end
        if singleLayerSupport
            k2 = NaN(size(layer{1}.k_temp)); % Not used
            Z_zeta2 = NaN;                 % Not used
        else
            k2 = layer{m2}.k_temp;
            Z_zeta2 = layer{m2}.Z_zeta_o;
        end
        if strcmp(layer{m}.media,'solid') || strcmp(layer{m}.media,'viscoelastic')
            void = zeros(n > 0,1,numel(omega),prec);
            D1{m} = cat(1,void,...
                          D2_(n,k2,R_i,Z_zeta2,Z_r_s2,options),...    
                          D1_(n,k2,R_i,Z_zeta2,Z_r_s2,rho2,omega,options,gzeta_i));
        else
            D1{m} = cat(1,D1_(n,k2,R_i,Z_zeta2,Z_r_s2,rho2,omega,options,gzeta_i),...    
                          D2_(n,k2,R_i,Z_zeta2,Z_r_s2,options));
        end
    end
    switch layer{m}.media
        case 'fluid'
            rho = layer{m}.rho;
            switch nextMedia
                case 'void'
                    if SHBC
                        H1{m} = dp_dr_s_(n,k,R_i,Rt_m,Z_zeta,rho,omega,isSphere,isOuterDomain,nu_a,gzeta_i);
                    elseif SSBC
                        H1{m} = p_(n,k,R_i,Rt_m,Z_zeta,isSphere,isOuterDomain,nu_a);
                    elseif IBC
                        H1{m} = p_(n,k,R_i,Rt_m,Z_zeta,isSphere,isOuterDomain,nu_a) + reshape(zs,1,1,[]).*dp_dr_s_(n,k,R_i,Rt_m,Z_zeta,rho,omega,isSphere,isOuterDomain,nu_a,gzeta_i);
                    end
                case 'fluid'
                        k2 = layer{m2}.k_temp;
                        rho2 = layer{m2}.rho;
                        Z_zeta2 = layer{m2}.Z_zeta_o;
                        zeta2 = k2.*R_i;
                        gzeta2  = {g_(n,1,zeta2,nu_a),g_(n,2,zeta2,nu_a)};
                        dp_dr_s  = dp_dr_s_(n,k, R_i,Rt_m,Z_zeta, rho, omega, isSphere,isOuterDomain,nu_a,gzeta_i);
                        dp_dr_s2 = dp_dr_s_(n,k2,R_i,Rt_m2,Z_zeta2,rho2,omega,isSphere2,false,nu_a,gzeta2);
                        p  = p_(n,k, R_i,Rt_m, Z_zeta, isSphere, isOuterDomain,nu_a);
                        p2 = p_(n,k2,R_i,Rt_m2,Z_zeta2,isSphere2,false,        nu_a);
                        
                        H1{m} = cat(1,cat(2, dp_dr_s, -dp_dr_s2),...    
                                      cat(2, p, -p2));
                case {'solid','viscoelastic'}
                        G2 = layer{m2}.G_temp;
                        a2 = layer{m2}.a_temp;
                        b2 = layer{m2}.b_temp;
                        Z_xi2  = layer{m2}.Z_xi_o;
                        Z_eta2 = layer{m2}.Z_eta_o;
                        xi2 = a2.*R_i;
                        eta2 = b2.*R_i;
                        gxi2  = {g_(n,1,xi2,nu_a),g_(n,2,xi2,nu_a)};
                        geta2  = {g_(n,1,eta2,nu_a),g_(n,2,eta2,nu_a)};
                        dp_dr_s = dp_dr_s_(n,k,R_i,Rt_m,Z_zeta,rho,omega,isSphere,isOuterDomain,nu_a,gzeta_i);
                        p = p_(n,k,R_i,Rt_m,Z_zeta,isSphere,isOuterDomain,nu_a);
                        u_r2 = u_r_(n,a2,b2,R_i,Rt_m2,Z_xi2,Z_eta2,isSphere2,nu_a,gxi2,geta2);
                        sigma_rr2 = sigma_rr_(n,a2,b2,R_i,Rt_m2,Z_xi2,Z_eta2,isSphere2,G2,nu_a,gxi2,geta2);
                        sigma_rt2 = sigma_rt_(n,a2,b2,R_i,Rt_m2,Z_xi2,Z_eta2,isSphere2,G2,nu_a,gxi2,geta2);
                        void = zeros(size(sigma_rt2,1),size(p,2),size(sigma_rt2,3),prec);
                        
                        H1{m} = cat(1,cat(2, dp_dr_s,u_r2),...
                                      cat(2, p, sigma_rr2),... 
                                      cat(2, void,sigma_rt2));
            end
            if isSphere || isOuterDomain
                dofs(m) = 1;
            else
                dofs(m) = 2;
            end
        case {'solid','viscoelastic'}
            G = layer{m}.G_temp;
            a = layer{m}.a_temp;
            b = layer{m}.b_temp;

            if ~strcmp(nextMedia,'origin')
                Z_xi  = layer{m}.Z_xi_i;
                Z_eta = layer{m}.Z_eta_i;
                xi_i = a.*R_i;
                eta_i = b.*R_i;
                gxi  = {g_(n,1,xi_i, nu_a),g_(n,2,xi_i, nu_a)};
                geta = {g_(n,1,eta_i,nu_a),g_(n,2,eta_i,nu_a)};
            end
            switch nextMedia
                case 'void'
                    if SHBC
                        u_r = u_r_(n,a,b,R_i,Rt_m,Z_xi,Z_eta,isSphere,nu_a,gxi,geta);
                        u_t = u_t_(n,a,b,R_i,Rt_m,Z_xi,Z_eta,isSphere,nu_a,gxi,geta);
                        
                        H1{m} = cat(1,u_r,u_t);
                    elseif SSBC
                        sigma_rr = sigma_rr_(n,a,b,R_i,Rt_m,Z_xi,Z_eta,isSphere,G,nu_a,gxi,geta);
                        sigma_rt = sigma_rt_(n,a,b,R_i,Rt_m,Z_xi,Z_eta,isSphere,G,nu_a,gxi,geta);
                        
                        H1{m} = cat(1,sigma_rr,sigma_rt);
                    elseif IBC
                        error('Impedance boundary conditions can only be implemented on a fluid media')
                    end
                case 'fluid'
                        k2 = layer{m2}.k_temp;
                        rho2 = layer{m2}.rho;
                        Z_zeta2 = layer{m2}.Z_zeta_o;
                        zeta2 = k2.*R_i;
                        gzeta2  = {g_(n,1,zeta2,nu_a),g_(n,2,zeta2,nu_a)};
                        dp_dr_s2 = dp_dr_s_(n,k2,R_i,Rt_m2,Z_zeta2,rho2,omega,isSphere2,false,nu_a,gzeta2);
                        p2 = p_(n,k2,R_i,Rt_m2,Z_zeta2,isSphere2,false,nu_a);
                        u_r = u_r_(n,a,b,R_i,Rt_m,Z_xi,Z_eta,isSphere,nu_a,gxi,geta);
                        sigma_rr = sigma_rr_(n,a,b,R_i,Rt_m,Z_xi,Z_eta,isSphere,G,nu_a,gxi,geta);
                        sigma_rt = sigma_rt_(n,a,b,R_i,Rt_m,Z_xi,Z_eta,isSphere,G,nu_a,gxi,geta);
                        void = zeros(size(sigma_rt,1),size(p2,2),numel(a),prec);
                        
                        H1{m} = cat(1,cat(2, sigma_rt,void),...
                                      cat(2, sigma_rr, p2),...
                                      cat(2, u_r,dp_dr_s2));
                case {'solid','viscoelastic'}
                        G2 = layer{m2}.G_temp;
                        a2 = layer{m2}.a_temp;
                        b2 = layer{m2}.b_temp;
                        Z_xi2  = layer{m2}.Z_xi_o;
                        Z_eta2 = layer{m2}.Z_eta_o;
                        xi2 = a2.*R_i;
                        eta2 = b2.*R_i;
                        gxi2  = {g_(n,1,xi2,nu_a),g_(n,2,xi2,nu_a)};
                        geta2  = {g_(n,1,eta2,nu_a),g_(n,2,eta2,nu_a)};
                        u_r  = u_r_(n,a, b, R_i,Rt_m, Z_xi, Z_eta, isSphere, nu_a,gxi, geta);
                        u_r2 = u_r_(n,a2,b2,R_i,Rt_m2,Z_xi2,Z_eta2,isSphere2,nu_a,gxi2,geta2);
                        u_t  = u_t_(n,a, b, R_i,Rt_m, Z_xi, Z_eta, isSphere, nu_a,gxi, geta);
                        u_t2 = u_t_(n,a2,b2,R_i,Rt_m2,Z_xi2,Z_eta2,isSphere2,nu_a,gxi2,geta2);
                        sigma_rr  = sigma_rr_(n,a, b, R_i,Rt_m, Z_xi, Z_eta, isSphere, G, nu_a,gxi,geta);
                        sigma_rr2 = sigma_rr_(n,a2,b2,R_i,Rt_m2,Z_xi2,Z_eta2,isSphere2,G2,nu_a,gxi2,geta2);
                        sigma_rt  = sigma_rt_(n,a, b, R_i,Rt_m, Z_xi, Z_eta, isSphere, G, nu_a,gxi,geta);
                        sigma_rt2 = sigma_rt_(n,a2,b2,R_i,Rt_m2,Z_xi2,Z_eta2,isSphere2,G2,nu_a,gxi2,geta2);
                        
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
% if n == 1100
%     keyboard
% end
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
H = zeros(systemSize,prec); % global matrix
D = zeros(systemSize,1,prec); % righ hand side
for j = 1:length(omega)
% parfor j = 1:length(omega)
%     H = zeros(systemSize,prec); % global matrix
%     D = zeros(systemSize,1,prec); % righ hand side
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
    H2 = H*Pinv;
    Pinv2 = diag(1./max(abs(H2),[],2));
    H2 = Pinv2*H2;
    CC(j,:) = diag(Pinv).*(H2\(Pinv2*D));
    if any(isinf(CC(j,:))) || any(isnan(CC(j,:)))
        warning('e3Dss:singularK','The modal matrix, K, was singular.')
    end
    
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

function D1 = D1_(n,k,R,Z,Z_r_s,rho,omega,options,g)
zeta = k*R;
switch options.applyLoad
    case 'planeWave'
        D1 = (2*n+1)*1i^n/R*dbessel_s(n,zeta,1,Z,true,g); 
    case 'pointCharge'
        r_s = options.r_s;
        if R < r_s
            D1 = (2*n+1)*1i/R*k.*dbessel_s(n,zeta,1,Z,true,g).*Z_r_s{3,1};
        else
            D1 = (2*n+1)*1i/R*k.*Z_r_s{1,1}.*dbessel_s(n,zeta,3,Z,true,g);
        end
    case 'mechExcitation'
        D1 = zeros(size(k),class(k));
    case 'surfExcitation'
        D1 = zeros(size(k),class(k));
    case 'radialPulsation'
        if n == 0
            D1 = -(1/R+1i*k);
        else
            D1 = zeros(size(k),class(k));
        end
    case 'custom'
        warning('This is experimental')
        D1 = zeros(size(k),class(k));
    otherwise
        error('Not implemented')
end
D1 = 1./(rho.*omega.^2).*D1;
D1 = reshape(D1,1,1,numel(D1));

function D2 = D2_(n,k,R,Z,Z_r_s,options)
switch options.applyLoad
    case 'planeWave'  
        D2 = (2*n+1)*1i^n*Z{1,1};
    case 'pointCharge'
        r_s = options.r_s;
        if R < r_s
            D2 = (2*n+1)*1i*k.*Z{1,1}.*Z_r_s{3,1};
        else
            D2 = (2*n+1)*1i*k.*Z_r_s{1,1}.*Z{3,1};
        end
    case 'mechExcitation'
        D2 = (2*n+1)/(4*pi*R^2)*ones(size(k),class(k));
    case 'surfExcitation'
        theta_s = options.theta_s;
        D2 = (  legendre_(n-1,cos(theta_s(2))) - legendre_(n+1,cos(theta_s(2))) ...
              -(legendre_(n-1,cos(theta_s(1))) - legendre_(n+1,cos(theta_s(1)))))/2*ones(size(k),class(k));
    case 'radialPulsation' 
        if n == 0
            D2 = ones(size(k),class(k));
        else
            D2 = zeros(size(k),class(k));
        end
    case 'custom' 
        warning('This is experimental')
        theta_s = options.theta_s;
        Psi = @(theta) exp(-theta.^2/theta_s(1));
        integrand = @(theta) Psi(theta).*legendre_(n,cos(theta)).*sin(theta);
        D2 = (2*n+1)/2*integral(integrand,0,pi)*ones(numel(k),1,class(R));
    otherwise
        error('Not implemented')
end
D2 = reshape(-D2,1,1,numel(D2));
                    
function H = p_(n,k,R,Rt_m,Z,isSphere,isOuterDomain,nu_a)

zeta = k*R;
zetat_m = k*Rt_m;
if isSphere
    H = zeros(1,1,length(k),class(R));
    H(1,1,:) = w_(n,1,zetat_m,zeta,nu_a,0).*Z{1,1}; 
elseif isOuterDomain
    H = zeros(1,1,length(k),class(R));
    H(1,1,:) = w_(n,3,zetat_m,zeta,nu_a,0).*Z{3,1}; 
else
    H = zeros(1,2,length(k),class(R));

    H(1,1,:) = w_(n,1,zetat_m,zeta,nu_a,0).*Z{1,1};
    H(1,2,:) = w_(n,2,zetat_m,zeta,nu_a,0).*Z{2,1};
end

function H = dp_dr_s_(n,k,R,Rt_m,Z,rho,omega,isSphere,isOuterDomain,nu_a,gzeta)

zeta = k*R;
zetat_m = k*Rt_m;
if isSphere
    H = zeros(1,1,length(k),class(R));
    zeta_dj_n = dbessel_s(n,zeta,1,Z,true,gzeta); % = zeta*dj_n
    temp = -1./(rho*omega.^2);
    H(1,1,:) = temp.*w_(n,1,zetat_m,zeta,nu_a,0).*zeta_dj_n;
elseif isOuterDomain
    H = zeros(1,1,length(k),class(R));
    zeta_dh_n = dbessel_s(n,zeta,3,Z,true,gzeta); % = zeta*dh_n
    temp = -1./(rho*omega.^2);
    
    H(1,1,:) = temp.*w_(n,3,zetat_m,zeta,nu_a,0).*zeta_dh_n; 
else
    H = zeros(1,2,length(k),class(R));

    zeta_dj_n = dbessel_s(n,zeta,1,Z,true,gzeta); % = zeta*dj_n
    zeta_dy_n = dbessel_s(n,zeta,2,Z,true,gzeta); % = zeta*dy_n
    temp = -1./(rho*omega.^2);
    H(1,1,:) = temp.*w_(n,1,zetat_m,zeta,nu_a,0).*zeta_dj_n;
    H(1,2,:) = temp.*w_(n,2,zetat_m,zeta,nu_a,0).*zeta_dy_n;
end
H = H/R;

function H = u_r_(n,a,b,R,Rt_m,Z_xi,Z_eta,isSphere,nu_a,gxi,geta)

xi = a*R;
eta = b*R;
xit_m = a*Rt_m;
etat_m = b*Rt_m;
if isSphere
    if n == 0
        H = zeros(1,1,length(a),class(R));

        H(1,1,:) = w_(n,1,xit_m,xi,nu_a,0).*S_(1, 1, n, xi, eta, Z_xi,gxi);
    else
        H = zeros(1,2,length(a),class(R));

        H(1,1,:) = w_(n,1,xit_m,xi,nu_a,0).*S_(1, 1, n, xi, eta, Z_xi,gxi);
        H(1,2,:) = w_(n,1,etat_m,eta,nu_a,0).*T_(1, 1, n, eta, Z_eta,geta);
    end
else
    if n == 0
        H = zeros(1,2,length(a),class(R));

        H(1,1,:) = w_(n,1,xit_m,xi,nu_a,0).*S_(1, 1, n, xi, eta, Z_xi,gxi);
        H(1,2,:) = w_(n,2,xit_m,xi,nu_a,0).*S_(1, 2, n, xi, eta, Z_xi,gxi);
    else
        H = zeros(1,4,length(a),class(R));

        H(1,1,:) = w_(n,1,xit_m,xi,nu_a,0).*S_(1, 1, n, xi, eta, Z_xi,gxi);
        H(1,2,:) = w_(n,2,xit_m,xi,nu_a,0).*S_(1, 2, n, xi, eta, Z_xi,gxi);
        H(1,3,:) = w_(n,1,etat_m,eta,nu_a,0).*T_(1, 1, n, eta, Z_eta,geta);
        H(1,4,:) = w_(n,2,etat_m,eta,nu_a,0).*T_(1, 2, n, eta, Z_eta,geta);
    end
end
H = H/R;

function H = u_t_(n,a,b,R,Rt_m,Z_xi,Z_eta,isSphere,nu_a,gxi,geta)

if n == 0
    if isSphere
        H = zeros(0,1,length(a),class(R));
    else
        H = zeros(0,2,length(a),class(R));
    end
else
    xi = a*R;
    eta = b*R;
    xit_m = a*Rt_m;
    etat_m = b*Rt_m;
    if isSphere
        H = zeros(1,2,length(a),class(R));

        H(1,1,:) = w_(n,1,xit_m,xi,nu_a,0).*S_(2, 1, n, xi, eta, Z_xi,gxi);
        H(1,2,:) = w_(n,1,etat_m,eta,nu_a,0).*T_(2, 1, n, eta, Z_eta,geta);
    else
        H = zeros(1,4,length(a),class(R));

        H(1,1,:) = w_(n,1,xit_m,xi,nu_a,0).*S_(2, 1, n, xi, eta, Z_xi,gxi);
        H(1,2,:) = w_(n,2,xit_m,xi,nu_a,0).*S_(2, 2, n, xi, eta, Z_xi,gxi);
        H(1,3,:) = w_(n,1,etat_m,eta,nu_a,0).*T_(2, 1, n, eta, Z_eta,geta);
        H(1,4,:) = w_(n,2,etat_m,eta,nu_a,0).*T_(2, 2, n, eta, Z_eta,geta);
    end
    H = H/R;
end

function H = sigma_rr_(n,a,b,R,Rt_m,Z_xi,Z_eta,isSphere,G,nu_a,gxi,geta)

xi = a*R;
eta = b*R;
xit_m = a*Rt_m;
etat_m = b*Rt_m;
if isSphere
    if n == 0
        H = zeros(1,1,length(a),class(R));

        H(1,1,:) = w_(n,1,xit_m,xi,nu_a,0).*S_(5, 1, n, xi, eta, Z_xi,gxi);
    else
        H = zeros(1,2,length(a),class(R));

        H(1,1,:) = w_(n,1,xit_m,xi,nu_a,0).*S_(5, 1, n, xi, eta, Z_xi,gxi);
        H(1,2,:) = w_(n,1,etat_m,eta,nu_a,0).*T_(5, 1, n, eta, Z_eta,geta);
    end
else
    if n == 0
        H = zeros(1,2,length(a),class(R));

        H(1,1,:) = w_(n,1,xit_m,xi,nu_a,0).*S_(5, 1, n, xi, eta, Z_xi,gxi);
        H(1,2,:) = w_(n,2,xit_m,xi,nu_a,0).*S_(5, 2, n, xi, eta, Z_xi,gxi);
    else
        H = zeros(1,4,length(a),class(R));

        H(1,1,:) = w_(n,1,xit_m,xi,nu_a,0).*S_(5, 1, n, xi, eta, Z_xi,gxi);
        H(1,2,:) = w_(n,2,xit_m,xi,nu_a,0).*S_(5, 2, n, xi, eta, Z_xi,gxi);

        H(1,3,:) = w_(n,1,etat_m,eta,nu_a,0).*T_(5, 1, n, eta, Z_eta,geta);
        H(1,4,:) = w_(n,2,etat_m,eta,nu_a,0).*T_(5, 2, n, eta, Z_eta,geta);
    end
end
H = 2*reshape(G,1,1,[])/R^2.*H;

function H = sigma_rt_(n,a,b,R,Rt_m,Z_xi,Z_eta,isSphere,G,nu_a,gxi,geta)

if n == 0
    if isSphere
        H = zeros(0,1,length(a),class(R));
    else
        H = zeros(0,2,length(a),class(R));
    end
else
    xi = a*R;
    eta = b*R;
    xit_m = a*Rt_m;
    etat_m = b*Rt_m;
    if isSphere
        H = zeros(1,2,length(a),class(R));
        
        H(1,1,:) = w_(n,1,xit_m,xi,nu_a,0).*S_(7, 1, n, xi, eta, Z_xi,gxi);
        H(1,2,:) = w_(n,2,etat_m,eta,nu_a,0).*T_(7, 1, n, eta, Z_eta,geta);
    else
        H = zeros(1,4,length(a),class(R));
        
        H(1,1,:) = w_(n,1,xit_m,xi,nu_a,0).*S_(7, 1, n, xi, eta, Z_xi,gxi);
        H(1,2,:) = w_(n,2,xit_m,xi,nu_a,0).*S_(7, 2, n, xi, eta, Z_xi,gxi);
        H(1,3,:) = w_(n,1,etat_m,eta,nu_a,0).*T_(7, 1, n, eta, Z_eta,geta);
        H(1,4,:) = w_(n,2,etat_m,eta,nu_a,0).*T_(7, 2, n, eta, Z_eta,geta);
    end
    H = 2*reshape(G,1,1,[])/R^2.*H;
end







