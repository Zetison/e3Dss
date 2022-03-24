function media = getTerm_n(layer,m,n,xi,eta,A,B,Rt_m,omega,isSphere,isOuterDomain,applyLoad,besselIndices,nu_a,exponentShift,options)
% Note that in the case of isSphere and r = 0: 
% --- dpdz =: dpdr and dpdx = dpdy = dpdt = 0, with same convention for p_inc
% --- nabla p =: d2pdr2, d2pdt2 := 0
% --- u_z =: u_r and u_x = u_y = u_t = 0
% --- du_xdx =: du_rdr, du_ydy =: du_rdt, du_zdz =: du_tdr, 0 =: du_tdt
% --- sigma_11 =: sigma_rr, sigma_22 =: sigma_tt, sigma_33 =: sigma_pp, 0 =: sigma_rt
% Also note that dpdt and dp_incdt are scaled by csc(theta)
% Finally note that u_t, du_tdr and du_rdt are scaled by csc(theta)
media = [];

rho = layer{m}.rho;
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

if layer{m}.calc_p || layer{m}.calc_dpdr || layer{m}.calc_dpdt || layer{m}.calc_d2pdr2 || layer{m}.calc_d2pdr2 ...
        || layer{m}.calc_p_inc || layer{m}.calc_dp_incdr || layer{m}.calc_dp_incdt
    w = cell(1,3);
    wZt = cell(1,3);
    dwZt = cell(1,3);
    d2wZt = cell(1,3);
    for i = 1:numel(besselIndices)
        if besselIndices(i)
            wZt{i} = Zxi{i,1};
            if layer{m}.calc_dpdr
                dwZt{i} = dbessel_s(n,zeta,i,Zxi,false,gxi);
            end
            if layer{m}.calc_d2pdr2
                d2wZt{i} = d2bessel_s(n,zeta,i,Zxi,gxi);
            end
        end
    end
    if layer{m}.calc_p_inc && options.p_inc_fromSeries
        switch applyLoad
            case 'planeWave'
                s1 = exp(exponent_(1,nu,zeta,nu_a));
                media.p_inc = (2*n+1)*1i^n*Q0.*wZt{1}./s1;
            case 'pointCharge'
                s1 = exp(exponent_(1,nu,zeta,nu_a));
                s3 = exp(exponent_(3,nu,zeta,nu_a));
                Z_r_s = layer{m}.Z_r_s;
                r_s = options.r_s;
                s1_r_s = exp(exponent_(1,nu,a*r_s,nu_a));
                s3_r_s = exp(exponent_(3,nu,a*r_s,nu_a));
                media.p_inc = zeros(size(zeta),class(zeta));
                indices = repmat(r < r_s,numel(a),1);
                temp = (2*n+1)*1i*a*Q0.*wZt{1}.*Z_r_s{3,1}./s3_r_s./s1;
                media.p_inc(indices)  = temp(indices);
                temp = (2*n+1)*1i*a*Q0.*Z_r_s{1,1}.*wZt{3}./s1_r_s./s3;
                media.p_inc(~indices) = temp(~indices);
            otherwise
                media.p_inc = zeros(size(zeta),class(zeta));
        end
    end
    if layer{m}.calc_dp_incdr && options.p_inc_fromSeries
        switch applyLoad
            case 'planeWave'
                s1 = exp(exponent_(1,nu,zeta,nu_a));
                media.dp_incdr = (2*n+1)*1i^n*a*Q0.*dwZt{1}./s1;
            case 'pointCharge'
                s1 = exp(exponent_(1,nu,zeta,nu_a));
                s3 = exp(exponent_(3,nu,zeta,nu_a));
                Z_r_s = layer{m}.Z_r_s;
                r_s = options.r_s;
                s1_r_s = exp(exponent_(1,nu,a*r_s,nu_a));
                s3_r_s = exp(exponent_(3,nu,a*r_s,nu_a));
                media.dp_incdr = zeros(size(zeta),class(zeta));
                indices = repmat(r < r_s,numel(a),1);
                temp = (2*n+1)*1i*a*Q0.*dwZt{1}.*Z_r_s{3,1}./s3_r_s./s1;
                media.dp_incdr(indices)  = temp(indices);
                temp = (2*n+1)*1i*a*Q0.*Z_r_s{1,1}.*dwZt{3}./s1_r_s./s3;
                media.dp_incdr(~indices) = temp(~indices);
            otherwise
                media.dp_incdr = zeros(size(zeta),class(zeta));
        end
    end
    if layer{m}.calc_dp_incdt && options.p_inc_fromSeries
        switch applyLoad
            case 'planeWave'
                s1 = exp(exponent_(1,nu,zeta,nu_a));
                media.dp_incdt = (2*n+1)*1i^n*Q1.*wZt{1}./s1;
            case 'pointCharge'
                s1 = exp(exponent_(1,nu,zeta,nu_a));
                s3 = exp(exponent_(3,nu,zeta,nu_a));
                Z_r_s = layer{m}.Z_r_s;
                r_s = options.r_s;
                s1_r_s = exp(exponent_(1,nu,a*r_s,nu_a));
                s3_r_s = exp(exponent_(3,nu,a*r_s,nu_a));
                media.dp_incdt = zeros(size(zeta),class(zeta));
                indices = repmat(r < r_s,numel(a),1);
                temp = (2*n+1)*1i*a*Q1.*wZt{1}.*Z_r_s{3,1}./s3_r_s./s1;
                media.dp_incdt(indices)  = temp(indices);
                temp = (2*n+1)*1i*a*Q1.*Z_r_s{1,1}.*wZt{3}./s1_r_s./s3;
                media.dp_incdt(~indices) = temp(~indices);
            otherwise
                media.dp_incdt = zeros(size(zeta),class(zeta));
        end
    end
    for i = 1:numel(besselIndices)
        if besselIndices(i)
            w{i} = w_(n,i,xit_m,zeta,nu_a,exponentShift);
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

A_p = rho.*omega.^2.*A;
if isOuterDomain
    if layer{m}.calc_p_0
        cs = exp(exponent_(3,nu,xit_m,nu_a)-exponentShift);
        if isinf(cs)
            warning('e3Dss:infWeight','A weight evaluation was too large')
        end
        h_n_0  = 1i^(-n-1)./a;
        media.p_0 = A_p*Q0.*h_n_0.*cs;
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
        media.p = A_p*Q0.*wZt{3};
    end
    if layer{m}.calc_dpdr
        media.dpdr = A_p.*a*Q0.*dwZt{3};
    end
    if layer{m}.calc_dpdt % calculate a scaled dpdt (scaled by csc)
        media.dpdt = A_p*Q1.*wZt{3};
    end
    if layer{m}.calc_d2pdr2
        media.d2pdr2 = A_p.*a.^2*Q0.*d2wZt{3};
    end
    if layer{m}.calc_d2pdt2
        media.d2pdt2 = A_p*Q2.*wZt{3};
    end
elseif isSphere
    if layer{m}.calc_p
        media.p = A_p*Q0.*wZt{1};
        if n == 0
            media.p(:,indices) = repmat(A_p,1,sum(indices)).*sOrigin_xi;
        else
            media.p(:,indices) = 0;
        end
    end
    if layer{m}.calc_dpdr
        media.dpdr = A_p.*a*Q0.*dwZt{1};
        if n == 1
            media.dpdr(:,indices) = repmat(a/3.*A_p,1,sum(indices)).*sOrigin_xi;
        else
            media.dpdr(:,indices) = 0;
        end
    end
    if layer{m}.calc_dpdt % calculate a scaled dpdt (scaled by csc)
        media.dpdt = A_p*Q1.*wZt{1};
        media.dpdt(:,indices) = 0;
    end
    if layer{m}.calc_d2pdr2
        media.d2pdr2 = A_p.*a.^2*Q0.*d2wZt{1};
        if n == 0
            media.d2pdr2(:,indices) = repmat(-a.^2.*A_p,1,sum(indices)).*sOrigin_xi;
        else
            media.d2pdr2(:,indices) = 0;
        end
    end
    if layer{m}.calc_d2pdt2
        media.d2pdt2 = A_p*Q2.*wZt{1};
        media.d2pdt2(:,indices) = 0;
    end
else
    if layer{m}.calc_p
        media.p = A_p(:,1)*Q0.*wZt{1} + A_p(:,2)*Q0.*wZt{2};
    end
    if layer{m}.calc_dpdr
        media.dpdr = A_p(:,1).*a*Q0.*dwZt{1} + A_p(:,2).*a*Q0.*dwZt{2};
    end
    if layer{m}.calc_dpdt % calculate a scaled dpdt (scaled by csc)
        media.dpdt = A_p(:,1)*Q1.*wZt{1} + A_p(:,2)*Q1.*wZt{2};
    end
    if layer{m}.calc_d2pdr2
        media.d2pdr2 = A_p(:,1).*a.^2*Q0.*d2wZt{1} + A_p(:,2).*a.^2*Q0.*d2wZt{2};
    end
    if layer{m}.calc_d2pdt2
        media.d2pdt2 = A_p(:,1)*Q2.*wZt{1} + A_p(:,2)*Q2.*wZt{2};
    end
end
r2 = r.^2;
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
    media.u_r = u_r;
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
    media.u_t = u_t;
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
    media.du_rdr = du_rdr;
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
    media.du_rdt = du_rdt;
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
    media.du_tdr = du_tdr;
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
    media.du_tdt = du_tdt;
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
    
	media.sigma_rr = sigma_rr;
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
    
	media.sigma_tt = sigma_tt;
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
    
	media.sigma_pp = sigma_pp;
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
    
	media.sigma_rt = sigma_rt;
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
    media.navier1 = dsigma_rr_dr + dsigma_rt_dt + 1./R.*(2*sigma_rr-sigma_tt-sigma_pp+sigma_rt.*cos(Theta));
    media.navier2 = dsigma_rt_dr + dsigma_diffr + 3./R.*sigma_rt.*sin(Theta);
    if isSphere
        if n == 1
            media.navier1(:,indices) = repmat(-omega.^2.*rho.*(a.*A(:,1).*sOrigin_xi - 2*b.*B(:,1).*sOrigin_eta)/3,1,sum(indices));
        else
            media.navier1(:,indices) = 0;
        end
        media.navier2(:,indices) = 0;
    end
end