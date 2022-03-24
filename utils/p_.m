function H = p_(n,k,R,Rt_m,Z,isSphere,isOuterDomain,nu_a,A,B)

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