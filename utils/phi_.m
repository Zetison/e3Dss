function phi = phi_(n,k,R,Rt_m,Z,isSphere,isOuterDomain,nu_a,A)

if nargin < 8
    zeta = k*R;
    zetat_m = k*Rt_m;
    if isSphere
        phi = zeros(1,1,length(k),class(R));
        phi(1,1,:) = w_(n,1,zetat_m,zeta,nu_a,0).*Z{1,1}; 
    elseif isOuterDomain
        phi = zeros(1,1,length(k),class(R));
        phi(1,1,:) = w_(n,3,zetat_m,zeta,nu_a,0).*Z{3,1}; 
    else
        phi = zeros(1,2,length(k),class(R));

        phi(1,1,:) = w_(n,1,zetat_m,zeta,nu_a,0).*Z{1,1};
        phi(1,2,:) = w_(n,2,zetat_m,zeta,nu_a,0).*Z{2,1};
    end
else
    if isOuterDomain
        phi = A*Q0.*wZt{3};
    elseif isSphere
        phi = A*Q0.*wZt{1};
        if n == 0
            phi(:,indices) = repmat(A,1,sum(indices)).*sOrigin_xi;
        else
            phi(:,indices) = 0;
        end
    else
        phi = A(:,1)*Q0.*wZt{1} + A(:,2)*Q0.*wZt{2};
    end
end