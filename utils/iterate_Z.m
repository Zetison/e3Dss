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