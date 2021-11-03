function indices = indices_(nu,z,nu_a)
% indices = and(nu > (z+32)/0.97, nu > nu_a);  % nu > (z+32)/0.97
indices = and(nu > z, and(abs(nu.*zeta23_(z./nu)) > log(sqrt(realmax(class(z)))), nu > nu_a));
% indices = and(nu > z, and(abs(nu*zeta23_(z/nu)) > log(realmax(class(z))*eps), nu > nu_a));

