function g = g_(n,i,x,nu_a)

nu = n + 1/2;
% indices1 = indices_(i,nu,x,nu_a);
% indices2 = indices_(i,nu+1,x,nu_a);
% g = zeros(size(x),class(x));
% if i == 1 % Scale for besselj
%     g(indices1) = nu*zeta23_(x(indices1)/nu);
%     g(indices2) = g(indices2) - (nu+1)*zeta23_(x(indices2)/(nu+1));
% elseif i == 2 % Scale for bessely
%     g(indices1) = -abs(real(nu*zeta23_(x(indices1)/nu)));
%     g(indices2) = g(indices2) + abs(real((nu+1)*zeta23_(x(indices2)/(nu+1))));
% else % Scale for besselh
%     g(indices1) = -nu*zeta23_(x(indices1)/nu);
%     g(indices2) = g(indices2) + (nu+1)*zeta23_(x(indices2)/(nu+1));
% end

% g = exp(g);
g = exp(exponent_(i,nu,x,nu_a)-exponent_(i,nu+1,x,nu_a));
if any(isinf(g))
	warning('e3Dss:infWeight','A weight evaluation was too large')
end