function g = g_(n,i,x,nu_a)

nu = n + 1/2;
indices1 = indices_(nu,x,nu_a);
indices2 = indices_(nu+1,x,nu_a);
g = zeros(size(x));
g(indices1) = (-1)^(i+1)*nu*zeta23_(x(indices1)/nu);
g(indices2) = g(indices2) - (-1)^(i+1)*(nu+1)*zeta23_(x(indices2)/(nu+1));

g = exp(g);
if any(isinf(g))
	warning('e3Dss:infWeight','A weight evaluation was too large')
end