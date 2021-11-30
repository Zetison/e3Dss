function w = w_(n,i,x,y,nu_a,exponentShift)

if isempty(x)
    w = [];
    return
end
if isempty(y)
    w = [];
    return
end
nu = n + 1/2;

w = exp(exponent_(i,nu,x,nu_a)-exponent_(i,nu,y,nu_a)-exponentShift);
if any(isinf(w))
	warning('e3Dss:infWeight','A weight evaluation was too large')
end
