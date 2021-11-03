function w = w_(n,i,x,y,nu_a)

if isempty(x)
    w = zeros(size(x));
    return
end
if isempty(y)
    w = zeros(size(y));
    return
end
nu = n + 1/2;
if any(size(x) ~= size(y))
    xx = repmat(x,1+size(y)-size(x));
else
    xx = x;
end
indices_x = indices_(nu,xx,nu_a);
indices_y = indices_(nu,y,nu_a);
w = zeros(size(xx));
w(indices_x) = (-1)^(i+1)*nu*zeta23_(xx(indices_x)/nu);
w(indices_y) = w(indices_y) - (-1)^(i+1)*nu*zeta23_(y(indices_y)/nu);

w = exp(w);
if any(isinf(w))
	warning('e3Dss:infWeight','A weight evaluation was too large')
end
