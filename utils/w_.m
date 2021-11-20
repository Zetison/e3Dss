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
% if any(size(x) ~= size(y))
%     xx = repmat(x,1+size(y)-size(x));
% else
%     xx = x;
% end
% indices_x = indices_(i,nu,xx,nu_a);
% indices_y = indices_(i,nu,y,nu_a);
% w = zeros(size(xx));
% if i == 1
%     w(indices_x) = nu*zeta23_(xx(indices_x)/nu);
%     w(indices_y) = w(indices_y) - nu*zeta23_(y(indices_y)/nu);
% elseif i == 2
%     w(indices_x) = -abs(real(nu*zeta23_(xx(indices_x)/nu)));
%     w(indices_y) = w(indices_y) + abs(real(nu*zeta23_(y(indices_y)/nu)));
% else
%     w(indices_x) = -nu*zeta23_(xx(indices_x)/nu);
%     w(indices_y) = w(indices_y) + nu*zeta23_(y(indices_y)/nu);
% end

w = exp(exponent_(i,nu,x,nu_a)-exponent_(i,nu,y,nu_a)-exponentShift);
if any(isinf(w))
	warning('e3Dss:infWeight','A weight evaluation was too large')
end
