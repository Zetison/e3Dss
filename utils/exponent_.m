function exponent = exponent_(i,nu,x,nu_a)

indices = indices_(i,nu,x,nu_a);
exponent = zeros(size(x),class(x));
if nu_a ~= -1
    if i == 1
        exponent(~indices) = -abs(imag(x(~indices)));
        exponent(indices) = nu*zeta23_(x(indices)/nu);
    elseif i == 2
        exponent(~indices) = -abs(imag(x(~indices)));
        exponent(indices) = -abs(real(nu*zeta23_(x(indices)/nu)));
    else
        exponent(~indices) = -1i*x(~indices);
        exponent(indices) = -nu*zeta23_(x(indices)/nu);
    end
end