function exponent = exponent_(i,nu,x,nu_a)

exponent = zeros(size(x),class(x));
if nu_a ~= -1
    if i == 1
        exponent = nu*zeta23_(x/nu);
    elseif i == 2
        exponent = -abs(real(nu*zeta23_(x/nu)));
    else
        zeta23 = zeta23_(x/nu);
        exponent = -nu.*zeta23;
        indices = imag(x) < 0;
        exponent(indices) = -abs(real(nu.*zeta23(indices)));
    end
end