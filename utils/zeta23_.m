function zeta = zeta23_(z)
% zeta = zeros(size(z),class(z));
% indices = real(z) <= 1;
% z1 = z(indices);
% z2 = z(~indices);
% sqrt1mz1 = sqrt(1-z1.^2);
% zeta(indices) = log((1+sqrt1mz1)./z1) - sqrt1mz1;
% zeta(~indices) = -1i*(sqrt(z2.^2-1) - asec(z2));

sqrt1mz1 = sqrt(1-z.^2);
zeta = log((1+sqrt1mz1)./z) - sqrt1mz1;
% zeta = log((1+sqrt(1-z.^2))./z) - sqrt(1-z.^2);
% zeta = (sqrt(z.^2-1) - asec(z))/(-1i);
end