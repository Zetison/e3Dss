function zeta = zeta23_(z)
zeta = zeros(size(z));
indices = z <= 1;
z1 = z(indices);
z2 = z(~indices);

zeta(indices) = log((1+sqrt(1-z1.^2))./z1) - sqrt(1-z1.^2);
zeta(~indices) = sqrt(z2.^2-1) - asec(z2);

% zeta = log((1+sqrt(1-z.^2))./z) - sqrt(1-z.^2);
end