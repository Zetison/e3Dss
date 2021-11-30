function P = legendre_(n, x)
% This function returns n'th Legendre polynomial evaluated at x

if iscolumn(x)
    x = x';
    xIsCol = true;
else
    xIsCol = false;
end
i = 1;
P = [zeros(size(x),class(x)); ones(size(x),class(x))];
if n < 0
    P = P(1,:);
    return
end
while i <= n
    P = legendreDerivs(i, x, P);
    i = i + 1;
end
P = P(2,:);
if xIsCol
    P = P.';
end