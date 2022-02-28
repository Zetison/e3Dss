function zeta = zeta_(z,zeta23)
prec = class(z);
twoThirds = getC(prec,'2/3');
if nargin < 2
    zeta = (3/2*zeta23_(z)).^twoThirds;
else
    zeta = (3/2*zeta23).^twoThirds;
end
    
end