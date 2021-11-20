function zeta = zeta_(z,zeta23)
if nargin < 2
    zeta = (3/2*zeta23_(z)).^(2/3);
else
    zeta = (3/2*zeta23).^(2/3);
end
    
end