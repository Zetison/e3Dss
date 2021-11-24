function PI = getPI(prec)
switch prec
    case 'single'
        PI = pi;
    case 'double'
        PI = pi;
    case 'sym'
        PI = vpa(pi,digits);
    case 'mp'
        PI = mp('pi');
end