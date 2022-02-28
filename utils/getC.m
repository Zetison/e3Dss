function C = getC(prec,constant)
% Get various constants in given prec(ision)

switch prec
    case 'single'
        C = single(str2num(constant)); %#ok: Since str2double does not handle fractions like 3/2
    case 'double'
        C = str2num(constant); %#ok: see comment above
    case 'sym'
        C = vpa(constant,digits);
    case 'mp'
        C = mp(constant);
end