function layer = setS135Parameters(prec)

if nargin < 1
    prec = 'double';
end

switch prec
    case 'double'
        rho_s = [7850, 7850, 7850]; % Density of solid
        rho_f = [1000, 1000, 1.2, 1.2]; % Density of fluids
        c = [1500, 1500, 340, 340];  % Speed of sound in fluid domains
        t = [0.008, 0.02, 0.05]; % The thickness of the sphere
        E = [210e9, 210e9, 210e9]; % Youngs modulus of elastic material
        nu = [0.3, 0.3, 0.3]; % Poisson ratio of elastic material
        R_o = [5, 3, 1, 0];
    case 'sym'
        rho_s = vpa(str2sym('[7850, 7850, 7850]')); % Density of solid
        rho_f = vpa(str2sym('[1000, 1000, 1.2, 1.2]')); % Density of fluids
        c = vpa(str2sym('[1500, 1500, 340, 340]'));  % Speed of sound in fluid domains
        t = vpa(str2sym('[0.008, 0.02, 0.05]')); % The thickness of the sphere
        E = vpa(str2sym('[210e9, 210e9, 210e9]')); % Youngs modulus of elastic material
        nu = vpa(str2sym('[0.3, 0.3, 0.3]')); % Poisson ratio of elastic material
        R_o = vpa(str2sym('[5, 3, 1, 0]'));
    case 'mp'
        rho_s = mp('[7850, 7850, 7850]'); % Density of solid
        rho_f = mp('[1000, 1000, 1.2, 1.2]'); % Density of fluids
        c = mp('[1500, 1500, 340, 340]');  % Speed of sound in fluid domains
        t = mp('[0.008, 0.02, 0.05]'); % The thickness of the sphere
        E = mp('[210e9, 210e9, 210e9]'); % Youngs modulus of elastic material
        nu = mp('[0.3, 0.3, 0.3]'); % Poisson ratio of elastic material
        R_o = mp('[5, 3, 1, 0]');
end
convertToLayerFormat
