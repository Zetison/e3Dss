function layer = setS35Parameters(prec)

if nargin < 1
    prec = 'double';
end

switch prec
    case 'double'
        rho_s = [7850, 7850]; % Density of solid
        rho_f = [1000, 1000, 1.2]; % Density of fluids
        c_f = [1500, 1500, 340];  % Speed of sound in fluid domains
        t = [0.008, 0.02]; % The thickness of the sphere
        E = [210e9, 210e9]; % Youngs modulus of elastic material
        nu = [0.3, 0.3]; % Poisson ratio of elastic material
        R_o = [5, 3, 0]; 
    case 'sym'
        rho_s = vpa('[7850, 7850]'); % Density of solid
        rho_f = vpa('[1000, 1000, 1.2]'); % Density of fluids
        c_f = vpa('[1500, 1500, 340]');  % Speed of sound in fluid domains
        t = vpa('[0.008, 0.02]'); % The thickness of the sphere
        E = vpa('[210e9, 210e9]'); % Youngs modulus of elastic material
        nu = vpa('[0.3, 0.3]'); % Poisson ratio of elastic material
        R_o = vpa('[5, 3, 0]'); 
end
convertToLayerFormat
