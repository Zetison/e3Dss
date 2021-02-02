function layer = setHetmaniukParameters(noDomains)
if nargin < 1
    noDomains = 2;
end
rho_f = [1000,1000]; % Density of outer fluid
rho_s = 7850; % Density of solid
c_f = [1500,1500];   % Speed of sound in outer (exterior) fluid domain
t = 0.05; % The thickness of the sphere
R_o = [1,0]; % Outer radius of shell
E = 2.0e11; % Youngs modulus of elastic material
nu = 0.3; % Poisson ratio of elastic material

convertToLayerFormat
layer = layer(1:noDomains);
