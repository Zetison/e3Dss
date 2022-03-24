function layer = setAyres1987arsParameters(fluid,alpha,beta)
if nargin < 1
    fluid = 'air';
end
if nargin < 3
    alpha = 0;
    beta = 0;
end
layer{1}.media = 'fluid';
layer{1}.R = 1;
if strcmp(fluid,'water')
    layer{1}.rho = 1.0*1e3;
    layer{1}.c = 1.493*1e5/100;
elseif strcmp(fluid,'air')
    layer{1}.rho = 0.0012*1e3;
    layer{1}.c = 0.334*1e5/100;
end

rho = 1.13*1e3;    % density of rubber
c_l = 1.4*1e5/100;
c_s = 0.94*1e4/100;

c_l = c_l*sqrt(1-1i*(alpha+2*beta)/(rho*c_l^2));
c_s = c_s*sqrt(1-1i*beta/(rho*c_s^2));
E = c_s^2*rho*(3*c_l^2-4*c_s^2)/(c_l^2-c_s^2);
nu = (c_l^2-2*c_s^2)/(2*(c_l^2-c_s^2));
layer{2}.E = E;
layer{2}.nu = nu;
layer{2}.lossFactor = 0;
layer{2}.media = 'solid';
layer{2}.R = 0;
layer{2}.rho = rho;