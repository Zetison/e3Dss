function layer = setS45Parameters()

R = 5;
R2 = 4;
t_steel2 = 0.02;
t_steel = 0.008;
t_coating = 0.01;
rho_water = 1000;
c_f_water = 1500;

layer{1}.media = 'fluid';
layer{1}.R = R+t_coating;
layer{1}.rho = rho_water;
layer{1}.c_f = c_f_water;

rho = 800;    % density of coating
layer{2}.lossFactor = 0.1;
layer{2}.media = 'solid';
layer{2}.R = R;
layer{2}.rho = rho;
if 1
    layer{2}.E = 0.260e7;
    layer{2}.nu = 0.460;
else
    c_l = 123; %#ok
    c_s = 33;
    E = c_s^2*rho*(3*c_l^2-4*c_s^2)/(c_l^2-c_s^2);
    nu = (c_l^2-2*c_s^2)/(2*(c_l^2-c_s^2));
    layer{2}.E = E;
    layer{2}.nu = nu;
end

rho = 7850;    % density of steel
layer{3}.media = 'solid';
layer{3}.lossFactor = 0.001;
layer{3}.R = R-t_steel;
layer{3}.rho = rho;  % density of steel
layer{3}.E = 210e9;
layer{3}.nu = 0.3;
if false
    layer{2}.lossFactor = 0; %#ok
    layer{3}.lossFactor = 0;
end

layer{4} = layer{1};
layer{4}.R = R2;

layer{5} = layer{3};
layer{5}.R = R2-t_steel2;

layer{6} = layer{1};
layer{6}.R = 0;
layer{6}.rho = 1.2;
layer{6}.c_f = 340;
