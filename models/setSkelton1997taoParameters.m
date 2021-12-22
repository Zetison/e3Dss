function layer = setSkelton1997taoParameters()

R = 1;
t_steel = 0.02;
t_coating = 0.02;

layer{1}.media = 'fluid';
layer{1}.R = R+t_coating;
layer{1}.rho = 1000;
layer{1}.c_f = 1500;

rho_c = 800;    % density of coating
layer{2}.lossFactor = 0.1;
layer{2}.media = 'solid';
layer{2}.R = R;
layer{2}.rho = rho_c;
if 1
    layer{2}.E = 0.260e7;
    layer{2}.nu = 0.460;
else
    c_l = 123; %#ok
    c_s = 33;
    E = c_s^2*rho_c*(3*c_l^2-4*c_s^2)/(c_l^2-c_s^2);
    nu = (c_l^2-2*c_s^2)/(2*(c_l^2-c_s^2));
    layer{2}.E = E;
    layer{2}.nu = nu;
end

rho_s = 7700;    % density of coating
layer{3}.media = 'solid';
layer{3}.lossFactor = 0.01;
layer{3}.R = R-t_steel;
layer{3}.rho = rho_s;  % density of steel
if 1
    layer{3}.E = 0.195e12;
    layer{3}.nu = 0.290;
else
    c_l = 5760; %#ok
    c_s = 3130;
    E = c_s^2*rho_s*(3*c_l^2-4*c_s^2)/(c_l^2-c_s^2);
    nu = (c_l^2-2*c_s^2)/(2*(c_l^2-c_s^2));
    layer{3}.E = E;
    layer{3}.nu = nu;
end
if false
    layer{2}.lossFactor = 0; %#ok
    layer{3}.lossFactor = 0;
end