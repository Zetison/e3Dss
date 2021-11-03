function layer = setSkelton1997taoParameters()

if 1
    R = 1;
    t_steel = 0.02;
    t_coating = 0.02;

    layer{1}.media = 'fluid';
    layer{1}.R_i = R+t_coating;
    layer{1}.rho = 1000;
    layer{1}.c_f = 1500;

    rho_c = 800;    % density of coating
    eta = 0.1;  	% loss factor
    c_l = 123*sqrt(1-1i*eta);
    c_s = 33*sqrt(1-1i*eta);
    E = c_s^2*rho_c*(3*c_l^2-4*c_s^2)/(c_l^2-c_s^2);
    nu = (c_l^2-2*c_s^2)/(2*(c_l^2-c_s^2));
    layer{2}.media = 'solid';
    layer{2}.R_i = R;
    layer{2}.rho = 800;
    layer{2}.E = E; % 0.260e7
    layer{2}.nu = nu; % 0.460

    rho_s = 7700;   % density of steel
    eta = 0.01;  	% loss factor
    layer{3}.media = 'solid';
    layer{3}.R_i = R-t_steel;
    layer{3}.rho = rho_s;
    layer{3}.E = 0.195e12*sqrt(1-1i*eta);
    layer{3}.nu = 0.290*sqrt(1-1i*eta);
else
    R = 1;
    t_steel = 0.02;
    t_coating = 0.02;

    layer{1}.media = 'fluid';
    layer{1}.R_i = R+t_coating;
    layer{1}.rho = 1000;
    layer{1}.c_f = 1500;

    rho_c = 800;    % density of coating
    c_l = 123;
    c_s = 33;
    E = c_s^2*rho_c*(3*c_l^2-4*c_s^2)/(c_l^2-c_s^2);
    nu = (c_l^2-2*c_s^2)/(2*(c_l^2-c_s^2));
    layer{2}.media = 'solid';
    layer{2}.R_i = R;
    layer{2}.rho = 800;
    layer{2}.E = E; % 0.260e7
    layer{2}.nu = nu; % 0.460

    rho_s = 7700;   % density of steel
    layer{3}.media = 'solid';
    layer{3}.R_i = R-t_steel;
    layer{3}.rho = rho_s;
    layer{3}.E = 0.195e12;
    layer{3}.nu = 0.290;
end