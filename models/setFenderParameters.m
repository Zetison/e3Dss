function layer = setFenderParameters()

rho_f = [1.026e3, 1.21]; % Density of fluids
rho_s = 2.7e3; % Density of solid
c = [1.5e3, 343]; % Speed of sound in fluid domains
c_s_L = 6.412e3;
c_s_T = 3.043e3;
c_s_1 = c_s_L;
c_s_2 = c_s_T;
E = c_s_2^2*rho_s*(3*c_s_1^2-4*c_s_2^2)/(c_s_1^2-c_s_2^2);
nu = (c_s_1^2-2*c_s_2^2)/(2*(c_s_1^2-c_s_2^2));
R_o = [1, 0];
t = 0.05*R_o(1);

convertToLayerFormat
