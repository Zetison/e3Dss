function layer = setChangParameters()

a = 1; % mid surface radius
h = 0.01*a;
R_o = a + h/2;
R_i = a - h/2;

layer{1}.media 	= 'fluid'; % Media; % solid or fluid (Navier equation or Helmholtz equation)
layer{1}.R_i = R_o;
layer{1}.c_f = 1460;
layer{1}.rho = 1000;

layer{2}.media 	= 'solid'; % Media; % solid or fluid (Navier equation or Helmholtz equation)
layer{2}.R_i = R_i;
layer{2}.E = 2.0e11;
layer{2}.nu = 0.3;
layer{2}.rho = 7800;
