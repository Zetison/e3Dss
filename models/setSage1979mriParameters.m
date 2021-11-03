function layer = setSage1979mriParameters()
a = 1;

layer{1}.media = 'fluid';
layer{1}.R_i = a;
layer{1}.rho = 1025;
layer{1}.c_f = 1531;
layer{1}.calc_p_0 = true; % Calculate the far field pattern

layer{2}.media = 'fluid';
layer{2}.R_i = 0;
layer{2}.rho = 1.293;
layer{2}.c_f = 346.2;