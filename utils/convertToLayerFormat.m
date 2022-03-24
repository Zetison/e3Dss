
M = (numel(rho_s)+numel(rho_f));
layer = cell(M,1);
for m = 1:M
    if mod(m,2)
        md2 = (m+1)/2;
        layer{m}.media 	= 'fluid'; % Media; % solid or fluid (Navier equation or Helmholtz equation)
        layer{m}.R = R_o(md2);
        layer{m}.c = c(md2);
        layer{m}.rho = rho_f(md2);
    else
        md2 = m/2;
        layer{m}.media 	= 'solid'; % Media; % solid or fluid (Navier equation or Helmholtz equation)
        layer{m}.R = R_o(md2)-t(md2);
        layer{m}.E = E(md2);
        layer{m}.nu = nu(md2);
        layer{m}.rho = rho_s(md2);
    end
end