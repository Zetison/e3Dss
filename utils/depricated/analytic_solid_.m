function u = analytic_solid_(v,options)
error('Depricated')
d_vec = options.d_vec;
if size(d_vec,2) > 1
    nPts = size(v,1);
    v = v(1,:); % Assume monostatic scattering
    options.d_vec = d_vec(:,1);
    data = e3Dss({[], v},options);
    u = [data(1).u_x, data(1).u_y, data(1).u_z];
    u = ones(nPts,1)*u;
else
    data = e3Dss({[], v},options);
    u = [data(1).u_x, data(1).u_y, data(1).u_z];
end
