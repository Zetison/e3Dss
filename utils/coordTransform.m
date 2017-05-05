function [r_m, theta_m, phi_m, A] = coordTransform(v, d_vec)
% This function first orthogonally transform the Cartesian coordinates in v
% to a cartesian coordinate system with k_vec/norm(k_vec) as the third unit
% vector. Finally, the spherical coordinates are computed in this new
% coordinate system.

if isrow(d_vec)
    d_vec = d_vec.';
end
n_vec = d_vec/norm(d_vec);

%% Create orthonormal basis (e_x_m, e_y_m, n_vec)
e1 = [1; 0; 0];

if isa(v,'sym')
    e1 = vpa(e1);
end
e_x_m = e1 - dot(e1,n_vec)*n_vec;
if sum(abs(e_x_m)) < 1e-15
    e_x_m = [0; 1; 0];
    if isa(v,'sym')
        e_x_m = vpa(e_x_m);
    end
end
e_x_m = e_x_m/norm(e_x_m);

e_y_m = cross(n_vec,e_x_m);

%% Create transformation matrix
A = [e_x_m, e_y_m, n_vec];

%% Do transformations
v_m = (A\v.').';

x_m = v_m(:,1);
y_m = v_m(:,2);
z_m = v_m(:,3);

%% Transform Cartesian coordinates to spherical coordinates
r_m = sqrt(x_m.^2 + y_m.^2 + z_m.^2);
theta_m = acos(z_m./r_m);
theta_m(r_m < eps) = 0;
phi_m = atan2(y_m, x_m);