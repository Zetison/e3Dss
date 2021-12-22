function layer = addCoating(layer,t,coatingLayer,i)
if nargin < 4
    i = 1;
end
if nargin < 3
    % Parameters from Skelton1997tao
    rho_c = 800;    % density of coating
    coatingLayer.lossFactor = 0.1;
    coatingLayer.media = 'solid';
    coatingLayer.R = layer{1}.R;
    coatingLayer.rho = rho_c;
    coatingLayer.E = 0.260e7;
    coatingLayer.nu = 0.460;
end
layer = [layer(i); {coatingLayer}; layer(i+1:end)];
layer{i}.R = layer{i}.R+t;