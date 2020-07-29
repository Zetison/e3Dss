function p = analytic_(v,layers,options)
error('Depricated, use e3Dss directly instead')
d_vec = options.d_vec;
if ~iscell(v)
    v = {v}; % assume points in outermost layer
end
if size(d_vec,2) > 1
    layers.X{1} = v{1}(1,:); % Assume monostatic scattering
    nPts = size(v,1);
    options.d_vec = d_vec(:,1);
    layers = e3Dss(layers,options);
    p = layers{1}.p;
    p = p*ones(nPts,1);
else
    layers.X = v;
    layers = e3Dss(v,options);
    p = layers{1}.p;
end
