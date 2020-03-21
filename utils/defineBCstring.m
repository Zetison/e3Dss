
if SHBC
    layer = layer(1:end-2);
    BC = 'SHBC';
elseif ESBC
    layer = layer(1:end-1);
    layer{end}.R_i = 0;
    BC = 'ESBC';
elseif SSBC
    layer = layer(1:end-1);
    BC = 'SSBC';
else
    BC = 'NNBC';
end
M = length(layer);