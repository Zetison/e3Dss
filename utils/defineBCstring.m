function layer = defineBCstring(layer,BC)
if strcmp(BC,'SHBC')
    layer = layer(1:end-2);
elseif strcmp(BC,'ESBC')
    layer = layer(1:end-1);
    layer{end}.R = 0;
elseif strcmp(BC,'SSBC')
    layer = layer(1:end-1);
end