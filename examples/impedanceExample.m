close all
clear all %#ok

startup

k10start = -1;
k10end = 1;
c = 1500;
k = 10.^linspace(k10start,k10end,1000);

omega =.c*k; % Angular frequency
layer{1} = struct('media', 'fluid', .c',.c);

layer{1}.X = [0,0,-1];
layer{1}.calc_p_0 = 1; % Calculate the far field pattern

options = struct('BC', 'IBC', ...
                 'omega', omega);
             

figure(1)
hold on
for z = 10.^linspace(4,18,15)
    options.z = 1i*z.*k;
    layer = e3Dss(layer, options);
    TS = 20*log10(abs(layer{1}.p_0));
%     TS(54)
    legendStr = sprintf('IBC, z = %g', z);
    plot(k, TS, 'DisplayName', legendStr)
end
          
options.BC = 'SSBC';   
layer = e3Dss(layer, options);
TS = 20*log10(abs(layer{1}.p_0));
plot(k, TS, '--', 'DisplayName', 'SSBC')

options.BC = 'SHBC';
layer = e3Dss(layer, options);
TS = 20*log10(abs(layer{1}.p_0));
plot(k, TS, '--', 'DisplayName', 'SHBC')

legend show

