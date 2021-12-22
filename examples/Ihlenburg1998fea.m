%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This study correspond to Figure9a, Figure 9b, Figure 10, Figure 11 and Figure 12 in Venas2019e3s
% Venas2019e3s is available at https://doi.org/10.1016/j.cma.2018.02.015 (open access version at http://hdl.handle.net/11250/2493754)
% It is based on the example in Ihlenburg1998fea Figure 5.2
% Ihlenburg1998fea is available at https://books.google.no/books?id=JrMPBwAAQBAJ

close all
clear all %#ok

startup
resultsFolder = [folderName '/Ihlenburg1998fea'];
if ~exist(resultsFolder, 'dir')
    mkdir(resultsFolder);
end

%% Ihlenburg (1998) example
Ieye = [1, 0;
        0, 1;
        0, 0];
ESBC = 0;
for i = 1:3
    if i == 1
        nFreqs = 2000;
        color = [0,70,147]/255;
        legendEntry = 'Sound-hard boundary condition';
    elseif i == 2
        nFreqs = 5000;
        color = [178,0,0]/255;
        legendEntry = 'Sound-soft boundary condition';
    else
        nFreqs = 5000;
        color = [59,124,37]/255;
        legendEntry = 'Neumann-Neumann boundary condition';
    end
    layer = setIhlenburgParameters();
    layer{1}.calc_p_0 = true; % Calculate the far field pattern
    
    SHBC = Ieye(i,1);
    SSBC = Ieye(i,2);
    defineBCstring

    k = linspace(2/nFreqs,2,nFreqs)'; % wave number
    omega = k*layer{1}.c_f;   % Wave number for outer fluid domain

    theta = 180*pi/180;
    d_vec = [1;0;0];
    options = struct('BC', BC,...
                     'd_vec', d_vec, ...
                     'omega', omega, ...
                     'P_inc', 1);
    R = layer{1}.R;
    if SHBC
        specialValues = [];
    else
        if 0
            layer{1}.X = R_o(1)*[cos(0),0,0];
            f = @(k)-objFunc(k,layer,options);
            specialValues = findExtremas(f, 2/nFreqs, 2, 100000)';
            layer{1}.X = R_o(1)*[cos(pi),0,0];
            f = @(k)-objFunc(k,layer,options);
            specialValues = [specialValues; findExtremas(f, 2/nFreqs, 2, 100000)'];
            delta = 1e-5*k(end);
            specialValues = sort([specialValues; (specialValues-delta); (specialValues+delta)]);
            save(['../miscellaneous/Ihlenburg_' BC '_extremas'], 'specialValues')
        else
            load(['../miscellaneous/Ihlenburg_' BC '_extremas'])
        end
        delta = 1e-4;
        specialValues = sort([specialValues; (specialValues-delta); (specialValues+delta)]);
    end
    k = unique(sort([k; specialValues]));
    omega = k*layer{1}.c_f;   % Wave number for outer fluid domain
    options.omega = omega;

    layer{1}.X = R*[cos(pi),0,0;
                      cos(0),0,0];
    layer = e3Dss(layer, options);

    figure(3)
    F = layer{1}.p_0;
    TS = 20*log10(abs(F));
    plot(k*R, TS(1,:),'DisplayName',legendEntry,'color',color)
    set(0,'defaulttextinterpreter','latex')
    hold on
    title('Ihlenburg (1998) example, $$\theta = 180^\circ$$')
    xlabel('$$k_1 R_{0,1}$$')
    xlim([0, max(k*R)])
    ylim([-50, 35])
    ylabel('TS [dB]')  
    legend('off');
    legend('show','location','southeast');
    savefig([resultsFolder '/Figure9a'])

    figure(4)
    F = layer{1}.p_0;
    TS = 20*log10(abs(F));
    plot(k*R, TS(2,:),'DisplayName',legendEntry,'color',color)
    set(0,'defaulttextinterpreter','latex')
    hold on
    title('Ihlenburg (1998) example, $$\theta = 0^\circ$$')
    xlabel('$$k_1 R_{0,1}$$')
    xlim([0, max(k*R)])
    ylim([-50, 35])
    ylabel('TS [dB]')  
    legend('off');
    legend('show','location','southeast');
    savefig([resultsFolder '/Figure9b'])
    
    
    if 1
        figure(40+i)
        nFreqs = 500;
        k = linspace(2/nFreqs,2,nFreqs)'; % wave number
        k = unique(sort([k; specialValues]));
        omega = k*layer{1}.c_f;   % Wave number for outer fluid domain
        options.omega = omega;

        createConvergencePlot('3D',layer,options,35, [])
        savefig([resultsFolder '/Figure' num2str(9+i)])
    end
end


function TS = objFunc(k,layer,options)

options.Display = 'none';
options.omega = k*layer{1}.c_f;
layer = e3Dss(layer, options);
TS = 20*log10(abs(layer{1}.p_0));

end
