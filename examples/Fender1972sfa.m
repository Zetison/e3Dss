%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This study correspond to Figure 13a, Figure 13b and Figure 14 in Venas2019e3s
% Venas2019e3s is available at https://doi.org/10.1016/j.cma.2018.02.015 (open access version at http://hdl.handle.net/11250/2493754)
% It is based on the example in Fender1972sfa Figure 2 and Figure 3
% Fender1972sfa is available at https://apps.dtic.mil/docs/citations/AD0752760

close all
clear all %#ok

startup
resultsFolder = [folderName '/Fender1972sfa'];
if ~exist(resultsFolder, 'dir')
    mkdir(resultsFolder);
end

%% Fender (1972) example
layer = setFenderParameters();
nFreqs = 3000;
% nFreqs = 2;
k = linspace(32/nFreqs,32,nFreqs)';
R_o = layer{1}.R;
c = layer{1}.c;
omega = k*c(1);

d_vec = -[0,0,1].';
options = struct('d_vec', d_vec, ...
                 'BC', 'NNBC', ...
                 'omega', omega);

if true
    SPL_Fender0 = importdata('models/Fender1972sfa/Fig2.csv');
    SPL_Fender180 = importdata('models/Fender1972sfa/Fig3.csv');
    theta = 0;

    layer{1}.calc_p = true;
    if 0
        layer{1}.X = [0,0,R_o*cos(0)];
        f = @(k)-objFunc(k,layer,options);
        specialValues = findExtremas(f, 2/nFreqs, 32, 100000)';
        layer{1}.X = [0,0,R_o*cos(pi)];
        f = @(k)-objFunc(k,layer,options);
        specialValues = [specialValues; findExtremas(f, 2/nFreqs, 32, 100000)'];
        delta = 1e-5*k(end);
        specialValues = sort([specialValues; (specialValues-delta); (specialValues+delta)]);
        save('miscellaneous/Fender_extremas', 'specialValues')
    else
        load('miscellaneous/Fender_extremas')
    end
    k = unique(sort([k; specialValues]));
    X = [0,0,R_o*cos(0);
         0,0,R_o*cos(pi)];
    layer{1}.X = X;
    options.omega = k*c(1);
    tic
    [layer,N_eps] = e3Dss(layer, options);
    toc
    p_inc = @(v) exp(1i*dot3(v,d_vec)*k.');
    p_tot = layer{1}.p + p_inc(X);
    SPL = 20*log10(abs(p_tot)); % sound pressure level
    figure(2)
    plot(k*R_o, SPL(1,:), SPL_Fender0(:,1), SPL_Fender0(:,2))
    set(0,'defaulttextinterpreter','latex')
    title('Predicted total pressure as a function of $$ka$$ at the surface of the shell, $$\theta = 0^\circ$$')
    xlabel('$$k_1 R_{0,1}$$')
    ylabel('Surface sound pressure level [dB]')
    ylim([-80 120])
    legend({'Present work', 'Reference Solution from Fender (1972)'})
    xlim(R_o*[k(1) k(end)])
    savefig([resultsFolder '/Figure13a'])

    %%%%%%%%
    figure(3)
    plot(k*R_o, SPL(2,:), SPL_Fender180(:,1), SPL_Fender180(:,2))
    set(0,'defaulttextinterpreter','latex')
    title('Predicted total pressure as a function of $$ka$$ at the surface of the shell, $$\theta = 180^\circ$$')
    xlabel('$$k_1 R_{0,1}$$')
    ylabel('Surface sound pressure level [dB]')
    ylim([-60 120])
    legend({'Present work', 'Reference Solution from Fender (1972)'})
    xlim(R_o*[k(1) k(end)])
    savefig([resultsFolder '/Figure13b'])
end

if 1
    nFreqs = 2000;
    k_max = 450;
    k = [linspace(1e-300,k_max/nFreqs,200), linspace((k_max+1)/nFreqs,k_max,nFreqs)];
    omega = k*c(1);   % Wave number for outer fluid domain

    x = k*R_o(1);

    E = layer{2}.E;
    nu = layer{2}.nu;
    rho_s = layer{2}.rho;
    R = layer{2}.R;
    R_o = layer{1}.R;
    c = layer{1}.c;

    K = E./(3*(1-2*nu));
    G = E./(2*(1+nu));
    c_s_1 = sqrt((3*K+4*G)./(3*rho_s));
    c_s_2 = sqrt(G./rho_s);

    Upsilon = min([R./c_s_1, R./c_s_2, R_o..c(1:end-1)]);

    NN = 0:600;
    N = zeros(size(x));
    for i = 1:length(x)
        z = omega(i)*Upsilon;
        temp = abs(sqrt(pi/2)./sqrt(z).*bessely(NN,z)) - 10^290;
        indices = find(temp > 0);
        N(i) = NN(indices(1));
    end
    figure(4)
    hold on
    set(0,'defaulttextinterpreter','latex')
    plot(x, N, 'DisplayName','$N=\lceil\upsilon\rceil$ as a function of $k_1 R_{0,1}$ for when $\mathrm{y}_\upsilon(\omega\Upsilon)=10^{290}$','color',[0,70,147]/255)

    k_max = 1000;
    k = linspace(k_max/nFreqs,k_max,nFreqs);
    omega = k*c(1);   % Wave number for outer fluid domain
    options.omega = omega;
    options.nu_a = -1;
    layer{1}.X = [0,0,R_o*cos(0);
                  0,0,R_o*cos(pi)];
    tic
    [layer,N_eps,flag] = e3Dss(layer, options);
    toc
    plot(k(~flag)*R_o, N_eps(~flag), 'DisplayName','$N=N_\varepsilon$ needed for convergence of $p_1^{(N_\varepsilon)}$ within machine epsilon precision','color',[178,0,0]/255)
    xlabel('$k_1R_{0,1}$')
    ylabel('$N$')
    leg1 = legend('show','Location','northwest');
    set(leg1,'Interpreter','latex');

    savefig([resultsFolder '/Figure14'])
end

function SPL = objFunc(k,layer,options)

options.Display = 'none';
options.omega = k*layer{1}.c;
layer = e3Dss(layer, options);
p_inc = @(v) exp(1i*dot3(v,options.d_vec)*k);
p_tot = layer{1}.p + p_inc(layer{1}.X);
SPL = 20*log10(abs(p_tot)); % sound pressure level

end

