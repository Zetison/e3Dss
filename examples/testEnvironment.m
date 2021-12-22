close all
clear all %#ok

startup

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The default example
% When calling the e3Dss function without any extra options, the default
% set of parameters will be chosen. In particular, the rigid scattering of
% a unit sphere impinged by a plane wave traveling along the z-axis is
% simulated at f=1kHz using c_f = 1500m/s as the speed of sound. 
% The result is plotted in the far field as a function of the polar angle.

alpha = 0;
beta = 0;
r = 1;
v = zeros(length(alpha)*length(beta)*length(r),3);
counter = 1;
for l = 1:length(r)
    for j = 1:length(beta)
        for i = 1:length(alpha)
            v(counter,:) = r(l)*([cos(beta(j))*cos(alpha(i)), cos(beta(j))*sin(alpha(i)), sin(beta(j))]);
            counter = counter + 1;
        end
    end
end
nu_a = 100;
% nu_a = -1;
options = struct('BC', 'SHBC', ...
               'd_vec', [1;0;0], ...
               'N_max', Inf, ...
               'P_inc', 1, ...
               'omega', 1.524, ...
           'applyLoad', 'planeWave', ...
             'Display', 'none', ...
             'nu_a', nu_a, ...
    'p_inc_fromSeries', 0);
% General parameters in layer i
layer{1}.R = 5.075;
layer{1}.rho = 1000;
layer{1}.c_f = 1524;

layer{1}.media 	= 'fluid'; % Media; % solid or fluid (Navier equation or Helmholtz equation)
layer{1}.X     	= layer{1}.R*[-1,0,0];       % Evaluation points
% layer{1}.X     	= [-1,0,0];       % Evaluation points
layer{1}.calc_p_0 = true; % Calculate the far field pattern

layer = e3Dss(layer,options); % Compute solution

% Plot real part of the scattered pressure at the surface of the unit sphere
plot(alpha, real(layer{1}.p_0))
xlim([0,pi])
legend('Modulus of scattered field')
xlabel('$$\theta$$, polar angle','interpreter','latex')
ylabel('$$\mathrm{real}(p_0)$$, real part of scattered pressure','interpreter','latex')


layer{1}.p_0
layer{1}.X = [-1,0,0];       % Evaluation points
layer = e3Dss(layer,options); % Compute solution
layer{1}.p_0


