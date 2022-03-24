close all
clear all %#ok

startup

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The default example
% When calling the e3Dss function without any extra options, the default
% set of parameters will be chosen. In particular, the rigid scattering of
% a unit sphere impinged by a plane wave traveling along the z-axis is
% simulated at f=1kHz using c = 1500m/s as the speed of sound. 
% The result is plotted in the far field as a function of the polar angle.

theta_arr = linspace(0,pi,2000)'; % Set of angles used for plotting
X = [sin(theta_arr), zeros(size(theta_arr)), cos(theta_arr)]; % Evaluate physical location of plotting points


% General parameters in layer i
layer{1}.media 	= 'fluid'; % Media; % solid or fluid (Navier equation or Helmholtz equation)
layer{1}.X     	= X;       % Evaluation points
layer{1}.calc_p_0 = true; % Calculate the far field pattern

layer = e3Dss(layer); % Compute solution

% Plot real part of the scattered pressure at the surface of the unit sphere
plot(theta_arr, real(layer{1}.p_0))
xlim([0,pi])
legend('Modulus of scattered field')
xlabel('$$\theta$$, polar angle','interpreter','latex')
ylabel('$$\mathrm{real}(p_0)$$, real part of scattered pressure','interpreter','latex')