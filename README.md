# Exact 3D scattering solutions
Exact 3D scattering solutions for spherical symmetric scatterers computes the solution to scattering problems on multilayered spherical (elastic or fluid) shells impinged by a plane wave or a wave due to a point source. The code is written in MATLAB and may for that reason be improved compuationally.

## Boundary conditions
The following boundary conditions can be used for the innermost layer:
- 'IBC' (Impedance Boundary Condition) simulates a impedance boundary condition (Robin boundary condition)
- 'SHBC' (Sound Hard Boundary Condition) simulates a rigid boundary
- 'SSBC' (Sound Soft Boundary Condition) simulates the innermost layer having zero pressure
- 'NNBC' (Neumann-Neumann boundary condition) simulates full acoustic structure interaction

## Loads (Incident wave type)
The following loads can be used
- 'planeWave' simulates a plane incident wave
- 'pointCharge' simulates an incident wave from a point charge (the location of the point source is given by options.d_vec*options.r_s)
- 'mechExcitation' simulates a point force mechanical excitation at options.d_vec*options.r_s with amplitude P_inc (options.r_s must be at layer{m}.R for a given layer m)
- 'surfExcitation' simulates a surface excitation over the region theta \in options.theta_s and phi \in [0, 2*pi] at r_s with amplitude P_inc (direction of z-axis given by options.d_vec)
- 'radialPulsation' simulates a spherical symmetric wave originating from infinity

## Parameters (with default values)
Global parameters and options
```Matlab
options.d_vec            = [0;0;1];      % Direction of the incident wave
options.omega            = 2*pi*1e3;     % Angular frequency (can be an array of angular frequencies)
options.P_inc            = 1;            % Amplitude of incident wave at the origin. P_inc can be given as a function handle P_inc(omega) where omega is the angular frequency
options.N_max            = Inf;          % Upper limit for the number of terms in the series
options.Display          = 'final';      % Print options ('final', 'iter' or 'none')
options.BC               = 'SHBC';       % Boundary condition on the inermost layer 'SSBC' (Sound soft boundary condition), 'NNBC' (Neumann-Neumann boundary condition) 
options.applyLoad        = 'planeWave';  % Incident wave type
options.r_s              = NaN;          % Radius to source location for point charge incident waves
options.p_inc_fromSeries = false;        % Calculate p_inc by series expansions (in terms of Bessel functions)
options.nu_a             = 100;          % Value of nu at which scaled asymptotic expansions are used in bessel_c (set nu_a = -1 to turn off scaling)
options.z                = 1;            % Impedance for an impedance boundary condition
options.Eps              = eps;          % Parameter for series truncation. The summation are terminated whenever the relative contribution of the given term is less then Eps. 
                                         % If vector input is given for either X or omega, the maximal relative contribution of the given term is compared with Eps
options.debug            = false;        % Print additional warning messages
options.saveRelTermMax   = false;        % Save the maximum of the relative terms added to the series
options.prec             = 'double';     % Precision of the calculations (default: 'double'). 
                                         % Both 'sym' (requires Symbolic Math Toolbox in MATLAB) and 'mp' (requires Advanpix toolbox: https://www.advanpix.com/) are supported, 
                                         % with arbitrary precision altered by Digits and mp.Digits, respectively

% General parameters in layer m
layer{m}.media       = 'fluid'; % Media; % fluid or solid/viscoelastic (Helmholtz equation or Navier equation)
layer{m}.R           = 1;       % Inner radius of layer
layer{m}.X           = [0,0,1]; % Evaluation points
layer{m}.rho         = 1000;    % Mass density
layer{m}.lossFactor  = 1000;    % Hysteretic loss factor (values around 0.001 for lightly damped materials, values around 0.01 for moderately damped materials and values around 0.1 for heavily damped materials)
layer{m}.calc_err_dc = false;   % Calculate the errors for the displacement conditions
layer{m}.calc_err_pc = false;   % Calculate the errors for the pressure conditions

% Parameters in layer m for options{i}.media = 'fluid'
layer{m}.c_f          	= 1500;       % Speed of sound
layer{m}.calc_p_0       = false;      % Toggle calculation of the far field pattern
layer{m}.calc_p       	= false;      % Toggle calculation of the scattered pressure
layer{m}.calc_dp      	= false(1,3); % Toggle calculation of the three components of the gradient of the pressure
layer{m}.calc_p_laplace	= false;      % Toggle calculation of the Laplace operator of the scattered pressure fields
layer{m}.calc_errHelm	= false;      % Toggle calculation of the errors for the Helmholtz equation
layer{m}.calc_p_inc	    = false;      % Toggle calculation of the incident pressure
layer{m}.calc_dp_inc	= false;      % Toggle calculation of the three components of the gradient of the incident pressure

% Parameters in layer m for options{i}.media = 'solid' or 'viscoelastic'
layer{m}.E                = 200e9;      % Youngs modulus
layer{m}.nu               = 0.3;        % Poisson ratio
layer{m}.calc_u           = false(1,3); % Toggle calculation of the three components of the displacement
layer{m}.calc_du          = false(3,3); % Toggle calculation of the three cartesian derivatives of the three components of the displacement [du_xdx du_xdy du_xdz; 
                                        %                                                                                                    du_ydx du_ydy du_ydz; 
                                        %                                                                                                    du_zdx du_zdy du_zdz]
layer{m}.calc_sigma       = false(1,6); % Toggle calculation of the six components of the stress field (cartesian coordinates) [sigma_xx sigma_yy sigma_zz sigma_yz sigma_xz sigma_xy]
layer{m}.calc_sigma_s     = false(1,6); % Toggle calculation of the six components of the stress field (spherical coordinates) [sigma_rr sigma_tt sigma_pp sigma_tp sigma_rp sigma_rt]
layer{m}.calc_err_navier  = false(1,2); % Toggle calculation of the errors for the Navier equation
```

It is recomended to build upon examples in the examples folder. In particular, a minimal example is given in examples/theDefaultExample.m

Some suggested extensions include
% Why has the error in the displacement condition worsen slightly?
% Expand to cylinder case?
% Add functionality to compute resonance frequencies
% Enable p_inc in multiple domains


