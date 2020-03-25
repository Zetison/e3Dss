# Exact 3D scattering solutions
Exact 3D scattering solutions for spherical symmetric scatterers computes the solution to scattering problems on multilayered spherical (elastic) shells impinged by a plane wave or a wave due to a point source. The code is written in MATLAB and may for that reason be improved compuationally.

## Boundary conditions
The following boundary conditions can be used for the innermost layer:
- 'SHBC' (Sound Hard Boundary Condition) simulates a rigid boundary
- 'SSBC' (Sound Soft Boundary Condition) simulates the innermost layer having p=0
- 'NNBC' (Neumann-Neumann boundary condition) simulates full acoustic structure interaction

## Loads (Incident wave type)
The following loads can be used for the outermost layer
- 'planeWave' simulates a plane incident wave
- 'pointCharge' simulates an incident wave from a point charge (the location of the point source is given by options.d_vec*options.r_s)
- 'mechExcitation' simulates a point force mechanical excitation at options.d_vec*options.r_s with amplitude P_inc (options.r_s must be at layer{m}.R_i for a given layer m)
- 'surfExcitation' simulates a surface excitation over the region theta \in options.theta_s and phi \in [0, 2*pi] at r_s with amplitude P_inc (direction of z-axis given by options.d_vec)
- 'radialPulsation' simulates a spherical symmetric wave originating from infinity

## Parameters (with default values)
Global parameters and options
```Matlab
options.d_vec         = [0;0;1];      % Direction of the incident wave
options.omega         = 2*pi*1e3;     % Angular frequency (can be a vector)
options.Eps           = eps;          % Parameter for series truncation. The summation are terminated whenever the relative contribution of the given term is less then Eps. 
                                      % If vector input is given for either X or omega, the maximal relative contribution of the given term is compared with Eps
options.N_max         = inf;          % Upper limit for the number of terms in the series
options.prec          = 'double';     % Precision of the calculations (default: 'double'). 
                                      % Both 'sym' and 'mp' are supported, with arbitrary precision altered by Digits and mp.Digits respectively
options.applyLoad     = 'planeWave';  % Incident wave type
options.r_s           = 2;            % Radius to source location for point charge incident waves
options.P_inc         = 1;            % Amplitude of incident wave at the origin. P_inc can be given as a function handle P_inc(omega) where omega is the angular frequency
options.BC            = 'SHBC';       % Boundary condition on the inermost layer 'SSBC' (Sound soft boundary condition), 'NNBC' (Neumann-Neumann boundary condition) 

% General parameters in layer m
layer{m}.media            = 'fluid'; % Media; % fluid or solid/viscoelastic (Helmholtz equation or Navier equation)
layer{m}.R_i              = 1;       % Inner radius of layer
layer{m}.X                = [0,0,1]; % Evaluation points
layer{m}.rho              = 1000;    % Mass density
layer{m}.calc_errDispCond = false;   % Calculate the errors for the displacement conditions
layer{m}.calc_errPresCond = false;   % Calculate the errors for the pressure conditions

% Parameters in layer m for options{i}.media = 'fluid'
layer{m}.c_f          	= 1500;       % Speed of sound
layer{m}.calc_p_0       = false;      % Toggle calculation of the far field pattern
layer{m}.calc_p       	= false;      % Toggle calculation of the scattered pressure
layer{m}.calc_dp      	= false(1,3); % Toggle calculation of the three components of the gradient of the pressure
layer{m}.calc_p_laplace	= false;      % Toggle calculation of the Laplace operator of the scattered pressure fields
layer{m}.calc_errHelm	= false;      % Toggle calculation of the errors for the Helmholtz equation

% Parameters in layer m for options{i}.media = 'solid' or 'viscoelastic'
layer{m}.E            = 200e9;      % Youngs modulus
layer{m}.nu           = 0.3;        % Poisson ratio
layer{m}.calc_u       = false(1,3); % Toggle calculation of the three components of the displacement
layer{m}.calc_du      = false(3,3); % Toggle calculation of the three cartesian derivatives of the three components of the displacement [du_xdx du_ydx du_zdx; 
                                    %                                                                                                    du_xdy du_ydy du_zdy; 
                                    %                                                                                                    du_xdz du_ydz du_zdz]
layer{m}.calc_sigma   = false(1,6); % Toggle calculation of the six components of the stress field (cartesian coordinates) [sigma_xx sigma_yy sigma_zz sigma_yz sigma_xz sigma_xy]
layer{m}.calc_sigma_s = false(1,6); % Toggle calculation of the six components of the stress field (spherical coordinates) [sigma_rr sigma_tt sigma_pp sigma_tp sigma_rp sigma_rt]
layer{m}.calc_errNav  = false;      % Toggle calculation of the errors for the Navier equation
```

It is recomended to build upon examples in the examples folder. In particular, a minimal example is given in examples/theDefaultExample.m


