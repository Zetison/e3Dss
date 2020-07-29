function [xx,yy,zz] = sphere(varargin)
%SPHERE Generate sphere.
%   [X,Y,Z] = SPHERE(N) generates three (N+1)-by-(N+1)
%   matrices so that SURF(X,Y,Z) produces a unit sphere.
%
%   [X,Y,Z] = SPHERE uses N = 20.
%
%   SPHERE(N) and just SPHERE graph the sphere as a SURFACE
%   and do not return anything.
%
%   SPHERE(AX,...) plots into AX instead of GCA.
%
%   See also ELLIPSOID, CYLINDER.

%   Clay M. Thompson 4-24-91, CBM 8-21-92.
%   Copyright 1984-2002 The MathWorks, Inc. 

% ***NOTE***: Edited for usage in Paraview

% Parse possible Axes input
narginchk(0,2);
[cax,args,nargs] = axescheck(varargin{:});

n = 20;
if nargs > 0
    n = args{1}; 
end

phi = linspace(0,2*pi-1e-10,2*n+1);
theta = linspace(0,pi,n+1).';
sintheta = sin(theta); 
sinphi = sin(phi); 

x = sintheta*cos(phi);
y = sintheta*sinphi;
z = cos(theta)*ones(1,2*n+1);

if nargout == 0
    cax = newplot(cax);
    surf(x,y,z,'parent',cax)
else
    xx = x; yy = y; zz = z;
end
