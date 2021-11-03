function dZ = dbessel_s(n,z,i,Z,scaled,g)
% Returns the derivative of the n'th spherical bessel function of kind i
% evaluated at every element in z
if nargin < 5
    scaled = false;
end
if scaled
    dZ = n*Z{i,1} - z.*g{i}.*Z{i,2};
else
    if isa(z,'sym') || isa(z,'mp')
        tiny = realmin('double');
    else
        tiny = eps;
    end
    if i == 1
        indices = logical(abs(z) < tiny);
        z(indices) = 1; % Avoid symbolic division by zero. These values will be replaced anyways
    end
    dZ = n./z.*Z{i,1} - g{i}.*Z{i,2};
    if i == 1
        if n == 1
            dZ(indices) = ones(1,class(z))/3;
        else
            dZ(indices) = 0;
        end
    end
end
if ~(isa(z,'sym') || isa(z,'mp'))
    if any(any(abs(dZ) > 1e290))
        error('e3Dss:infBessel','A Bessel function evaluation was too large')
    end
end