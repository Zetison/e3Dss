function h = dhankel_s(n,z,type)
%Returns the derivative of the n'th spherical bessel function of kind "type" 
%evaluated at every element in z

if type == 1
    h = dbessel_s(n,z,1,NaN) + 1i*dbessel_s(n,z,2,NaN);
else
    h = dbessel_s(n,z,1,NaN) - 1i*dbessel_s(n,z,2,NaN);
end