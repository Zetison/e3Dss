function S = S_(tIdx, i, n, xi, eta, Z, g)

switch tIdx
    case 1
        S = n*Z{i,1} - xi.*g{i}.*Z{i,2};    
    case 2 % Th_21, Th_31
        S = Z{i,1};           
    case 3
        S = (n^2-xi.^2-n).*Z{i,1} + 2*xi.*g{i}.*Z{i,2};   
    case 4 % T_41
        S = (n-1).*Z{i,1} - xi.*g{i}.*Z{i,2};            
    case 5 % T_11
        S = (n^2-n-0.5*eta.^2).*Z{i,1} + 2*xi.*g{i}.*Z{i,2};
    case 6 % T_31
        S = (n-0.5*eta.^2+xi.^2).*Z{i,1} - xi.*g{i}.*Z{i,2};
    case 7 % T_41
        S = (n-1).*Z{i,1} - xi.*g{i}.*Z{i,2};
    case 8
        S = (n^3-3*n^2+2*n-n/2*eta.^2+2*xi.^2).*Z{i,1} ... 
           +(-n^2-n-6+0.5*eta.^2).*xi.*g{i}.*Z{i,2};
    case 9
        S = (n^2-3*n+2-xi.^2).*Z{i,1} + 4*xi.*g{i}.*Z{i,2};
end