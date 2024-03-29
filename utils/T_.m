function T = T_(tIdx, i, n, eta, Z, g)

switch tIdx  
    case 1
        T = -n*(n+1)*Z{i,1};
    case 2
        T = -(n+1)*Z{i,1} + eta.*g{i}.*Z{i,2};   
    case 3 % -T_13
        T = -n*(n+1)*((n-1)*Z{i,1} - eta.*g{i}.*Z{i,2});   
    case 4
        T = (eta.^2-n^2+1).*Z{i,1} - eta.*g{i}.*Z{i,2};        
    case 5
        T = -n*(n+1)*((n-1)*Z{i,1} - eta.*g{i}.*Z{i,2});
    case 6 % -T_33
        T = -n*(n+1)*Z{i,1};
    case 7 % -T_43
        T = -(n^2-1-0.5*eta.^2).*Z{i,1} - eta.*g{i}.*Z{i,2}; 
    case 8
        T = n*(n+1)*((-n^2+3*n-2+eta.^2).*Z{i,1} - 4*eta.*g{i}.*Z{i,2});   
    case 9
        T = (-n^3+2*n^2+n-2+(n/2-1)*eta.^2).*Z{i,1} ...
           +(n^2+n+2-0.5*eta.^2).*eta.*g{i}.*Z{i,2};        
end