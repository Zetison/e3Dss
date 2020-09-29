function p_0 = fAnalyticPulsation_(v,C_n,y,k)
error('Depricated use functions in getAnalyticSolutions.m')
p_0 = zeros(size(v,1),1);
for i = 1:size(y,1)
    p_0 = p_0 + C_n(i)/(4*pi)*exp(-1i*k*dot3(v./repmat(norm2(v),1,3),y(i,:).'));
end