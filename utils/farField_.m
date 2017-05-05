function p = farField_(v,options)

options.calc_farField = true;
data = e3Dss(v,options);
p = data(1).p;