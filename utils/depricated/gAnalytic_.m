function dp = gAnalytic_(v,options)
error('Depricated, use e3Dss directly instead')
data = e3Dss(v,options);
dp = [data(1).dpdx, data(1).dpdy, data(1).dpdz];

