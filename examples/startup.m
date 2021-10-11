addpath ..
addpath ../utils
addpath ../models
homeDir = expanduser('~');
folderName = [homeDir '/results/e3Dss'];
if ~exist(folderName, 'dir')
    error('The folder in which results should be stored does not exist. Please make such a folder and alter the variable folderName in startup.m accordingly.')
end

if ~isfile('../miscellaneous/U_pol.mat')
    i_max = 100;
    U_pol = cell(i_max,1);
    for i = 1:i_max
        U_pol{i} = U_p(i-1,U_pol);
    end
    save('../miscellaneous/U_pol.mat','U_pol')
%     load('../miscellaneous/U_pol.mat')
    
    U_p(1,U_pol)-[-5,0,3,0]/24
    U_p(2,U_pol)-[385,0,-462,0,81,0,0]/1152
    U_p(3,U_pol)-[-425425,0,765765,0,-369603,0,30375,0,0,0]/414720
end


function U = U_p(k,U_pol)
if k == 0
    U = 1;
else
    U = conv([-0.5,0,0.5,0,0], polyder(U_pol{k}));
    if U(1) == 0
        U = U(2:end);
    end
    U = U + 1/8*polyint(conv([-5,0,1],U_pol{k}));
end
end