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
    U_11 = [-1556514560957130010757145625/31553580184752881664,0, ...
                   3424332034105686023665720375/10517860061584293888,0, ...
                   -365967912305800454531456575/389550372651270144,0, ...
                   201734750392525792544487385/129850124217090048,0, ...
                   -11694306169843138084657687/7213895789838336,0, ...
                   164293183874328160710877/148434069749760,0, ...
                   -23186185730591085896097833/46756731971174400,0, ...
                   527174389121818780771231/3710851743744000,0, ...
                   -329641577686894230674187/13469017440256000,0, ...
                   241770821762631191867/107752139522048,0, ...
                   -26416375998266454375/314460325543936,0, ...
                   1212400457192925/2199023255552,0,0,0,0,0,0,0,0,0,0,0];
    (U_p(11,U_pol)-U_11)./U_11
                   

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