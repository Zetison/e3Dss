addpath utils
addpath models
addpath examples
addpath examples/tests
homeDir = expanduser('~');
addpath([homeDir '/Documents/MATLAB/AdvanpixMCT-4.8.6.14636'])
folderName = [homeDir '/results/e3Dss'];
if ~exist(folderName, 'dir')
    error('The folder in which results should be stored does not exist. Please make such a folder and alter the variable folderName in startup.m accordingly.')
end
if ~isfile('miscellaneous/U_pol_double.mat') || ~isfile('miscellaneous/U_pol_mp.mat') || ~isfile('miscellaneous/U_pol_sym.mat')
    for prec = {'double', 'mp', 'sym'}
        switch prec{1}
            case 'double'
                i_max = 64;
            case 'mp'
                mp.Digits(34);
                i_max = mp('128');
            case 'sym'
                digits(32);
                i_max = vpa('128');
        end
        U_pol = cell(i_max,1);
        u_k = zeros(double(i_max),1,class(i_max));
        v_k = zeros(double(i_max),1,class(i_max));
        for i = 1:i_max
            i
            if isa(i_max,'sym')
                U_pol{i} = vpa(U_p(i-1,U_pol));
                u_k(i) = vpa(u_K(i-1));
                v_k(i) = vpa(v_K(i-1));
            else
                U_pol{i} = U_p(i-1,U_pol);
                u_k(i) = u_K(i-1);
                v_k(i) = v_K(i-1);
            end
        end
        save(['miscellaneous/U_pol_' prec{1} '.mat'],'U_pol','u_k','v_k')
    %     load(['miscellaneous/U_pol_' prec{1} '.mat'],'U_pol','u_k','v_k')
    
        % Check results against known values
        U_p1 = '[-5/24,0,3/24,0]';
        U_p2 = '[385/1152,0,-462/1152,0,81/1152,0,0]';
        U_p3 = '[-425425/414720,0,765765/414720,0,-369603/414720,0,30375/414720,0,0,0]';
        U_p11 = '[-1556514560957130010757145625/31553580184752881664,0, 3424332034105686023665720375/10517860061584293888,0, -365967912305800454531456575/389550372651270144,0, 201734750392525792544487385/129850124217090048,0, -11694306169843138084657687/7213895789838336,0,164293183874328160710877/148434069749760,0, -23186185730591085896097833/46756731971174400,0, 527174389121818780771231/3710851743744000,0, -329641577686894230674187/13469017440256000,0, 241770821762631191867/107752139522048,0, -26416375998266454375/314460325543936,0,1212400457192925/2199023255552,0,0,0,0,0,0,0,0,0,0,0]';

        switch prec{1}
            case 'mp'
                U_p(1,U_pol)-mp(U_p1)
                U_p(2,U_pol)-mp(U_p2)
                U_p(3,U_pol)-mp(U_p3)
                U_11 = mp(U_p11);
            case 'sym'
                U_p(1,U_pol)-vpa(str2sym(U_p1))
                U_p(2,U_pol)-vpa(str2sym(U_p2))
                U_p(3,U_pol)-vpa(str2sym(U_p3))
                U_11 = vpa(str2sym(U_p11));
            case 'double'
                U_p(1,U_pol)-str2num(U_p1)
                U_p(2,U_pol)-str2num(U_p2)
                U_p(3,U_pol)-str2num(U_p3)
                U_11 = str2num(U_p11);
        end
        U_11_ref = U_11;
        U_11_ref(U_11_ref == 0) = 1;
        (U_p(11,U_pol)-U_11)./U_11_ref
    end
end


function U = U_p(k,U_pol)
if k == 0
    U = ones(1,class(k));
else
    U = conv([-0.5,0,0.5,0,0], polyder(U_pol{k}));
    if U(1) == 0
        U = U(2:end);
    end
    U = U + 1/8*polyint(conv([-5,0,1],U_pol{k}));
end
end


function u = u_K(k)
u = prod((2*k+1):2:(6*k-1))./216.^k./factorial(k);
end

function v = v_K(k)
v = (6*k+1)./(1-6*k).*u_K(k);
end

function Q = polyder(P)

p = numel(P)-1;
if p == 0
    Q = zeros(1,class(P));
else
    Q = (p+1-(1:p)).*P(1:end-1);
end
end

function C = conv(A,B)

p = numel(A)-1;
q = numel(B)-1;
C = zeros(1,p+q+1,class(B));
for i = 1:p+1
    for j = 1:q+1
        k = i+j-1;
        C(k) = C(k) + A(i)*B(j);
    end
end
end
