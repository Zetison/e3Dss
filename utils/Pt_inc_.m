function p = Pt_inc_(t,z,omega_c,k_c,P_inc,type,terms)
if nargin < 7
    terms = 1;
end
switch type
    case 1
        if length(z) > 1 && length(t) == 1
            p = zeros(size(z));
            indices = logical((omega_c*t-k_c*z > 0).*(omega_c*t-k_c*z < 2*pi));

            p(indices) = 4/(3*sqrt(3))*(sin(omega_c*t-k_c*z(indices))-1/2*sin(2*(omega_c*t-k_c*z(indices))));
        elseif length(z) == 1 && length(t) > 1
            p = zeros(size(t));
            indices = logical((omega_c*t-k_c*z > 0).*(omega_c*t-k_c*z < 2*pi));

            p(indices) = 4/(3*sqrt(3))*(sin(omega_c*t(indices)-k_c*z)-1/2*sin(2*(omega_c*t(indices)-k_c*z)));
        end
        p = P_inc*p;
    case 2
        if length(z) > 1 && length(t) == 1
            p = zeros(size(z));
            indices = logical((omega_c*t-k_c*z > 0).*(omega_c*t-k_c*z < 2*pi));

            p(indices) = 1/2*(-cos(omega_c*t-k_c*z(indices))+cos(2*(omega_c*t-k_c*z(indices))));
        elseif length(z) == 1 && length(t) > 1
            p = zeros(size(t));
            indices = logical((omega_c*t-k_c*z > 0).*(omega_c*t-k_c*z < 2*pi));

            p(indices) = 1/2*(-cos(omega_c*t(indices)-k_c*z)+cos(2*(omega_c*t(indices)-k_c*z)));
        end
        p = P_inc*p;
    case 3
        b = -ones(terms-1,1);
        A = zeros(terms-1);
        for i = 1:terms-1
            A(i,:) = (2:terms).^(2*i-1);
        end
        a = A\b;
        a = [1; a];
        if length(z) > 1 && length(t) == 1
            p = zeros(size(z));
            indices = logical((omega_c*t-k_c*z > 0).*(omega_c*t-k_c*z < 2*pi));
            for i = 1:terms
                p(indices) = sum(sin(omega_c*t-k_c*z(indices))-1/2*sin(2*(omega_c*t-k_c*z(indices))));
            end
        elseif length(z) == 1 && length(t) > 1
            p = zeros(size(t));
            indices = logical((omega_c*t-k_c*z > 0).*(omega_c*t-k_c*z < 2*pi));

            for i = 1:terms
                p(indices) = p(indices) + a(i)*sin(i*(omega_c*t(indices)-k_c*z));
            end
        end
        p = P_inc*p;
    case 4
        [ys,Fs] = audioread('miscellaneous/sonar_raw.wav'); % From https://freesound.org/s/28693/
        T = 16;
        N = T*Fs;
        p = zeros(N,1);
        idx = 186894;
        ysCut = ys(idx:end,1);
        p(1:numel(ysCut)) = ysCut/max(abs(ysCut));
        audiowrite('miscellaneous/sonar.wav',p,Fs);
        p = audioread('miscellaneous/sonar.wav');
        if isrow(t)
            p = p.';
        end
end
        


