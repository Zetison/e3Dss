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
    case 5
        % calculate an LFM chirp
        %
        % 
        %            fc : carrier (centre) frequency (hertz)
        %     bandwidth : bandwidth of the chirp  (hertz)
        % pulseDuration : duration of the pulse (s)
        %            fs : sampling frequency (hertz)
                
        Fs = 44100; % sampling rate
        B_p = 1000; % bandwidth
        pulseDuration = 1; % pulse duration
        
        T = 16; %60/f_c
        % N = 2^11*ss;
        N = T*Fs;
        B = N/T; % bandwidth
        
        dt = 1/Fs;
%         t = [0:dt:pulseDuration];
        fc = omega_c/(2*pi);
        t_s = 0.06*pulseDuration;

        p = sin(2*pi*((fc-B_p/2)*t+B_p/(2*pulseDuration).*t.^2)).*smoothing(t,0,pulseDuration,t_s);
    case 6
        % FMCW
        Fs = 44100; % sampling rate
        B_p = 200; % bandwidth
        pulseDuration = 6.85-6.3; % pulse duration
        
        T = 16; %60/f_c
        % N = 2^11*ss;
        N = T*Fs;
        B = N/T; % bandwidth
        
        dt = 1/Fs;
        fc = omega_c/(2*pi);
        t_s = 0.06*pulseDuration;
        p = zeros(size(t));
        indices = t < pulseDuration;
        p(indices) = sin(2*pi*((fc-B_p/2)*t(indices)+B_p/(2*pulseDuration).*t(indices).^2)).*smoothing(t(indices),0,pulseDuration,t_s);
        
        t_sep = 0.05;
        indices = t > pulseDuration+t_sep;
        f_CW = 3.5e3;
%         f_CW = (fc-B_p/2)+B_p/pulseDuration;
        pulseDuration2 = 1; % pulse duration
        t_s = 0.3*pulseDuration2;
        p(indices) = sin(2*pi*f_CW.*t(indices)).*smoothing(t(indices),pulseDuration+t_sep,pulseDuration+pulseDuration2+t_sep,t_s);
    case {7,8}
        Fs = 44100; % sampling rate
        B_p = 200; % bandwidth
        pulseDuration = 6.85-6.3; % pulse duration
        
        T = 16; %60/f_c
        % N = 2^11*ss;
        N = T*Fs;
        B = N/T; % bandwidth
        
        pulseDuration = 1; % pulse duration
        
        fc = omega_c/(2*pi);
%         f_CW = (fc-B_p/2)+B_p/pulseDuration;
        t_s = 0.3*pulseDuration;
        p = sin(2*pi*fc.*t).*smoothing(t,0,pulseDuration,t_s);
        
end
        
function s = smoothing(t,t_start,t_end,t_s)


s = ones(size(t));
indices = and(t_start < t, t < t_s+t_start);
s(indices) = (1+cos(pi+pi*(t(indices)-t_start)/t_s))/2;
indices = t > t_end - t_s;
s(indices) = (1-cos(pi*(t(indices)-t_end)/t_s))/2;
indices = t > t_end;
s(indices) = 0;
