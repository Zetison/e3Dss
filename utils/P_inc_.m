function p = P_inc_(omega,omega_c,P_inc,type,terms)
if nargin < 5
    terms = 1;
end
switch type
    case 1
        p = 4/sqrt(3)*omega_c^3./((omega.^2-omega_c^2).*(omega.^2-4*omega_c^2)).*(1-exp(-2*pi*1i*omega/omega_c));
        indices = logical( (abs(abs(omega)-abs(omega_c))/abs(omega_c) < 10*eps) ...
                          +(abs(abs(omega)-abs(2*omega_c))/abs(2*omega_c) < 10*eps));
        p(indices) = 4/(3*sqrt(3))*1i*pi./omega(indices).*exp(1i*pi*omega(indices)/omega_c);
        p = P_inc*p;
    case 2
        p = -3/2*1i*omega*omega_c^2./((omega.^2-omega_c^2).*(omega.^2-4*omega_c^2)).*(1-exp(-2*pi*1i*omega/omega_c));
        indices = logical( (abs(abs(omega)-abs(omega_c))/abs(omega_c) < 10*eps) ...
                          +(abs(abs(omega)-abs(2*omega_c))/abs(2*omega_c) < 10*eps));
        p(indices) = pi/(2*omega_c)*exp(1i*pi*omega(indices)/omega_c);
        p = P_inc*p;
    case 3
        if terms > 1
            error('not implemented')
        end
        p = omega_c*(1-exp(2*pi*1i*omega/omega_c))./(omega_c.^2-omega.^2);
        indices = logical( (abs(abs(omega)-abs(omega_c))/abs(omega_c) < 10*eps));
        p(indices) = 1i*pi/omega_c;
        p = P_inc*p;
    case 4
        pt = audioread('../miscellaneous/sonar.wav');
        T = 16;
        Fs = 44100;
        N = T*Fs;
        p = ifft(pt,N);
        p = p(1:(numel(p)/2));
        p = p(2:end);
        if isrow(omega)
            p = p.';
        end
        
% [pt fs]=wavread('signal_name.wav');
% nf=1024; %number of point in DTFT
% Y = fft(pt,nf);
% f = fs/2*linspace(0,1,nf/2+1);
% plot(f,abs(Y(1:nf/2+1)));
end