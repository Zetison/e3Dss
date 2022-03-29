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
        pt = audioread('miscellaneous/sonar.wav');
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
    case {5,6,7}                
        Fs = 44100; % sampling rate
        
        T = 16; %60/f_c
        
        dt = 1/Fs;
        t = (0:dt:(T-dt)).';

        pt = Pt_inc_(t,0,omega_c,NaN,P_inc,type,terms);
        
        N = T*Fs;
        p = T*ifft(pt,N);
        p = p(1:(numel(p)/2));
        p = p(2:end);
        if isrow(omega)
            p = p.';
        end
    case 8        
        pulseDuration = 1; % pulse duration
        
        fc = omega_c/(2*pi);
        t_s = 0.3*pulseDuration;
        t_p = pulseDuration;
        p = integral1(omega,t_s,fc,t_s) - integral1(omega,0,fc,t_s) ...
            + integral3(omega,fc,t_s,t_p) ...
            + integral2(omega,t_p,fc,t_s,t_p) - integral2(omega,t_p-t_s,fc,t_s,t_p); 
        % p = integral(@(t) 1/2*(1 + cos(pi*(1 + t/t_s))).*sin(2*pi*f_CW*t).*exp(1i*omega(1)*t),0,t_s)
        %    +integral(@(t) sin(2*pi*f_CW*t).*exp(1i*omega(1)*t),t_s,t_p-t_s)
        %    +integral(@(t) 1/2*(1 - cos(pi*(t - t_p)/t_s)).*sin(2*pi*f_CW*t).*exp(1i*omega(1)*t),t_p-t_s,t_p)
end


function I = integral1(omega,t,f,s)
p = 1;
I = 1/4*exp(1i*t*omega).*((2*(2*pi*f*cos(2*pi*f*t) - 1i*omega*sin(2*pi*f*t)))./(omega.^2 - 4*pi^2*f^2) ...
    + (s*(pi*(2*f*s - 1)*cos((pi*t*(2*f*s - 1))/s) - 1i*s*omega*sin((pi*t*(2*f*s - 1))/s)))./(pi^2*(1 - 2*f*s)^2 - s^2*omega.^2)  ...
    + (s*(pi*(2*f*s + 1)*cos((pi*t*(2*f*s + 1))/s) - 1i*s*omega*sin((pi*t*(2*f*s + 1))/s)))./((2*pi*f*s + s*(-omega) +pi).*(2*pi*f*s + s*omega+pi)));

I(omega == 2*pi*f) = -((16*f^2*s^2 - 1)*(exp(4*1i*pi*f*t) - 4*1i*pi*f*t) + 4*1i*f*s*(16*f^2*s^2 + exp(4*1i*pi*f*t) - 1)*sin(pi*t/s) - 16*f^2 *s^2*exp(4*1i*pi*f*t)*cos(pi*t/s))/(16*pi*f*(4*f*s - 1)*(4*f*s + 1));
I(omega == -2*pi*f) = (exp(-4*1i*pi*f*t)*(-1i*(16*f^2*s^2 - 1)*(4*pi*f*t*exp(4*1i*pi*f*t) - 1i) + 4*1i*f*s*(1 + (16*f^2*s^2 - 1)*exp(4*1i*pi*f*t))*sin(pi*t/s) + 16*f^2*s^2*cos(pi*t/s)))/(16*pi*f*(4*f*s - 1)*(4*f*s + 1));
I(omega == pi*(2*f+1/s))  = 1/32*((2*s*exp(2*1i*pi*t*(2*f + 1/s)) - 4*1i*pi*t*(2*f*s + 1))/(2*pi*f*s + pi) + (16*s*exp(1i*pi*t*(2*f + 1/s))*(2*f*s*cos(2*pi*f*t) - 1i*(2*f*s + 1)*sin(2*pi*f*t)))/(4*pi*f*s + pi) + (exp(1i*pi*t*(2*f + 1/s))*((1 - 2*f*s)*cos((pi*t*(2*f*s - 1))/s) + 1i*(2*f*s + 1)*sin((pi*t*(2*f*s - 1))/s)))/(pi*f));
I(omega == -pi*(2*f+1/s)) = 1/32*((2*s*exp(-(2*1i*pi*t*(2*f*s + 1))/s) + 4*1i*pi*t*(2*f*s + 1))/(2*pi*f*s + pi) + (16*s*exp(-1i*pi*t*(2*f + 1/s))*(2*f*s*cos(2*pi*f*t) + 1i*(2*f*s + 1)*sin(2*pi*f*t)))/(4*pi*f*s+ pi) + (exp(-1i*pi*t*(2*f + 1/s))*((1 - 2*f*s)*cos((pi*t*(2*f*s - 1))/s) - 1i*(2*f*s + 1)*sin((pi*t*(2*f*s - 1))/s)))/(pi*f));
function I = integral2(omega,t,f,s,p)
I = 1/4*exp(1i*t*omega).*((2*(2*pi*f*cos(2*pi*f*t) - 1i*omega*sin(2*pi*f*t)))./(omega.^2 - 4*pi^2*f^2) ...
    + (s*(pi*(2*f*s - 1)*cos((pi*(t*(2*f*s - 1) + p))/s) - 1i*s*omega*sin((pi*(t*(2*f*s - 1) + p))/s)))./(pi^2*(1 - 2*f*s)^2 - s^2*omega.^2) ...
    + (s*(pi*(2*f*s + 1)*cos((pi*(p - t*(2*f*s + 1)))/s) + 1i*s*omega*sin((pi*(p - t*(2*f*s + 1)))/s)))./((2*pi*f*s + s*(-omega) + pi).*(2*pi*f*s + s*omega+pi)));

I(omega == 2*pi*f) = (4*1i*f*s*(16*f^2*s^2 + exp(4*1i*pi*f*t) - 1)*sin((pi*(p - t))/s) + 16*f^2*s^2*exp(4*1i*pi*f*t)*cos((pi*(p - t))/s) - (16*f^2*s^2 - 1)*(exp(4*1i*pi*f*t) - 4*1i*pi*f*t))/(16*pi*f*(4*f*s - 1)*(4*f*s + 1));
I(omega == -2*pi*f) = -(exp(-4*1i*pi*f*t)*(4*1i*f*s*(1 + (16*f^2*s^2 - 1)*exp(4*1i*pi*f*t))*sin((pi*(p - t))/s) - 16*f^2*s^2*cos((pi*(p - t))/s) + (16*f^2*s^2 - 1)*(1 + 4*1i*pi*f*t*exp(4*1i*pi*f*t))))/(16*pi*f*(4*f*s - 1)*(4*f*s + 1));
I(omega == pi*(2*f+1/s))  = 1/32*(-(2*1i*((s*exp((2*1i*pi*t*(2*f*s + 1))/s) + 2*1i*pi*t*(2*f*s + 1))*sin(pi*p/s) + (2*pi*t*(2*f*s + 1) + 1i*s*exp((2*1i*pi*t*(2*f*s + 1))/s))*cos(pi*p/s)))/(2*pi*f*s + pi) + (exp((1i*pi*t*(2*f*s + 1))/s)*((1 - 2*f*s)*cos((pi*(t*(2*f*s - 1) + p))/s) + 1i*(2*f*s + 1)*sin((pi*(t*(2*f*s - 1) + p))/s)))/(pi*f) + (16*s*exp((1i*pi*t*(2*f*s + 1))/s) *(2*f*s*cos(2*pi*f*t) - 1i*(2*f*s + 1)*sin(2*pi*f*t)))/(4*pi*f*s+ pi));
I(omega == -pi*(2*f+1/s)) = 1/32*((2*s*exp(-(2*1i*pi*t*(2*f*s + 1))/s)*(cos(pi*p/s) + 1i*sin(pi*p/s)))/(2*pi*f*s + pi) + (exp(-(1i*pi*t*(2*f*s + 1))/s)*((1 - 2*f*s)*cos((pi*(t*(2*f*s - 1) + p))/s) - 1i*(2*f*s + 1)*sin((pi*(t*(2*f*s - 1) + p))/s)))/(pi*f) + (16*s*exp(-(1i*pi*t*(2*f*s + 1))/s)*(2*f*s*cos(2*pi*f*t) + 1i*(2*f*s + 1)*sin(2*pi*f*t)))/(4*pi*f*s+ pi) + 4*t*sin(pi*p/s) + 4*1i*t*cos(pi*p/s));

function I = integral3(omega,f,s,p)
I = (exp(1i*omega*(p - s)).*(2*pi*f*cos(2*pi*f*(p - s)) - 1i*omega*sin(2*pi*f*(p - s))) + exp(1i*s*omega).*(-2*pi*f*cos(2*pi*f*s) + 1i*omega*sin(2*pi*f*s)))./(omega.^2 - 4*pi^2*f^2);
I(omega == 2*pi*f) = (4*1i*pi*f*(p - 2*s) - exp(4*1i*pi*f*(p - s)) + exp(4*1i*pi*f*s))/(8*pi*f);
I(omega == -2*pi*f) = (-4*1i*pi*f*(p - 2*s) - exp(-4*1i*pi*f*(p - s)) + exp(-4*1i*pi*f*s))/(8*pi*f);
