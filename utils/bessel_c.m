function f = bessel_c(nu,z,type)
% See continuedFractions.pdf p 348
if 0%1 %abs(x) > 1e3 % use assymptotic expansion
%     P = zeros(size(z));
%     k = 0;
%     while true
%         P_k = (-1)^k*C(nu,2*k)./(2*z).^(2*k);
%         if max(abs(P_k./P)) > eps
%             P = P + P_k;
%             k = k + 1;
%         else
%             break
%         end
%     end
%     Q = zeros(size(z));
%     k = 0;
%     while true
%         Q_k = (-1)^k*C(nu,2*k+1)./(2*z).^(2*k+1);
%         if max(abs(Q_k./Q)) > eps
%             Q = Q + Q_k;
%             k = k + 1;
%         else
%             break
%         end
%     end
%     if type == 1 % besselj
%         f = P.*cos(z-nu*pi/2-pi/4) - Q.*sin(z-nu*pi/2-pi/4);
%     else % bessely
%         f = P.*sin(z-nu*pi/2-pi/4) + Q.*cos(z-nu*pi/2-pi/4);
%     end
%     f = f.*sqrt(2/pi./z);
    if type == 1 % besselj
        f = real(hankel_c(nu,z));
    else % bessely
        f = imag(hankel_c(nu,z));
    end
%     keyboard
else
    if type == 1 % besselj
        f = besselj(nu,z);
    else % bessely
        f = bessely(nu,z);
    end
end
return
close all
W = @(nu,x) bessel_c(nu,x,1).*bessel_cDeriv(nu,x,2) - bessel_c(nu,x,2).*bessel_cDeriv(nu,x,1);
x = 10.^linspace(0,30,10000);
nu = 1;
loglog(x,abs(W(nu,x)-2./(pi*x))./abs(W(nu,x)))
hold on
loglog(x,abs(W(nu,x)))
loglog(x,2./(pi*x))
legend('error','W','2/pix')
% loglog(x,abs(bessel_c(nu,x,1).*bessel_cDeriv(nu,x,2)))
% loglog(x,abs(bessel_c(nu,x,2).*bessel_cDeriv(nu,x,1)))



function c = C(nu,k)

c = gamma(nu+k+1/2)/factorial(k)/gamma(nu-k+1/2);