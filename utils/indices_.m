function indices = indices_(i,nu,z,nu_a)
% indices = and(nu > abs(z), and(abs(nu.*real(zeta23_(z./nu))) > log(sqrt(realmax(class(z)))), nu > nu_a));
% indices = and(nu > abs(z), and(abs(nu.*real(zeta23_(z./nu))) > log(realmax(class(z))*eps), nu > nu_a));
if nu_a == -1
    indices = false(size(z));
else
    if i == 2
        indices = and(and(abs(nu.*real(zeta23_(z./nu))) > log(sqrt(realmax(class(z)))), nu > nu_a), abs(z) < nu);
    else
        indices = and(and(abs(nu.*zeta23_(z./nu)) > log(sqrt(realmax(class(z)))), nu > nu_a), abs(z) < nu);
    end
end
% if numel(nu) ~= numel(z)
%     indices = repmat(nu > nu_a,size(z));
% else
%     indices = nu > nu_a;
% end

