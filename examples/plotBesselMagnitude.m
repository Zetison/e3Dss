close all
clear all %#ok

startup
resultsFolder = [folderName '/BesselFunctionsForLargN'];
if ~exist(resultsFolder, 'dir')
    mkdir(resultsFolder);
end
tic
npts = 41;
% npts = 200;
useHP = 0;
if useHP
    prec = 'mp';
    npts = mp(npts);
    mp.Digits(100)
else
    prec = 'double';
end
if true
    z = 10;
    n = linspace(1,1000,1000).';
    nu = n+1/2;
    J = besselj(nu,z);
    Y = bessely(nu,z);
    T = (n - bessely(nu+1,z)./bessely(nu,z)*z).^2;
    T = bessely(nu+1,z)./bessely(nu,z);
    H = J + 1i*Y;
    semilogy(nu,abs(T),'DisplayName','|T|')
    hold on
%     semilogy(nu,abs(J),'DisplayName','|J|')
%     semilogy(nu,abs(Y),'DisplayName','|Y|')
    semilogy(nu,abs(n/z),'DisplayName','n/z')
    legend show
    yLim = ylim;
    semilogy(1.6*z*[1,1],yLim)
    return
end

if false % This is too time consuming :(
    noKappa = mp('3');
    kappa = 10.^linspace(mp('1'),mp('1.01'),noKappa).';
    fractions = zeros(noKappa,3*2,prec);
    for i = 1:noKappa
        bdry = 10000;
        bdry = mp(bdry);
        N_max = 2*bdry;
        n_arr = linspace(mp('0'),N_max,round(kappa(i)));
        ReZ = bdry*linspace(mp('-1'),mp('1'),round(kappa(i)));
        ImZ = bdry*linspace(mp('-1'),mp('1'),round(kappa(i)));
        [X,Y,NN] = ndgrid(ReZ,ImZ,n_arr);
        NU = NN + 1/2;
        Z = X+1i*Y;
        zeta23 = zeta23_(Z./NU);
        temp1 = besselj(NU,Z)./exp(-NU.*zeta23);
        temp2 = bessely(NU,Z)./exp(abs(real(NU.*zeta23)));
        temp = exp(NU.*zeta23);
        indices = imag(Z) < 0;
        temp(indices) = exp(abs(real(NU(indices).*zeta23(indices))));
        temp3 = besselh(NU,1,Z)./temp;
        fractions(i,1) = min(abs(temp1(:)));
        fractions(i,2) = min(abs(temp2(:)));
        fractions(i,3) = min(abs(temp3(:)));
        fractions(i,4) = max(abs(temp1(:)));
        fractions(i,5) = max(abs(temp2(:)));
        fractions(i,6) = max(abs(temp3(:)));
        fprintf('Completed %d out of %d', i, noKappa)
    end
    semilogy(kappa,fractions)
    return
end


noXiKnots = npts;
noEtaKnots = npts;
noZetaKnots = npts;
noVisElems  = (noXiKnots-1)*(noEtaKnots-1)*(noZetaKnots-1);
visElements = zeros(noVisElems,8,prec);
eVis = 1;
for k = 1:noZetaKnots-1
    for j = 1:noEtaKnots-1
        for i = 1:noXiKnots-1
            visElements(eVis,1) = i   +   (j-1)*noXiKnots +   (k-1)*noEtaKnots*noXiKnots;
            visElements(eVis,2) = i+1 +   (j-1)*noXiKnots +   (k-1)*noEtaKnots*noXiKnots;
            visElements(eVis,3) = i+1 +       j*noXiKnots +   (k-1)*noEtaKnots*noXiKnots;
            visElements(eVis,4) = i   +       j*noXiKnots +   (k-1)*noEtaKnots*noXiKnots;
            visElements(eVis,5) = i   +   (j-1)*noXiKnots +       k*noEtaKnots*noXiKnots;
            visElements(eVis,6) = i+1 +   (j-1)*noXiKnots +       k*noEtaKnots*noXiKnots;
            visElements(eVis,7) = i+1 +       j*noXiKnots +       k*noEtaKnots*noXiKnots;
            visElements(eVis,8) = i   +       j*noXiKnots +       k*noEtaKnots*noXiKnots;

            eVis = eVis + 1;
        end
    end
end 
bdry = 10000;
% bdry = 1000;
if useHP
    bdry = mp(bdry);
    % N_max = 20*bdry;
    N_max = 2*bdry;
    n_arr = round(N_max*linspaceHP(mp('0'),mp('1'),noZetaKnots).^2);
    % n_arr = linspace(0,noZetaKnots,noZetaKnots);
    ReZ = bdry*linspace(mp('-1'),mp('1'),noXiKnots);
    ImZ = bdry*linspace(mp('-1'),mp('1'),noEtaKnots);
else
    % N_max = 20*bdry;
    N_max = 2*bdry;
%     n_arr = round(N_max*linspaceHP(0,1,noZetaKnots).^2);
    n_arr = round(linspace(0,N_max,noZetaKnots));
    % n_arr = linspace(0,noZetaKnots,noZetaKnots);
    ReZ = bdry*linspace(-1,1,noXiKnots);
    ImZ = bdry*linspace(-1,1,noEtaKnots);
end
[X,Y,NN] = ndgrid(ReZ,ImZ,n_arr);
NU = NN + 1/2;
Z = X+1i*Y;
zeta23 = zeta23_(Z./NU);
ZETA = zeta_(Z./NU,zeta23);
if 0
    [Xi,Yi,Zi] = ndgrid((1:noXiKnots)/noXiKnots,(1:noEtaKnots)/noEtaKnots,(1:noZetaKnots)/noZetaKnots);
else
    Xi = X;
    Yi = Y;
    Zi = NN;
end
VTKdata.nodes = [Xi(:), Yi(:), Zi(:)];

VTKdata.visElements = visElements;

testField = zeros(size(VTKdata.nodes,1),9,prec);
temp = besselj(NU,Z);
if 0
    temp2 = log10(abs(temp));
    temp2(temp2 == Inf) = 3.082547155599167e+02;
    temp2(temp2 == -Inf) = -3.082547155599167e+02;
    imagesc(ImZ,n_arr,reshape(temp2(2,:,:),noEtaKnots,noZetaKnots).')
    set(gca,'YDir','normal')
    return
end
testField(:,1) = temp(:);

temp = bessely(NU,Z);
testField(:,2) = temp(:);

temp = besselh(NU,1,Z);
% temp = besselj(NU,Z) + 1i*bessely(NU,Z);
testField(:,3) = temp(:);

temp = exp(-NU.*zeta23); % nu*zeta23_(x/nu)
testField(:,4) = temp(:);

temp = exp(abs(real(NU.*zeta23))); % -abs(real(nu*zeta23_(x/nu)))
testField(:,5) = temp(:);

% temp = exp(NU.*zeta23); % -nu*zeta23_(x/nu)
% temp = exp(NU.*zeta23); % -nu*zeta23_(x/nu)
% temp = exp(-NU.*zeta23) + 1i*exp(abs(real(NU.*zeta23))); % -nu*zeta23_(x/nu)
% temp = exp(-real(NU.*zeta23)).*(exp(-1i*imag(NU.*zeta23)) + 1i*exp(abs(real(NU.*zeta23)) + real(NU.*zeta23))); % -nu*zeta23_(x/nu)
temp = exp(NU.*zeta23);
indices = imag(Z) < 0;
temp(indices) = exp(abs(real(NU(indices).*zeta23(indices))));
testField(:,6) = temp(:);

temp = exp( (NN+0.5).*zeta23_(Z./(NN+0.5)) - (NN+1.5).*zeta23_(Z./(NN+1.5)))./(Z./(2*NN));
testField(:,7) = temp(:);
temp = exp( -abs(real((NN+0.5).*zeta23_(Z./(NN+0.5)))) + abs(real((NN+1.5).*zeta23_(Z./(NN+1.5)))))./((2*NN)./abs(Z));
testField(:,8) = temp(:);
temp = exp( -(NN+0.5).*zeta23_(Z./(NN+0.5)) + (NN+1.5).*zeta23_(Z./(NN+1.5)))./((2*NN)./Z);
testField(:,9) = temp(:);

temp = testField(:,1)./testField(:,4);
testField(:,10) = temp(:);
temp = testField(:,2)./testField(:,5);
testField(:,11) = temp(:);
temp = testField(:,3)./testField(:,6);
testField(:,12) = temp(:);

if useHP
    for i = 1:3
        j = i+9;
        fprintf('%20g < |Z_n^{(%d)}/s_n^{(%d)}| < %g\n', min(abs(testField(:,j))), i, i, max(abs(testField(:,j))))
    end
    s1 = sprintf('%.3e',min(abs(testField(:,10:12)),[],'all'));
    s2 = sprintf('%.3e',max(abs(testField(:,10:12)),[],'all'));
    fprintf('\\num{%s} \\lessapprox |Z_n^{(i)}(z)/s_n^{(i)}(z)| \\lessapprox \\num{%s}\n', s1, s2)
    s1 = sprintf('%.3e',max(abs(testField(:,1:3)),[],'all'));
    fprintf('\\max_{i,n,z}|Z_n^{(i)}(z)| \\approx \\num{%s}\n', s1)
    toc
    return
end

% testField(isinf(testField(:,1)),1) = realmax;
% testField(isinf(testField(:,2)),2) = realmax;
% testField(:,1:2) = log10(abs(testField(:,1:2)));

testField(:,1:6) = log10(abs(testField(:,1:6)));
testField(:,10:12) = log10(abs(testField(:,10:12)));
testField(testField == Inf) = 3.082547155599167e+02;
testField(testField == -Inf) = -3.082547155599167e+02;
VTKdata.testField = testField;
VTKoptions = struct('name',[resultsFolder '/besselj_3D'], 'celltype', 'VTK_HEXAHEDRON','plotTestField',1); 
                
                
makeVTKfile(VTKdata, VTKoptions);
% 
% 
% noVisElems  = (noXiKnots-1)*(noEtaKnots-1);
% visElements = zeros(noVisElems,4);
% eVis = 1;
% for j = 1:noEtaKnots-1
%     for i = 1:noXiKnots-1
%         visElements(eVis,1) = i   +   (j-1)*noXiKnots;
%         visElements(eVis,2) = i+1 +   (j-1)*noXiKnots;
%         visElements(eVis,3) = i+1 +       j*noXiKnots;
%         visElements(eVis,4) = i   +       j*noXiKnots;
%         eVis = eVis + 1;
%     end
% end
% VTKdata.visElements = visElements;
% VTKdata.nodes = [Xi(:), Yi(:), zeta_(X(:)+1i*Y(:))];
% VTKoptions = struct('name',[resultsFolder '/besselj_3D'], 'celltype', 'VTK_HEXAHEDRON','plotTestField',1); 
% 
