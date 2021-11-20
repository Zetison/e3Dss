close all
clear all %#ok

startup
resultsFolder = [folderName '/BesselFunctionsForLargN'];
if ~exist(resultsFolder, 'dir')
    mkdir(resultsFolder);
end
npts = 401;
npts = 100;
noXiKnots = npts;
noEtaKnots = npts;
noZetaKnots = npts;
noVisElems  = (noXiKnots-1)*(noEtaKnots-1)*(noZetaKnots-1);
visElements = zeros(noVisElems,8);
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
N_max = 20*bdry;
% N_max = bdry;
n_arr = round(N_max*linspace(0,1,noZetaKnots).^2);
ReZ = bdry*linspace(-1,1,noXiKnots);
ImZ = bdry*linspace(-1,1,noEtaKnots);
[X,Y,NN] = ndgrid(ReZ,ImZ,n_arr);
NU = NN + 1/2;
Z = X+1i*Y;
zeta23 = zeta23_(Z./NU);
ZETA = zeta_(Z./NU,zeta23);
if false
    [Xi,Yi,Zi] = ndgrid((1:noXiKnots)/noXiKnots,(1:noEtaKnots)/noEtaKnots,(1:noZetaKnots)/noZetaKnots);
else
    Xi = X;
    Yi = Y;
    Zi = NN;
end
VTKdata.nodes = [Xi(:), Yi(:), Zi(:)];

VTKdata.visElements = visElements;

testField = zeros(size(VTKdata.nodes,1),9);
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
testField(:,3) = temp(:);
temp = exp(-2/3*NU.*ZETA.^(3/2));
testField(:,4) = temp(:);
temp = exp(abs(2/3*NU.*real(ZETA.^(3/2))));
testField(:,5) = temp(:);
temp = exp(2/3*NU.*ZETA.^(3/2));
testField(:,6) = temp(:);

temp = exp( (NN+0.5).*zeta23_(Z./(NN+0.5)) - (NN+1.5).*zeta23_(Z./(NN+1.5)))./(Z./(2*NN));
testField(:,7) = temp(:);
temp = exp( -abs(real((NN+0.5).*zeta23_(Z./(NN+0.5)))) + abs(real((NN+1.5).*zeta23_(Z./(NN+1.5)))))./((2*NN)./abs(Z));
testField(:,8) = temp(:);
temp = exp( -(NN+0.5).*zeta23_(Z./(NN+0.5)) + (NN+1.5).*zeta23_(Z./(NN+1.5)))./((2*NN)./Z);
testField(:,9) = temp(:);

% testField(isinf(testField(:,1)),1) = realmax;
% testField(isinf(testField(:,2)),2) = realmax;
% testField(:,1:2) = log10(abs(testField(:,1:2)));

testField(:,1:6) = log10(abs(testField(:,1:6)));
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
