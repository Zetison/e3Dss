close all
clear all %#ok

startup

if false
    layer = setHetmaniukParameters();
    BC = 'SSBC';

    % k = (6:0.05:18)';
    k = 18;
    omega = k*layer{1}.c_f;
    f = omega/(2*pi);
    d_vec = [0,0,1].';
    applyLoad = 'mechExcitation';
    % applyLoad = 'planeWave';
    options = struct('applyLoad', applyLoad, ...
                     'd_vec', d_vec, ...
                     'r_s', layer{1}.R, ...
                     'BC', BC, ...
                     'omega', omega, ...
                     'Display','iter', ...
                     'nu_a', 100, ...
                     'debug','on', ...
                     'saveRelTermMax', true, ... 
                     'P_inc', -1);
    theta = linspace(0,pi,1000).';
    layer{1}.X = layer{1}.R*[cos(theta),zeros(size(theta)),sin(theta)];
    % layer{1}.X = layer{1}.R*[0,0,1];
    layer{1}.calc_p = true;
    [layer,~,~,relTermMaxArr] = e3Dss(layer, options);

    figure(12)
    plot(theta, real(layer{1}.p),'DisplayName','Exact')
    title('Figure 12 in Hetmaniuk2012raa')
    hold on
    % real_p_Hetmaniuk = importdata('models/Hetmaniuk2012raa/Figure12.csv');
    % plot(real_p_Hetmaniuk(:,1), real_p_Hetmaniuk(:,2),'DisplayName','Hetmaniuk2012raa')
    xlabel('theta')
    xlim([theta(1), theta(end)])
    % xlim([f(1), f(end)])
    % ylim([-200, 200])
    ylabel('Real part of pressure')
    legend('show');
end

if true
    a = 6371e3;
    noLayers = 3;
    radii = getPREMprofiles([]);
    % radii = [];
    layer = cell(noLayers,1);
    R = sort(unique([radii,linspace(0,a,noLayers)]));
    noLayers = numel(R);
    [~, rho, c_l, c_s, lossFactor] = getPREMprofiles(R(1:end-1));
    E = c_s.^2.*rho.*(3*c_l.^2-4*c_s.^2)./(c_l.^2-c_s.^2);
    nu = (c_l.^2-2*c_s.^2)./(2*(c_l.^2-c_s.^2));
    layer{1}.media = 'fluid';
    layer{1}.R = a;
    layer{1}.rho = 1.225;
    layer{1}.c_f = 343;
    for i = 2:noLayers
        j = noLayers-i+1;
        if c_s(j) == 0
            layer{i}.media = 'fluid';
            layer{i}.c_f = c_l(j);
        else
            layer{i}.media = 'solid';
            layer{i}.E = E(j);
            layer{i}.nu = nu(j);
        end
        layer{i}.R = R(j);
        layer{i}.rho = rho(j);
        layer{i}.lossFactor = lossFactor(j,:);
    end

    layer{1}.X     	= a*[-1,0,0];       % Evaluation points
    layer{1}.calc_p_0 = true; % Calculate the far field pattern
    f_max = 1e2;
    f_max = 1;
    nFreq = 10;
    f = linspace(f_max/nFreq,f_max,nFreq);
    omega = 2*pi*f;
%     applyLoad = 'mechExcitation';
    applyLoad = 'planeWave';

    options = struct('applyLoad', applyLoad, ...
                     'd_vec', [1,0,-1].'/sqrt(2), ...
                     'r_s', 6368e3, ...
                     'Display','iter', ...
                     'BC','NNBC',...
                     'omega', omega);

    layer = e3Dss(layer,options); % Compute solution
    figure
    plot(f,abs(layer{1}.p_0))
    return
    % plot(r,rho)
    % for i = 1:numel(radii)
    %     hold on
    %     yLim = ylim;
    %     plot(radii(i)*[1,1],yLim,'--','color','black')
    % end
    r = linspace(0,a,10000);
    [radii, rho, c_l, c_s, lossFactor] = getPREMprofiles(r);
    figure
    plot((a-r)/1e3,rho,'DisplayName','$$\rho$$')
    hold on
    plot((a-r)/1e3,c_s,'DisplayName','$$c_s$$')
    xlim([0,1e3])
    legend('show','interpreter','latex')
    figure
    plot((a-r)/1e3,c_l,'DisplayName','$$c_l$$')
    xlim([0,1e3])
    legend('show','interpreter','latex')
end



function [radii, rho, c_l, c_s, lossFactor] = getPREMprofiles(r)
a = 6371e3;
radii = [0, 1221.5, 3480, 3630, 5600, 5701, 5771, 5971, 6151, 6291, 6346.6, 6356, 6368, 6371]*1e3;
if isempty(r)
    return
end
coeffs_rho = [13.0885,  0,          -8.8381,    0;
              12.5815,  -1.2638,    -3.6426,    -5.5281;
              7.9565,   -6.4761,    5.5283,     -3.0807;
              7.9565,   -6.4761,    5.5283,     -3.0807;
              7.9565,   -6.4761,    5.5283,     -3.0807;
              5.3197,   -1.4836,    0,          0;
              11.2494,  -8.0298,    0,          0;
              7.1089,   -3.8045,    0,          0;
              2.6910,   0.69274,    0,          0;
              2.6910,   0.69274,    0,          0;
              2.900,    0,          0,          0;
              2.600,    0,          0,          0;
              1.020,    0,          0,          0]*1e3;
coeffs_c_l = [11.2622,  0,          -6.3640,    0;
              11.0487,  -4.03262,   4.8023,     -13.5732;
              15.3891,  -5.3181,    5.5242,     -2.5514;
              24.9520,  -40.4673,   51.4832,    -26.6419;
              29.2766,   -23.6027,  5.5242,     -2.5514;
              19.0957,   -9.8672,    0,          0;
              39.7027,  -32.6166,    0,          0;
              20.3926,  -12.2569,    0,          0;
              0.8317,   7.2180,     0,          0;
%               3.5908,   4.6172,     0,          0;
              0.8317,   7.218,      0,          0;
%               3.5908,   4.6172,     0,          0;
              6.800,     0,          0,          0;
              5.800,     0,          0,          0;
              1.450,     0,          0,          0]*1e3;
coeffs_c_s = [3.6678,   0,          -4.4475,    0;
              0,        0,          0,          0;
              6.9254,   1.4672,    -2.0834,     0.9783;
              11.1671,  -13.7818,   17.4575,    -9.2777;
              22.3459,   -17.2473,   -2.0734,    0.9783;
              9.9839,   -4.9324,    0,          0;
              22.3512,  -18.5856,    0,          0;
              8.9496,   -4.4597,    0,          0;
              5.8582,    -1.4678,     0,          0;
%               -1.0839,   5.7176,     0,          0;
              5.8582,   -1.4678,      0,          0;
%               -1.0839,   5.7176,     0,          0;
              3.900,  	0,          0,          0;
              3.200,  	0,          0,          0;
              0,        0,          0,          0]*1e3;
if true % Use approximation for the effective isotropic velocities
	coeffs_c_l(9,:) = [4.1875,  3.9382, 0, 0]*1e3;
	coeffs_c_s(9,:) = [2.1519,  2.3481, 0, 0]*1e3;
end

Q_mu = [84.6, Inf, 312, 312, 312, 143, 143, 143, 80, 600, 600, 600, Inf];
Q_K =  [1327.7, 57823, 57823, 57823, 57823, 57823, 57823, 57823, 57823, 57823, 57823, 57823, 57823];
rho = zeros(size(r));
c_l = zeros(size(r));
c_s = zeros(size(r));
lossFactor = zeros(numel(r),2);
for i = 1:numel(radii)-1
    indices = and(radii(i) <= r, r < radii(i+1));
    x = r(indices)/a;
    for j = 1:4
        rho(indices) = rho(indices) + coeffs_rho(i,j)*x.^(j-1);
        c_l(indices) = c_l(indices) + coeffs_c_l(i,j)*x.^(j-1);
        c_s(indices) = c_s(indices) + coeffs_c_s(i,j)*x.^(j-1);
    end
    lossFactor(indices,1) = 1./Q_mu(i);
    lossFactor(indices,2) = 1./Q_K(i);
end

end
















































