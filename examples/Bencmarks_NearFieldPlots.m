close all
clear all %#ok

addpath ..
addpath ../utils
addpath ../models

startMatlabPool

ESBC = 0;
SSBC = 0;
for SHBC = 0 %[0, 1]
    for modelCell = {'IL'} %{'S5', 'S35', 'S135'}
        model = modelCell{1};
        switch model
            case 'S5'
                setS5Parameters
            case 'S35'
                setS35Parameters
            case 'S135'
                setS135Parameters
            case 'IL'
                setIhlenburgParameters
        end
        R_i = R_o - t;
        defineBCstring
        k_arr = [1.892808062465205, 1.8929];
        plotRealPart = 1;
        for k = 2 % 30/R_o(1)
            omega = k*c_f(1);
            % rho_f(end) = 0;
            R_a = 2*5; %R_o(1);
            d_vec = [1, 0, 0]';  
            usePointChargeWave = false;
            options = struct('d_vec', d_vec, ...
                             'omega', omega, ...
                             'R_i', R_i, ...
                             'R_o', R_o, ...
                             'P_inc', P_inc, ...
                             'E', E, ...
                             'nu', nu, ...
                             'rho_s', rho_s, ...
                             'rho_f', rho_f, ...
                             'c_f', c_f, ...
                             'SHBC', SHBC, ...
                             'ESBC', ESBC, ...
                             'SSBC', SSBC, ...
                             'plotTimeOscillation', 0, ...
                             'plotInTimeDomain', 0, ...
                             'usePointChargeWave', usePointChargeWave, ...
                             'usePlaneWave', ~usePointChargeWave, ...
                             'R_a', R_a);

%                     folderName = ['../results/paraviewResults/' model '/'];
            folderName = ['otherFunctions/e3Dss/results/paraviewResults/' model '/'];
            if ~exist(folderName, 'dir')
                mkdir(folderName);
            end
            vtfFileName = [folderName BC '_ka_' num2str(k*R_o(1))];

            extraPts = 40;

            createParaviewFiles_exact3(extraPts, vtfFileName, options)
        end
    end
end

% %         solid = getSphericalShellData(5, 4.992,'Zaxis');
% %         solid = getSphericalShellData(3, 2.98,'Zaxis');
%         solid = getSphericalShellData(1, 0.95,'Zaxis');
%         [nodes, ~, visElements] = buildVisualization3dMesh_new3(solid.knots{1}, solid.knots{2}, solid.knots{3}, 200, 200, 0, solid);
%         VTKoptions = struct('name','otherFunctions/e3Dss/results/sources/S1_solid', 'celltype', 'VTK_HEXAHEDRON'); 
%         VTKdata.nodes = nodes;
%         VTKdata.visElements = visElements;
%         VTKdata.omega = 1;
%         makeVTKfile(VTKdata, VTKoptions);