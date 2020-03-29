function createParaviewFiles_e3Dss(extraPts, vtfFileName, layer, options)

ESBC = options.ESBC;
SSBC = options.SSBC;
R_a = options.R_a;
d_vec = options.d_vec;

M = length(layer);
plotTimeOscillation = options.plotTimeOscillation;
plotInTimeDomain = options.plotInTimeDomain;
plotDisplacementVectors = 0;
if ~plotInTimeDomain
    omega = options.omega;
end

if isfield(options, 'computeForSolidDomain')
    computeForSolidDomain = options.computeForSolidDomain;
else
    computeForSolidDomain = 0;
end
s = 1.3;
h = 2*s*R_a/(round(extraPts*s*R_a)-1);
tic
nodes = cell(M,1);
visElements = cell(M,1);
R_i = Inf;
for m = 1:M
    R_o = R_i;
    R_i = layer{m}.R_i;
    switch layer{m}.media
        case 'fluid'
            layer{m}.calc_p = true;
            layer{m}.calc_dp = plotDisplacementVectors*[1,1,1];
            if m == 1
                [visElements{m}, nodes{m}] = meshRectangleWcircHole([-s*R_a, -R_a],[s*R_a, R_a],R_i,h);
                nodes{m} = [nodes{m}, zeros(size(nodes{m},1),1)];
            elseif R_i == 0
                [visElements{m}, nodes{m}] = mesh2DDisk(R_o,h);
                nodes{m} = [nodes{m}, zeros(size(nodes{m},1),1)];
            else
                [visElements{m}, nodes{m}] = mesh2DDonut(R_o,R_i,h);
                nodes{m} = [nodes{m}, zeros(size(nodes{m},1),1)];
            end
        case 'solid'
            if computeForSolidDomain
                layer{m}.calc_sigma = true(1,6);
                layer{m}.calc_sigma_s = [1,0,0,0,0,0];

                layer{m}.calc_u = plotDisplacementVectors*[1,1,1]; 

                layer{m}.calc_du = false(3,3); 
                if ESBC && m == M
                    [visElements{m}, nodes{m}] = mesh2DDisk(R_o,h);
                    nodes{m} = [nodes{m}, zeros(size(nodes{m},1),1)];
                else
                    [x,y,z] = sphere(round(2*pi*R_o/h));
                    surfObj = surf2patch(x,y,z,'triangles');
                    visElements{m} = [surfObj.faces, surfObj.faces+size(surfObj.vertices,1)];
                    nodes{m} = [R_i*surfObj.vertices; R_o*surfObj.vertices];
                end
            end
    end
end
if plotInTimeDomain
    [visTemp, temp_nodes] = meshRectangle([-1.3*R_a, 0],[1.3*R_a, 0.5*R_a],round(extraPts*1.3*R_a));
    nodes{1} = [nodes{1}; temp_nodes(:,1), R_a*ones(size(temp_nodes,1),1), temp_nodes(:,2)];
    visElements{1} = [visElements{1}; visTemp+max(max(visElements{1}))];
end
for m = 1:M
    if ~isempty(nodes{m})
        layer{m}.X = nodes{m};
    end
end
npts = 0;
for i = 1:length(nodes)
    npts = npts + size(nodes{i},1);
end

disp(['Time spent building mesh ' num2str(toc) ' seconds.'])

P_inc = options.P_inc;
if plotTimeOscillation
    N = 30;
    N_fine = 30;
    type = 1;
    startIdx = 1;
    options.P_inc = 1;
elseif plotInTimeDomain
    npts
    f_c = options.f_c;
    N = options.N; % 200
    N_fine = 2*N;
    if N > 2000
        switch options.applyLoad
            case 'planeWave'
                startIdx = 1900; % 900
            case {'pointCharge','mechExcitation','surfExcitation','radialPulsation'}
                startIdx = 2000;
        end
    else
        startIdx = 1;
    end
    T = options.T; %N/M;
    B = N/T; % bandwidth
    
    f_R = B/2;
    df = 1/T;
    f = linspace(0,f_R-df,N/2);
    omega = 2*pi*f;
    options.omega = omega(2:end);
    type = 1;
    totnpts = npts*N*4
    omega_c = 2*pi*f_c;
    options.P_inc = @(omega) P_inc_(omega,omega_c,P_inc,type);
%     keyboard
else
    N = 1;
    N_fine = 1;
    type = 1;
    startIdx = 1;
    options.P_inc = 1;
end
switch options.applyLoad
    case {'planeWave','radialPulsation'}
        m_s = 1;
    case {'pointCharge','mechExcitation','surfExcitation'}
        if ~isfield(options,'r_s')
            options.r_s = 2*layer{1}.R_i;
        end
        r_s = options.r_s;
        m_s = 1;
        for m = 1:M
            if r_s < layer{m}.R_i
                m_s = m_s + 1;
            else
                break
            end
        end
end
tic
layer = e3Dss(layer, options);
disp(['Time spent on computing solution: ' num2str(toc) ' seconds.'])
[pathstr,filename] = fileparts(vtfFileName);

for m = 1:length(nodes)
    if ~(SSBC && m == length(nodes))
        tic
        clear VTKdata
        if strcmp(layer{m}.media,'solid') && computeForSolidDomain
            toggleJacobianMatrix = layer{m}.calc_du;
            if ESBC && m == M
                VTKoptions = struct('name',[pathstr '/solid' num2str(m) filename], 'celltype', 'VTK_TRIANGLE', 'plotTimeOscillation', plotTimeOscillation, ...
                            'plotSphericalRadialDisplacement',0,'plotDisplacementVectors',plotDisplacementVectors,'plotSphericalStress_rr',1,'plotVonMisesStress',1); 
            else
                VTKoptions = struct('name',[pathstr '/solid' num2str(m) filename], 'celltype', 'VTK_WEDGE', 'plotTimeOscillation', plotTimeOscillation, ...
                            'plotSphericalRadialDisplacement',0, 'plotDisplacementVectors',plotDisplacementVectors,'plotSphericalStress_rr',1,'plotVonMisesStress',0); 
            end
                    
            if VTKoptions.plotDisplacementVectors || VTKoptions.plotSphericalRadialDisplacement
                displacement = zeros(size(nodes{m},1),3,length(omega));
                for i = 2:length(omega)
                    displacement(:,:,i) = [layer{m}.u_x(:,i-1) layer{m}.u_y(:,i-1) layer{m}.u_z(:,i-1)];
                end
                if ~plotInTimeDomain
                    VTKdata.displacement = [layer{m}.u_x layer{m}.u_y layer{m}.u_z];
                    
                    temp = VTKdata.displacement;
                    VTKdata.displacement = zeros([size(VTKdata.displacement), N]);
                    for i = 1:N
                        t = (i-1)/N*2*pi/omega;
                        VTKdata.displacement(:,:,i) = real(temp*exp(-1i*omega*t));
                    end
                end
            end
            if any(toggleJacobianMatrix)
                VTKdata.jacobian = zeros(size(nodes{m},1),6,length(omega));
                
                if ~plotInTimeDomain
                    VTKdata.jacobian = [layer{m}.du_xdx layer{m}.du_xdy layer{m}.du_xdz layer{m}.du_ydx layer{m}.du_ydy layer{m}.du_ydz layer{m}.du_zdx layer{m}.du_zdy layer{m}.du_zdz];
                    
                    temp = VTKdata.jacobian;
                    VTKdata.jacobian = zeros([size(VTKdata.jacobian), N]);
                    for i = 1:N
                        t = (i-1)/N*2*pi/omega;
                        VTKdata.jacobian(:,:,i) = real(temp*exp(-1i*omega*t));
                    end
                end
            end
            if VTKoptions.plotVonMisesStress || VTKoptions.plotSphericalStress_rr 
                VTKdata.stress = zeros(size(nodes{m},1),6,length(omega));
                
                if plotInTimeDomain
                    for i = 2:length(omega)
                        VTKdata.stress(:,:,i) = [layer{m}.sigma_xx(:,i-1) layer{m}.sigma_yy(:,i-1) layer{m}.sigma_zz(:,i-1) layer{m}.sigma_yz(:,i-1) layer{m}.sigma_xz(:,i-1) layer{m}.sigma_xy(:,i-1)];
                    end
                    VTKdata.stress = 2/T*real(fft(VTKdata.stress,N_fine,3));
                    temp = VTKdata.stress;
                    VTKdata.stress(:,:,1:N_fine-startIdx+1) = temp(:,:,startIdx:end);
                    VTKdata.stress(:,:,N_fine-startIdx+2:end) = temp(:,:,1:startIdx-1);
                else
                    VTKdata.stress = [layer{m}.sigma_xx layer{m}.sigma_yy layer{m}.sigma_zz layer{m}.sigma_yz layer{m}.sigma_xz layer{m}.sigma_xy];
                    
                    temp = VTKdata.stress;
                    VTKdata.stress = zeros([size(VTKdata.stress), N]);
                    for i = 1:N
                        t = (i-1)/N*2*pi/omega;
                        VTKdata.stress(:,:,i) = real(temp*exp(-1i*omega*t));
                    end
                end
            end
        elseif strcmp(layer{m}.media,'fluid')
            VTKoptions = struct('name',[pathstr '/fluid' num2str(m) filename], 'celltype', 'VTK_TRIANGLE', ...
                             'plotTimeOscillation', plotTimeOscillation, 'plotDisplacementVectors', plotDisplacementVectors, ...
                             'plotSphericalRadialDisplacement',0, 'plotTotField', 1, 'plotScalarField', 1, 'plotSPL', 1);
            if ~(plotInTimeDomain || plotTimeOscillation)
                VTKoptions.plotTotFieldAbs = 1;
            end
            if VTKoptions.plotDisplacementVectors 
                displacement = zeros(size(nodes{m},1),3,length(omega));
            end
            VTKdata.totField = zeros(size(nodes{m},1),1,length(omega));
            rho = layer{m}.rho;
            if m == m_s     
                c_f = layer{m}.c_f;
                if plotInTimeDomain
                    VTKdata.P_inc = zeros(size(nodes{m},1),1,length(omega));
                    for i = 2:length(omega)     
                        k = omega(i)/c_f;   
                        switch options.applyLoad
                            case 'pointCharge'
                                r_s = options.r_s;
                                x_s = r_s*d_vec.';
                                r = @(v) norm2(repmat(x_s,size(v,1),1)-v);
                                Phi_k = @(v) exp(1i*k*r(v))./(4*pi*r(v));     
                                p_inc = @(v) 4*pi*P_inc_(omega(i),omega_c,P_inc,type)*Phi_k(v);
                                gp_inc = @(v) 4*pi*P_inc_(omega(i),omega_c,P_inc,type)*elementProd(Phi_k(v).*(1i*k - 1./r(v))./r(v), repmat(-r_s*d_vec.',size(v,1),1)-v);
                            case 'planeWave'
                                k_vec = options.d_vec*k;
                                p_inc = @(v) P_inc_(omega(i),omega_c,P_inc,type)*exp(1i*dot3(v, k_vec));
                                gp_inc = @(v) P_inc_(omega(i),omega_c,P_inc,type)*exp(1i*dot3(v, k_vec))*1i*k_vec';
                            case 'radialPulsation'
                                R_i = layer{m}.R_i;
                                p_inc = @(v) P_inc_(omega(i),omega_c,P_inc,type)*exp(-1i*k*(norm2(v)-R_i))./norm2(v);
                                gp_inc = @(v) P_inc_(omega(i),omega_c,P_inc,type)*elementProd(exp(-1i*k*(norm2(v)-R_i))./norm2(v).*(1i*k - 1./norm2(v))./norm2(v), v);
                        end
                        if strcmp(options.applyLoad,'pointExcitation') || strcmp(options.applyLoad,'surfExcitation') || strcmp(options.applyLoad,'mechExcitation')
                            VTKdata.totField(:,:,i) = layer{m}.p(:,i-1);
                            if VTKoptions.plotDisplacementVectors 
                                gScalarField = [layer{m}.dpdx(:,i-1) layer{m}.dpdy(:,i-1) layer{m}.dpdz(:,i-1)];
                                displacement(:,:,i) = gScalarField/(rho*omega(i)^2);
                            end
                        else
                            VTKdata.P_inc(:,:,i) = p_inc(nodes{m});
                            VTKdata.totField(:,:,i) = layer{m}.p(:,i-1) + VTKdata.P_inc(:,:,i);
                            if VTKoptions.plotDisplacementVectors 
                                gScalarField = [layer{m}.dpdx(:,i-1) layer{m}.dpdy(:,i-1) layer{m}.dpdz(:,i-1)];
                                displacement(:,:,i) = (gScalarField+gp_inc(nodes{m}))/(rho*omega(i)^2);
                            end
                        end
                    end
                else
                    k = omega/c_f; 
                    switch options.applyLoad
                        case 'pointCharge'
                            r_s = options.r_s;
                            x_s = r_s*d_vec.';
                            r = @(v) norm2(repmat(x_s,size(y,1),1)-v);
                            Phi_k = @(v) exp(1i*k*r(v))./(4*pi*r(v));     
                            p_inc = @(v) 4*pi*P_inc*Phi_k(v);
                            gp_inc = @(v) 4*pi*P_inc*elementProd(Phi_k(v).*(1i*k - 1./r(v))./r(v), repmat(-r_s*d_vec.',size(v,1),1)-v);
                        case 'planeWave'
                            k_vec = options.d_vec*k;
                            p_inc = @(v) P_inc*exp(1i*dot3(v, k_vec));
                            gp_inc = @(v) P_inc*exp(1i*dot3(v, k_vec))*1i*k_vec';
                        case 'radialPulsation'
                            R_i = layer{m}.R_i;
                            p_inc = @(v) P_inc*exp(-1i*k*(norm2(v)-R_i))./norm2(v);
                            gp_inc = @(v) P_inc*elementProd(exp(-1i*k*(norm2(v)-R_i))./norm2(v).*(1i*k - 1./norm2(v))./norm2(v), v);
                    end
                    VTKdata.P_inc = p_inc(nodes{m});
                    VTKdata.scalarField = layer{m}.p(:,1);
                    if strcmp(options.applyLoad,'pointExcitation') || strcmp(options.applyLoad,'surfExcitation') || strcmp(options.applyLoad,'mechExcitation')
                        VTKdata.totField = layer{m}.p(:,1);
                        VTKdata.displacement = [layer{m}.dpdx(:,1) layer{m}.dpdy(:,1) layer{m}.dpdz(:,1)]/(rho*omega^2);
                    else
                        VTKdata.totField = layer{m}.p(:,1) + VTKdata.P_inc;
                        VTKdata.displacement = ([layer{m}.dpdx(:,1) layer{m}.dpdy(:,1) layer{m}.dpdz(:,1)]+gp_inc(nodes{m}))/(rho*omega^2);
                    end
                end
                VTKdata = rmfield(VTKdata,'P_inc');
            else
                VTKoptions.plotScalarField = false;
                VTKoptions.plotP_inc = false;
                if plotInTimeDomain
                    for i = 2:length(omega)
                        if VTKoptions.plotDisplacementVectors 
                            gScalarField = [layer{m}.dpdx(:,i-1) layer{m}.dpdy(:,i-1) layer{m}.dpdz(:,i-1)];
                            displacement(:,:,i) = gScalarField/(rho*omega(i)^2);
                        end
                        VTKdata.totField(:,:,i) = layer{m}.p(:,i-1);
                    end
                else
                    VTKdata.totField = layer{m}.p;
                    VTKdata.displacement = [layer{m}.dpdx(:,1) layer{m}.dpdy(:,1) layer{m}.dpdz(:,1)]/(rho*omega^2);
                end
            end
            if plotInTimeDomain
                VTKoptions.plotScalarField = false;
                VTKoptions.plotSPL = false;
                VTKdata.totField = 2/T*real(fft(VTKdata.totField,N_fine,3));
                temp = VTKdata.totField;
                VTKdata.totField(:,:,1:N_fine-startIdx+1) = temp(:,:,startIdx:end);
                VTKdata.totField(:,:,N_fine-startIdx+2:end) = temp(:,:,1:startIdx-1);
            elseif plotTimeOscillation
                temp = VTKdata.totField;
                VTKdata.totField = zeros([size(VTKdata.totField), N]);
                for i = 1:N
                    t = (i-1)/N*2*pi/omega;
                    VTKdata.totField(:,:,i) = real(temp*exp(-1i*omega*t));
                end
            else
                VTKdata.totFieldAbs = abs(VTKdata.totField);
                VTKdata.totField = real(VTKdata.totField);
            end
        end
        if plotInTimeDomain && VTKoptions.plotDisplacementVectors 
            VTKdata.displacement = 2/T*real(fft(displacement,N_fine,3));
            temp = VTKdata.displacement;
            VTKdata.displacement(:,:,1:N_fine-startIdx+1) = temp(:,:,startIdx:end);
            VTKdata.displacement(:,:,N_fine-startIdx+2:end) = temp(:,:,1:startIdx-1);
        end
        VTKdata.nodes = nodes{m};
        VTKdata.visElements = visElements{m};
        VTKdata.omega = omega;
        if plotInTimeDomain
            VTKoptions.T = T;
        end
        VTKoptions.N = N_fine;
        disp(['Time spent on assembling data ' num2str(toc) ' seconds.'])
        tic
        if ~isempty(nodes{m})
            makeVTKfile(VTKdata, VTKoptions);
        end
        disp(['Time spent on making VTK file ' num2str(toc) ' seconds.'])
    end
end
       