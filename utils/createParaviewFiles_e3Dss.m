function tasks = createParaviewFiles_e3Dss(extraPts, vtfFileName, layer, options)

R_a = options.R_a;
applyLoad = options.applyLoad;
M = length(layer);
plotTimeOscillation = options.plotTimeOscillation;
plotInTimeDomain = options.plotInTimeDomain;
plotDisplacementVectors = options.plotDisplacementVectors;
compDisplacementDers = options.compDisplacementDers;
if ~plotInTimeDomain
    omega = options.omega;
end
plotP_inc = ~(strcmp(applyLoad,'mechExcitation') || strcmp(applyLoad,'surfExcitation'));

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
R = Inf;
for m = 1:M
    R_o = R;
    R = layer{m}.R;
    switch layer{m}.media
        case 'fluid'
            layer{m}.calc_p = true;
            layer{m}.calc_dp = plotDisplacementVectors*[1,1,1];
            layer{m}.calc_p_inc = plotP_inc;
            layer{m}.calc_dp_inc = or(plotDisplacementVectors*[1,1,1],plotP_inc);
            if m == 1
                [visElements{m}, nodes{m}] = meshRectangleWcircHole([-s*R_a, -R_a],[s*R_a, R_a],R,h);
                nodes{m} = [nodes{m}, zeros(size(nodes{m},1),1)];
            elseif R == 0
                [visElements{m}, nodes{m}] = mesh2DDisk(R_o,h);
                nodes{m} = [nodes{m}, zeros(size(nodes{m},1),1)];
            else
                [visElements{m}, nodes{m}] = mesh2DDonut(R_o,R,h);
                nodes{m} = [nodes{m}, zeros(size(nodes{m},1),1)];
            end
        case 'solid'
            if computeForSolidDomain
                layer{m}.calc_sigma = true(1,6);
                layer{m}.calc_sigma_s = [1,0,0,0,0,0];

                layer{m}.calc_u = plotDisplacementVectors*[1,1,1]; 

                layer{m}.calc_du = compDisplacementDers*true(3,3); 
                if layer{M}.R == 0 && m == M
                    [visElements{m}, nodes{m}] = mesh2DDisk(R_o,h);
                    nodes{m} = [nodes{m}, zeros(size(nodes{m},1),1)];
                else
                    [x,y,z] = sphere(round(2*pi*R_o/h));
                    [st,sp] = size(x);
                    [~, A] = orthogonalTransform([x(:),y(:),z(:)], options.d_vec);
                    X = (A*[x(:),y(:),z(:)].').';
                    surfObj = surf2patch(reshape(X(:,1),st,sp),reshape(X(:,2),st,sp),reshape(X(:,3),st,sp),'triangles');
                    visElements{m} = [surfObj.faces, surfObj.faces+size(surfObj.vertices,1)];
                    nodes{m} = [R*surfObj.vertices; R_o*surfObj.vertices];
                end
            end
    end
end
if plotInTimeDomain
    [visTemp, temp_nodes] = meshRectangle([-s*R_a, 0],[s*R_a, 0.5*R_a],round(extraPts*s*R_a));
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
if ~strcmp(options.Display,'none')
    disp(['Time spent building mesh ' num2str(toc) ' seconds.'])
end

P_inc = options.P_inc;
if plotTimeOscillation
    N = 30;
    N_fine = 30;
    type = 1;
    startIdx = 1;
    options.P_inc = 1;
elseif plotInTimeDomain
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
    if ~strcmp(options.Display,'none')
        totnpts = npts*N*4
        npts
    end
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
switch applyLoad
    case {'planeWave','radialPulsation'}
        m_s = 1;
    case {'pointCharge','mechExcitation','surfExcitation'}
        if ~isfield(options,'r_s')
            options.r_s = 2*layer{1}.R;
        end
        r_s = options.r_s;
        m_s = 1;
        for m = 1:M
            if r_s < layer{m}.R
                m_s = m_s + 1;
            else
                break
            end
        end
end
tic
layer = e3Dss(layer, options);
if ~strcmp(options.Display,'none')
    disp(['Time spent on computing solution: ' num2str(toc) ' seconds.'])
end
[pathstr,filename] = fileparts(vtfFileName);

for m = 1:length(nodes)
    tic
    clear VTKdata
    if strcmp(layer{m}.media,'solid') && computeForSolidDomain
        toggleJacobianMatrix = layer{m}.calc_du;
        VTKoptions = struct('name',[pathstr '/solid' num2str(m) filename], 'celltype', 'VTK_WEDGE', 'plotTimeOscillation', plotTimeOscillation, ...
                    'plotSphericalRadialDisplacement',0, 'plotDisplacementVectors',plotDisplacementVectors,'plotSphericalStress_rr',1,'plotVonMisesStress',1,...
                     'plotdu_xdx',layer{m}.calc_du(1,1),...
                     'plotdu_xdy',layer{m}.calc_du(1,2),...
                     'plotdu_xdz',layer{m}.calc_du(1,3),...
                     'plotdu_ydx',layer{m}.calc_du(2,1),...
                     'plotdu_ydy',layer{m}.calc_du(2,2),...
                     'plotdu_ydz',layer{m}.calc_du(2,3),...
                     'plotdu_zdx',layer{m}.calc_du(3,1),...
                     'plotdu_zdy',layer{m}.calc_du(3,2),...
                     'plotdu_zdz',layer{m}.calc_du(3,3)); 
        if layer{M}.R == 0 && m == M
            VTKoptions.celltype = 'VTK_TRIANGLE';
        end
        if VTKoptions.plotDisplacementVectors || VTKoptions.plotSphericalRadialDisplacement
            VTKdata.displacement = transferFields(layer{m},{'u{1}','u{2}','u{3}'});
        end
        if any(toggleJacobianMatrix)
            VTKdata.jacobian = transferFields(layer{m},{'du{1,1}','du{1,2}','du{1,3}','du{2,1}','du{2,2}','du{2,3}','du{3,1}','du{3,2}','du{3,3}'});
        end
        if VTKoptions.plotVonMisesStress || VTKoptions.plotSphericalStress_rr 
            VTKdata.stress = transferFields(layer{m},{'sigma{1}','sigma{2}','sigma{3}','sigma{4}','sigma{5}','sigma{6}'});
        end
    elseif strcmp(layer{m}.media,'fluid')
        VTKoptions = struct('name',[pathstr '/fluid' num2str(m) filename], 'celltype', 'VTK_TRIANGLE', 'plotTimeOscillation', plotTimeOscillation, ...
                            'plotDisplacementVectors', plotDisplacementVectors, 'plotSphericalRadialDisplacement',0, 'plotTotField', 1, 'plotScalarField', 1, 'plotSPL', ~plotTimeOscillation);
        if ~(plotInTimeDomain || plotTimeOscillation)
            VTKoptions.plotTotFieldAbs = 1;
        end
        if VTKoptions.plotDisplacementVectors 
            rho = layer{m}.rho;
            VTKdata.displacement = transferFields(layer{m},{'dp{1}','dp{2}','dp{3}'},rho*omega.^2);
        end
        VTKdata.totField = transferFields(layer{m},{'p'});
        if m == m_s
            if plotP_inc
                VTKdata.P_inc = transferFields(layer{m},{'p_inc'});
            end
            if ~(strcmp(applyLoad,'pointExcitation') || strcmp(applyLoad,'surfExcitation') || strcmp(applyLoad,'mechExcitation'))
                if VTKoptions.plotDisplacementVectors && plotP_inc
                    VTKdata.displacement = VTKdata.displacement + transferFields(layer{m},{'dp_inc{1}','dp_inc{2}','dp_inc{3}'},rho*omega.^2);
                end
                VTKdata.scalarField = VTKdata.totField;
                VTKdata.totField = VTKdata.totField + VTKdata.P_inc;
            end
        else
            VTKoptions.plotScalarField = false;
            VTKoptions.plotP_inc = false;
        end
        if plotInTimeDomain
            VTKoptions.plotScalarField = false;
            VTKoptions.plotSPL = false;
        elseif ~plotTimeOscillation
            VTKdata.totFieldAbs = abs(VTKdata.totField);
        end
    end
    VTKdata.nodes = nodes{m};
    VTKdata.visElements = visElements{m};
    VTKdata.omega = omega;
    if plotInTimeDomain
        VTKoptions.T = T;
    end
    VTKoptions.N = N_fine;
    if ~strcmp(options.Display,'none')
        disp(['Time spent on assembling data ' num2str(toc) ' seconds.'])
    end
    if nargout > 0
        tasks(m).nodes = nodes;
        tasks(m).visElements = visElements;
        tasks(m).omega = omega;
        for fieldName = fieldnames(VTKdata).'
            tasks(m).(fieldName{1}) = VTKdata.(fieldName{1});
        end
    else
        tic
        if ~isempty(nodes{m})
            makeVTKfile(VTKdata, VTKoptions);
        end
        if ~strcmp(options.Display,'none')
            disp(['Time spent on making VTK file ' num2str(toc) ' seconds.'])
        end
    end
end


function VTKfield = transferFields(layer_m,fieldsToTransfer,scale)
if plotInTimeDomain
    eval(['temp = layer_m.' fieldsToTransfer{1} ';']);
    VTKfield = zeros(size(temp,1),1,numel(omega));
    for ii = 1:numel(fieldsToTransfer)
        eval(['temp = layer_m.' fieldsToTransfer{ii} ';']);
        if nargin > 2
            temp = temp./scale(2:end);
        end
        VTKfield(:,ii,2:end) = reshape(temp,size(temp,1),1,size(temp,2));
    end
    VTKfield = fftField(VTKfield);
else
    eval(['temp = layer_m.' fieldsToTransfer{1} ';']);
    VTKfield = zeros(size(temp,1),numel(fieldsToTransfer));
    for ii = 1:numel(fieldsToTransfer)
        eval(['VTKfield(:,ii) = layer_m.' fieldsToTransfer{ii} ';']);
    end
    if nargin > 2
        VTKfield = VTKfield./scale;
    end
    VTKfield = makeDynamic(VTKfield, VTKoptions, omega);
end
end
     
function VTKfield = fftField(VTKfield)
VTKfield = 2/T*real(fft(VTKfield,N_fine,3));
temp = VTKfield;
VTKfield(:,:,1:N_fine-startIdx+1) = temp(:,:,startIdx:end);
VTKfield(:,:,N_fine-startIdx+2:end) = temp(:,:,1:startIdx-1);
end


end
