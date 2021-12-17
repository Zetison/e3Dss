clear all
close all

startup
folderName = [homeDir '/results/e3Dss/'];
resultsFolder = [folderName '/S35d'];
if ~exist(resultsFolder, 'dir')
    mkdir(resultsFolder);
end

layer = cell(6,1);
layer{1}.R_i = 5;
layer{2}.R_i = 4.8;
layer{3}.R_i = 3;
layer{4}.R_i = 2.8;
layer{5}.R_i = 0;
layer{6}.R_i = 0;
R_i = 1.3*layer{1}.R_i;
layer = addCoating(layer,0.2);
for m = 1:numel(layer)
    R_o = R_i;
    R_i = layer{m}.R_i;
    [x,y,z] = sphere(100);
    surfObj = surf2patch(x,y,z,'triangles');
    visElements = [surfObj.faces, surfObj.faces+size(surfObj.vertices,1)];
    nodes = [R_i*surfObj.vertices; R_o*1.001*surfObj.vertices];        
    VTKdata.nodes = nodes;
    VTKdata.visElements = visElements;
    VTKoptions = struct('name',[resultsFolder '/domain' num2str(m)], 'celltype', 'VTK_WEDGE'); 
    makeVTKfile(VTKdata, VTKoptions);

end