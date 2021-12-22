clear all
close all

startup
folderName = [homeDir '/results/e3Dss/'];
resultsFolder = [folderName '/S35d'];
if ~exist(resultsFolder, 'dir')
    mkdir(resultsFolder);
end

layer = cell(6,1);
layer{1}.R = 5;
layer{2}.R = 4.8;
layer{3}.R = 3;
layer{4}.R = 2.8;
layer{5}.R = 0;
layer{6}.R = 0;
R = 1.3*layer{1}.R;
layer = addCoating(layer,0.2);
for m = 1:numel(layer)
    R_o = R;
    R = layer{m}.R;
    [x,y,z] = sphere(100);
    surfObj = surf2patch(x,y,z,'triangles');
    visElements = [surfObj.faces, surfObj.faces+size(surfObj.vertices,1)];
    nodes = [R*surfObj.vertices; R_o*1.001*surfObj.vertices];        
    VTKdata.nodes = nodes;
    VTKdata.visElements = visElements;
    VTKoptions = struct('name',[resultsFolder '/domain' num2str(m)], 'celltype', 'VTK_WEDGE'); 
    makeVTKfile(VTKdata, VTKoptions);

end