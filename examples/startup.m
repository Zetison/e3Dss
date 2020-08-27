addpath ..
addpath ../utils
addpath ../models
folderName = '../../../results/e3Dss';
if ~exist(folderName, 'dir')
    error('The folder in which results should be stored does not exist. Please make such a folder and alter the variable folderName in startup.m accordingly.')
end