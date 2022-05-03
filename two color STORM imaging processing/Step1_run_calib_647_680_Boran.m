%Folder of original dax file, folder of bin file, bin file name, subfolder of beads tform in original dax file. 
% this code can split channel, fit molecules, and then assign the molecules
% to two channels
clear all;
close all;
addpath('C:\Users\Boran\Desktop\matlab-storm\Functions\Analysis\');
addpath('C:\Users\Boran\Desktop\matlab-storm\Functions\Calibration\');
addpath('C:\Users\Boran\Desktop\matlab-storm\Functions\DataTypes\');
addpath('C:\Users\Boran\Desktop\matlab-storm\Functions\IO\');
addpath('C:\Users\Boran\Desktop\matlab-storm\Functions\Misc\');
%c:\MinGW64\mingw\bin\;c:\MinGW64\mingw\lib\;c:\Users\Hazen\lib\;
global daoSTORMexe;
daoSTORMexe='path=C:\Python27\;C:\Users\Hazen\storm-analysis\storm_analysis\windows_dll\; && set PYTHONPATH=%PYTHONPATH%;C:\Users\Hazen\storm-analysis\;  && python.exe C:\Users\Hazen\storm-analysis\storm_analysis\daostorm_3d\mufit_analysis.py';
dataPath = 'P:\20190214 clathrin ms 647 add rab 680\DIV34\axon\';
analysisPath = dataPath;
pixelSize = 167; 
threshold=10000;

dataFiles = dir(fullfile(strcat([dataPath],'*con*.dax')));

for ii=1:length(dataFiles)
    function_split_two_channel_storm3(dataPath, dataFiles(ii).name)
end

dataFiles = dir(fullfile(strcat([dataPath  'split\'],'*con*.dax')));

for ii=1:length(dataFiles)
    analysisLabel = '';
    desiredFrame = 1;
    
        parametersName = ['daoSTORMParameters_' analysisLabel '.xml'];
        
        WriteDaoSTORMParameters([analysisPath parametersName], 'max_frame',desiredFrame,...
        'iterations', 1, ... % Do not fit overlapping molecules
        'threshold', threshold, ...
        'sigma', 1.6, 'frame_step',500,... % The best guess for the PSF size in terms of pixels
        'pixel_size', pixelSize, 'drift_correction', true);

%% Run analysis of all fov

        daoSTORM([dataPath 'split\' dataFiles(ii).name], ... % Path to dax to analyze
        [analysisPath parametersName], ... % Path to daoSTORM configuration file
        'overwrite', true, ... % Overwrite any existing analysis
        'numParallel', 1, ... % Only process one file at a time
        'savePath', [analysisPath 'split\'], ... % Location of the mlist files
        'mListType', [analysisLabel 'mlist'], ... % A label to mark the analysis for each file
        'outputInMatlab', true); % Display the output of daoSTORM in matlab
  
%         mListFiles = dir(fullfile([dataPath dataFiles(ii).name(1:end-4) '_mlist.bin']));
       
  
end

for ii=5%:length(dataFiles)
    run_calib_647_680([dataPath, 'split\'], dataFiles(ii).name(1:end-6));
end

