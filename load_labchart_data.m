%% Clear workspace
clear all force hidden
close all force hidden

%% Load the MATLAB file

%Select file
[file, path, filterindex] = ...
    uigetfile('*.mat', 'Select an MATLAB Data File', 'MultiSelect', 'off');

if (file == 0)
    return;
end

% Filename;
% Load file

filename = [path,file];
load (filename);

%% Construct paths
scriptsPath = pwd;
[homeDir, ~, ~] = fileparts(scriptsPath);
[parentDir, ~, ~] = fileparts(homeDir);
eeglabDir = fullfile(parentDir, 'eeglab14_1_2b');
dataDir = fullfile(homeDir, 'data_block1');
outputDir = fullfile(homeDir, 'output');

% Get the base of the file name
% fileBase = strrep(filename, '.mat', '');
[~, fileBase, ~] = fileparts(filename);

%% Save the data
dataFileName = fullfile(outputDir, [fileBase, '.mat']);
save(dataFileName, 'data_block1', '-v7.3');

%% Open EEGLAB
eeglab
