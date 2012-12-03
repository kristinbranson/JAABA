% set up the JAABA path

jlabelpath = fileparts(mfilename('fullpath'));
% Initialize all the paths.
baseDir = fileparts(jlabelpath);
addpath(fullfile(baseDir,'misc'));
addpath(fullfile(baseDir,'filehandling'));
addpath(fullfile(jlabelpath,'compute_perframe_features'));
addpath(fullfile(jlabelpath,'larva_compute_perframe_features'));
addpath(fullfile(baseDir,'perframe','params'));