% set up the JAABA path

jlabelpath = fileparts(mfilename('fullpath'));
% Initialize all the paths.
addpath(jlabelpath);  % in case we ever want to cd out of this dir
baseDir = fileparts(jlabelpath);
addpath(fullfile(baseDir,'misc'));
addpath(fullfile(baseDir,'filehandling'));
addpath(fullfile(jlabelpath,'larva_compute_perframe_features'));
addpath(fullfile(jlabelpath,'compute_perframe_features'));
addpath(fullfile(baseDir,'perframe','params'));
addpath(fullfile(baseDir,'tests'));
