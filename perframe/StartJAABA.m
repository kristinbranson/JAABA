% JAABA start up script.
jlabelpath = fileparts(mfilename('fullpath'));
% Initialize all the paths.
baseDir = fileparts(jlabelpath);
addpath(fullfile(baseDir,'misc'));
addpath(fullfile(baseDir,'filehandling'));
addpath(fullfile(jlabelpath,'compute_perframe_features'));
addpath(fullfile(jlabelpath,'larva_compute_perframe_features'));
addpath(fullfile(baseDir,'perframe','params'));

try
  c=parcluster;
  if (c.NumWorkers>2) && (matlabpool('size')<1)
    matlabpool('open',c.NumWorkers-1);  % BJA: must save one for frame cache thread
  end
end
% Start JAABA.
JLabel();
