% JAABA start up script.

% Initialize all the paths.
baseDir = fileparts(pwd);
addpath(fullfile(baseDir,'misc'));
addpath(fullfile(baseDir,'filehandling'));
jlabelpath = fileparts(which('JLabel'));
addpath(fullfile(jlabelpath,'compute_perframe_features'));

try
  c=parcluster;
  if (c.NumWorkers>2) && (matlabpool('size')<1)
    matlabpool('open',c.NumWorkers-1);  % BJA: must save one for frame cache thread
  end
end
% Start JAABA.
JLabel();
