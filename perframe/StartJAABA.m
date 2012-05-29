% JAABA start up script.

% Initialize all the paths.
baseDir = fileparts(pwd);
addpath(fullfile(baseDir,'misc'));
addpath(fullfile(baseDir,'filehandling'));
jlabelpath = fileparts(which('JLabel'));
addpath(fullfile(jlabelpath,'compute_perframe_features'));

try
  if matlabpool('size')<1
    c=parcluster;
    matlabpool('open',c.NumWorkers-1);  % BJA: must save one for cache_thread
  end
end
% Start JAABA.
JLabel();
