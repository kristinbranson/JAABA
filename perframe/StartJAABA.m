% JAABA start up script.

% Initialize all the paths.
baseDir = fileparts(pwd);
addpath(fullfile(baseDir,'misc'));
addpath(fullfile(baseDir,'filehandling'));
jlabelpath = fileparts(which('JLabel'));
addpath(fullfile(jlabelpath,'compute_perframe_features'));

try
  if matlabpool('size')<1
    matlabpool open;
  end
end
% Start JAABA.
JLabel();