function [success,msgs] = CheckExpDir(expdir)

success = true;
msgs = {};

if ~exist(expdir,'dir'),
  success = false;
  msgs{end+1} = sprintf('Directory %s does not exist',expdir);
  return;
end

% check for features file
if ~exist(fullfile(expdir,'features.mat'),'file'),
  success = false;
  msgs{end+1} = 'Features file features.mat missing';
end