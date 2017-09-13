function setuppaths()

subdirs_j = {'misc','filehandling','perframe'};
subdirs = {'oly','analysis'};

p = fileparts(mfilename('fullpath'));
pj = fileparts(fileparts(p));
for i = 1:numel(subdirs),
  addpath(fullfile(p,subdirs{i}));
end
for i = 1:numel(subdirs_j),
  addpath(fullfile(pj,subdirs_j{i}));
end
addpath(genpath(p));

%addpath(genpath('..'));
