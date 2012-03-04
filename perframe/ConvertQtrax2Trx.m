function [outexpdir] = ConvertQtrax2Trx(experiment_name,rootoutputdir,varargin)

% rootoutputdir = '/groups/branson/home/bransonk/behavioranalysis/code/Jdetect/data/eric';
% experiment_subdir = 'EH111019';
% experiment_name = 'Agg_111019_GMR_38G02_AE_01_UAS_dTrpA1_2_0002_G_29_p02';

%% set parameters

% default parameters
rootdatadir = '/groups/branson/home/bransonk/behavioranalysis/code/Jdetect/data/eric';
experiment_subdir = 0;
featext = '_1_feat.mat';
roiext = '_1_roi.mat';
movieext = '.mp4';
trxfilestr = 'trx.mat';
moviefilestr = 'movie.mp4';
perframedirstr = 'perframe';
doflipud = false;
dofliplr = false;
rot = 0;
%movieinfo = struct('fps',30,'width',144,'height',144);
movieinfo = [];

[rootdatadir,experiment_subdir,featext,roiext,movieext,trxfilestr,moviefilestr,perframedirstr,...
  doflipud,dofliplr,rot,movieinfo] = ...
  myparse(varargin,...
  'rootdatadir',rootdatadir,...
  'experiment_subdir',experiment_subdir,...
  'featext',featext,...
  'roiext',roiext,...
  'movieext',movieext,...
  'trxfilestr',trxfilestr,...
  'moviefilestr',moviefilestr,...
  'perframedirstr',perframedirstr,...
  'doflipud',doflipud,...
  'dofliplr',dofliplr,...
  'rot',rot,...
  'movieinfo',movieinfo);

if isempty(experiment_subdir),
  m = regexp(experiment_name,'^Agg_(\d{6})','once','tokens');
  if ~isempty(m),
    experiment_subdir = m{1};
  end
elseif ~ischar(experiment_subdir),
  experiment_subdir = '';
end

%% names of files

featname = fullfile(rootdatadir,experiment_subdir,[experiment_name,featext]);
roiname = fullfile(rootdatadir,experiment_subdir,[experiment_name,roiext]);
moviename = fullfile(rootdatadir,experiment_subdir,[experiment_name,movieext]);

if ~exist(featname,'file'),
  error('Feature file %s does not exist',featname);
end
if ~exist(roiname,'file'),
  error('ROI file %s does not exist',roiname);
end
if ~exist(moviename,'file'),
  error('Movie file %s does not exist',moviename);
end

% make sure the root output directory exists
if ~exist(rootoutputdir,'dir'),
  mkdir(rootoutputdir);
end

% create the experiment directory
outexpdir = fullfile(rootoutputdir,experiment_name);
if ~exist(outexpdir,'dir'),
  mkdir(outexpdir);
end

% create the trx file
outmatname = fullfile(outexpdir,trxfilestr);
cadabra2ctrax(featname,roiname,moviename,outmatname,doflipud,dofliplr,rot,movieinfo);

% copy over the movie
if ~exist(moviename,'file'),
  error('Movie %s does not exist',moviename);
end
outmoviename = fullfile(outexpdir,moviefilestr);
if isunix && makesoftlink,
  cmd = sprintf('ln -s %s %s',moviename,outmoviename);
  unix(cmd);
else
  [success,msg] = copyfile(moviename,outmoviename);
  if ~success,
    error('Error copying file %s to %s: %s',moviename,outmoviename,msg);
  end
end

% save the per-frame features
perframedir = fullfile(outexpdir,perframedirstr);
if ~exist(perframedir,'dir'),
  mkdir(perframedir);
end

feat = load(featname);
obj = [feat.fly_feat.obj1,feat.fly_feat.obj2];

% allocate variables that will be saved
data = cell(1,2);
units = parseunits('unit');

% frames for which fly is tracked
firstframes = nan(1,2);
lastframes = nan(1,2);
for fly = 1:2,
  % (x,y) = (0,0) for untracked frames
  badframes = obj(fly).pos_x == 0 & obj(fly).pos_y == 0;
  lastframe = find(~badframes,1,'last');
  firstframe = find(~badframes,1);
  if isempty(lastframe), 
    firstframe = 1;
    lastframe = 0;
  end
  firstframes(fly) = firstframe;
  lastframes(fly) = lastframe;
end

% joint features
fns_ignore = {'ind1_count'
  'ind2_count'
  'obj1'
  'obj2'
  'frame'
  'time'
  };
fns = setdiff(fieldnames(feat.fly_feat),fns_ignore);
nframesall = numel(feat.fly_feat.frame);
for i = 1:numel(fns),
  fn = fns{i};
  n = numel(feat.fly_feat.(fn));
  offend = ceil((nframesall-n)/2);
  offstart = nframesall-n-offend;
  for fly = 1:2,
    data{fly} = feat.fly_feat.(fn)(firstframes(fly)+offstart:lastframes(fly)-offend);
  end
  save(fullfile(perframedir,[fn,'.mat']),'data','units');
end

% per-fly features
fns = fieldnames(obj);
for i = 1:numel(fns),
  fn = fns{i};
  for fly = 1:2,
    data{fly} = obj(fly).(fn);
  end
  save(fullfile(perframedir,[fn,'.mat']),'data','units');
end
