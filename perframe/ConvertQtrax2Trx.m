function [outexpdir] = ConvertQtrax2Trx(experiment_name,rootoutputdir,varargin)

%rootoutputdir = '/groups/branson/home/bransonk/behavioranalysis/code/Jdetect/data/eric';
%experiment_subdir = 'EH111019';
%experiment_name = 'Agg_111019_GMR_38G02_AE_01_UAS_dTrpA1_2_0002_G_29_p02';

%% set parameters

% default parameters
rootdatadir = '/groups/branson/home/bransonk/behavioranalysis/code/Jdetect/data/eric';
experiment_subdir = 0;
featext = '_1_feat.mat';
roiext = '_1_roi.mat';
movieext = '.mp4';
trxfilestr = 'trx.mat';
moviefilestr = 'movie.mp4';
doflipud = false;
dofliplr = false;
rot = 0;
%movieinfo = struct('fps',30,'width',144,'height',144);
movieinfo = [];

[rootdatadir,experiment_subdir,featext,roiext,movieext,trxfilestr,moviefilestr,...
  doflipud,dofliplr,rot,movieinfo] = ...
  myparse(varargin,...
  'rootdatadir',rootdatadir,...
  'experiment_subdir',experiment_subdir,...
  'featext',featext,...
  'roiext',roiext,...
  'movieext',movieext,...
  'trxfilestr',trxfilestr,...
  'moviefilestr',moviefilestr,...
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
  