function [success,msg,fps] = ReadFPS_MAGATAnalyzer(varargin)

success = false;
msg = ''; %#ok<NASGU>

%% set parameters

[expfile,fps] = myparse(varargin,...
  'expfile','',...
  'fps',30);

%% check for files

if isempty(expfile),
  msg = 'Experiment file is empty';
  return;
end

if ~exist(expfile,'file'),
  msg = sprintf('Experiment file %s does not exist',expfile);
  return;
end

%% load in experiment data

try
  expdata = load(expfile);
catch ME,
  msg = getReport(ME);
  return;
end

%% read

timestamps = expdata.experiment_1.elapsedTime';
fps = 1/nanmedian(diff(timestamps));
success = true;
msg = sprintf('Read fps = %f from expfile %s',fps,expfile);