function [success,msg,pxpermm] = ReadArenaParameters_MAGATAnalyzer(varargin)

success = false;
msg = ''; %#ok<NASGU>
SCALE = 10;

%% set parameters

[expfile,pxpermm] = myparse(varargin,...
  'expfile','',...
  'pxpermm',1);

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

try
  pxpermm = 1/expdata.experiment_1.camcalinfo.realUnitsPerPixel/SCALE;
catch ME
  msg = getReport(ME);
  return;
end

msg = sprintf('Read pxpermm = %f from expfile %s',pxpermm,expfile);
success = true;