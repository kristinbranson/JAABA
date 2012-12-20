function FlyBowlJAABADetect(expdir,varargin)

[analysis_protocol,settingsdir,datalocparamsfilestr,...
  forcecompute,DEBUG] = ...
  myparse(varargin,'analysis_protocol','current',...
  'settingsdir','/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings',...
  'datalocparamsfilestr','dataloc_params.txt',...
  'forcecompute',false,...
  'debug',false);

if isunix && settingsdir(1) ~= '/',
  warning('settingsdir path must be the global path');
end

if ischar(forcecompute),
  forcecompute = str2double(forcecompute);
end
if ischar(DEBUG),
  DEBUG = str2double(DEBUG);
end

datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);
paramsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.jaabadetectparamsfilestr);
params = ReadParams(paramsfile);
%jaabadir = fullfile(settingsdir,analysis_protocol,params.jaabadir);

% 
% oldpwd = pwd;
% cd(params.jaabadir);

% loop through the classifier params files to process in order
if ~iscell(params.classifierparamsfiles),
  params.classifierparamsfiles = {params.classifierparamsfiles};
end

for i = 1:numel(params.classifierparamsfiles),
  
  classifierparamsfile = fullfile(settingsdir,analysis_protocol,params.classifierparamsfiles{i});
  JAABADetect(expdir,'classifierparamsfile',classifierparamsfile,...
    'forcecompute',forcecompute,...
    'debug',DEBUG,...
    'isrelativepath',true,...
    'fnsrelative',{'featureconfigfile','featureparamfilename'});
  
end
% 
% cd(oldpwd);