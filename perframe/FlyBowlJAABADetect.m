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
jaabadir = fullfile(settingsdir,analysis_protocol,params.jaabadir);

% 
% oldpwd = pwd;
% cd(params.jaabadir);

% loop through the classifier params files to process in order
if ~iscell(params.classifierparamsfiles),
  params.classifierparamsfiles = {params.classifierparamsfiles};
end
for i = 1:numel(params.classifierparamsfiles),

  classifierparamsfile = fullfile(settingsdir,analysis_protocol,params.classifierparamsfiles{i});
  if ~exist(classifierparamsfile,'file'),
    error('File %s does not exist',classifierparamsfile);
  end
  fid = fopen(classifierparamsfile,'r');
  if fid < 1,
    error('Could not open file %s for reading',classifierparamsfile);
  end
  classifierfiles = {};
  configfiles = {};
  while true,
    l = fgetl(fid);
    if ~ischar(l),
      break;
    end
    if isempty(l),
      continue;
    end
    if strcmp(l(1),'%'), continue; end
    ws = regexp(l,',','split');
    classifierfiles{end+1} = fullfile(settingsdir,analysis_protocol,ws{1}); %#ok<AGROW>
    configfiles{end+1} = fullfile(settingsdir,analysis_protocol,ws{2}); %#ok<AGROW>
  end
  fclose(fid);
  
  configparams = cell(1,numel(configfiles));
  % make the featureConfig and window features paths global
  for j = 1:numel(configfiles),
    configparams{j} = ReadXMLParams(configfiles{j});
    if ~isabspath(configparams{j}.file.featureconfigfile),
      configparams{j}.file.featureconfigfile = fullfile(jaabadir,configparams{j}.file.featureconfigfile);
    end
    if ~isabspath(configparams{j}.file.featureparamfilename),
      configparams{j}.file.featureparamfilename = fullfile(jaabadir,configparams{j}.file.featureparamfilename);
    end
  end
  
  JAABADetect(expdir,'classifierfiles',classifierfiles,...
    'configparams',configparams,...
    'forcecompute',forcecompute,...
    'debug',DEBUG);
  
end
% 
% cd(oldpwd);