function FlyBowlJAABADetect(expdir,varargin)

version = '0.2';

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

%% log file

if isfield(dataloc_params,'jaabadetect_logfilestr') && ~DEBUG,
  logfile = fullfile(expdir,dataloc_params.jaabadetect_logfilestr);
  logfid = fopen(logfile,'a');
  if logfid < 1,
    warning('Could not open log file %s\n',logfile);
    logfid = 1;
  end
else
  logfid = 1;
end

timestamp = datestr(now,'yyyymmddTHHMMSS');
real_analysis_protocol = GetRealAnalysisProtocol(analysis_protocol,settingsdir);

jaabapath = fileparts(mfilename('fullpath'));
versionfile = fullfile(jaabapath,'version.txt');
jaaba_version = 'unknown';
if exist(versionfile,'file'),
  versionfid = fopen(versionfile,'r');
  if versionfid > 0,
    while true,
      s = fgetl(versionfid);
      if ~ischar(s),
        break;
      end
      s = strtrim(s);
      if ~isempty(s),
        jaaba_version = s;
        break;
      end
    end
    fclose(versionfid);
  end
end

fprintf(logfid,'\n\n***\nRunning FlyBowlJAABADetect version %s (JAABA version %s), analysis_protocol %s (linked to %s) at %s\n',version,jaaba_version,analysis_protocol,real_analysis_protocol,timestamp);

%% main loop


%jaabadir = fullfile(settingsdir,analysis_protocol,params.jaabadir);

% 
% oldpwd = pwd;
% cd(params.jaabadir);

% loop through the classifier params files to process in order
if ~iscell(params.classifierparamsfiles),
  params.classifierparamsfiles = {params.classifierparamsfiles};
end

classifierinfo = [];
for i = 1:numel(params.classifierparamsfiles),
  
  classifierparamsfile = fullfile(settingsdir,analysis_protocol,params.classifierparamsfiles{i});
  classifierinfocurr = JAABADetect(expdir,'jablistfile',classifierparamsfile,...
    'forcecompute',forcecompute,...
    'debug',DEBUG);
  %'isrelativepath',true,...
  %'fnsrelative',{'featureparamfilename'});
  classifierinfo = structappend(classifierinfo,classifierinfocurr);

end
% 
% cd(oldpwd);

%% save info

if ~DEBUG,
  
  savefile = fullfile(expdir,dataloc_params.jaabadetectinfomatfilestr);
  if exist(savefile,'file'),
    try %#ok<TRYNC>
      delete(savefile);
    end
  end
  try
    save(savefile,'classifierinfo','timestamp','analysis_protocol','real_analysis_protocol','version');
  catch ME
    warning('Could not save jaabadetect info to file %s: %s',savefile,getReport(ME));
  end
  
end

%% close log

fprintf(logfid,'Finished running FlyBowlJAABADetect at %s.\n',datestr(now,'yyyymmddTHHMMSS'));

if logfid > 1,
  fclose(logfid);
end