function [classifierinfo,allScores] = JAABADetect(expdir,varargin)
% Run classifiers trained from JAABA on experiments
% JAABADetect(expdir,'jabfiles',jabfiles)
% JAABADetect(expdir,'jablistfile',jablistfile)
%
% expdir (String or cell array of strings) -- is the experiment directory 
% or a list of experiment directories. 
% jabfiles (Cell array of strings) -- locations of the jab files for the 
% behaviors that you want to detect.
% 
% Example usage:
% JAABADetect('testExp','jabfiles',{'ChaseClassifier.jab'});

if ~isdeployed,
  SetUpJAABAPath;
end

[blockSize,jabfiles,jablistfile,forcecompute,DEBUG] = ...
  myparse(varargin,'blockSize',10000,...
  'jabfiles',{},...
  'jablistfile',0,...
  'forcecompute',false,...
  'debug',false);

if ischar(forcecompute),
  forcecompute = str2double(forcecompute) ~= 0;
end
if ischar(DEBUG),
  DEBUG = str2double(DEBUG);
end
if ~iscell(expdir)
  expdir = {expdir};
end

if ischar(jablistfile),
  jabfiles = {};
  fid = fopen(jablistfile);
  while(true)
    l = fgetl(fid);
    if ~ischar(l),
      break; 
    end
    l = strtrim(l);
    if l(1) == '#',
      continue;
    end
    jabfiles{end+1} = l; 
  end
end

if ischar(jabfiles),
  jabfiles = {jabfiles};
end

nbehaviors = numel(jabfiles);

scorefilenames = cell(1,nbehaviors);
scoresasinputs = cell(1,nbehaviors);
jabts = zeros(1,nbehaviors);
behavior = cell(1,nbehaviors);
classifierinfo = [];
for ndx = 1:nbehaviors
  Q = loadAnonymous(jabfiles{ndx});
  if isstruct(Q)
    Q = Macguffin(Q);
  end
  Q.modernize(true);
  
  assert(numel(Q.behaviors.names)==2,...
    'JAABADetect is currently unsupported for multiclassifier projects.');
  
  [~,scorefilenames{ndx},~] = fileparts(Q.file.scorefilename{1});
  scoresasinputs{ndx} = Q.scoreFeatures;
  jabts(ndx) = Q.classifierStuff.timeStamp;
  behavior{ndx} = Q.behaviors.names{1};
  classifierinfocurr = struct('jabfile',jabfiles{ndx},...
    'behavior',behavior{ndx},...
    'scorefilename',scorefilenames{ndx},...
    'timestamp',Q.classifierStuff.timeStamp,...
    'jaaba_version',Q.version);
  classifierinfo = structappend(classifierinfo,classifierinfocurr);
end

E = false(nbehaviors);
for ndx = 1:nbehaviors
  for sndx = 1:numel(scoresasinputs{ndx})
    matchndx = find(strcmp(scoresasinputs{ndx}(sndx).scorefilename,scorefilenames));
    if numel(matchndx) == 0,
      continue;
    end
    if scoresasinputs{ndx}(sndx).ts ~= jabts(matchndx)
      error('Classifier for behavior %s depends on %s generated on %s. The matching classifier:%s was trained on:%s',...
        behavior{ndx},behavior{matchndx},datestr(scoresasinputs{ndx}(sndx).ts),jabfiles{matchndx},datestr(jabts(matchndx)));
    end
    E(matchndx,ndx) = 1;
  end
end

if any(E(:))
  order = graphtopoorder(sparse(E));
else
  order = 1:nbehaviors;
end
classifierinfo = classifierinfo(order);

data = JLabelData();
data.isInteractive = false;
allScores = cell(nbehaviors,numel(expdir));
for ndx = order(:)'
  
  fprintf('Opening project %s...\n',jabfiles{ndx});
  data.openJabFileNoExps(jabfiles{ndx},false);
  for expi = 1:numel(expdir)
    if ~forcecompute
      sfn = fullfile(expdir{expi},scorefilenames{ndx});
      if ~strcmp(sfn(end-3:end),'.mat')
        sfn = [sfn '.mat']; %#ok<*AGROW>
      end
      if exist(sfn,'file'),
        Q = load(sfn);
        if Q.timestamp == jabts(ndx),
          fprintf('Skipping experiment %d for behavior %s, predictions already exist\n',expi,behavior{ndx});
          continue;
        end
      end
    end
    try 
    [success,msg] = data.AddExpDir(expdir{expi});
    if ~success,
      error(msg);
    end
    fprintf('Added experiment %d for behavior %s\n',expi,behavior{ndx});

    fprintf('Predicting on experiment %d for behavior %s\n',expi,behavior{ndx});
    if DEBUG ~= 0,
      allScores{ndx,expi} = data.PredictSaveMovie(data.nexps,0);
    else
      allScores{ndx,expi} = data.PredictSaveMovie(data.nexps);
    end
    data.RemoveExpDirs(data.nexps);
    catch ME
      fprintf('Could not classify experiment %d for %s:\n%s\n',expi,behavior{ndx},getReport(ME));
      data.closeJabFile();
      data.openJabFileNoExps(jabfiles{ndx},false);
    end
  end
  data.closeJabFile();
end
