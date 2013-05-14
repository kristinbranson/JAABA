function JAABADetect(expdir,varargin)
% Run classifiers trained from JAABA on experiments
% scores = JAABADetectNew(expdir,'jabfiles',jabfiles)
% scores = JAABADetectNew(expdir,'jablistfile',jablistfile)
% expdir (String) -- is the experiment directory. 
% jabfiles (Cell array of Strings) -- locations of the jab files
% for the behaviors that you want to detect.
% Example usage:
% JAABADetectNew('testExp','jabfiles',{'ChaseClassifier.jab'});

[blockSize,jabfiles,jablistfile,forcecompute,DEBUG] = ...
  myparse(varargin,'blockSize',10000,...
  'jabfiles',{},...
  'jablistfile',0,...
  'forcecompute',false,...
  'debug',false);

if ~isdeployed,
  SetUpJAABAPath;
end

if ischar(forcecompute),
  forcecompute = str2double(forcecompute) ~= 0;
end
if ischar(DEBUG),
  DEBUG = str2double(DEBUG);
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

nbehaviors = numel(jabfiles);

scorefilenames = cell(1,nbehaviors);
scoresasinputs = cell(1,nbehaviors);
jabts = zeros(1,nbehaviors);
behavior = cell(1,nbehaviors);
for ndx = 1:nbehaviors
  Q = load(jabfiles{ndx},'-mat');
  [~,scorefilenames{ndx},~] = fileparts(Q.x.file.scorefilename);
  scoresasinputs{ndx} = Q.x.scoreFeatures;
  jabts(ndx) = Q.x.classifierStuff.timeStamp;
  behavior{ndx} = Q.x.behaviors.names{1};
end

E = false(nbehaviors);
for ndx = 1:nbehaviors
  for sndx = 1:numel(scoresasinputs{ndx})
    matchndx = find(strcmp(scoresasinputs{ndx}(sndx).scorefilename,scorefilenames));
    if numel(matchndx) == 0,
      continue;
    end
    if scoresasinputs{ndx}(sndx).ts ~= jabts(matchndx)
      error('Classifier for behavior %s depends on %s generated on %s. The matching classifier:%s was trained on:',...
        behavior{ndx},behavior{matchndx},datestr(scoresasinputs{ndx}(sndx).ts),jabfiles{matchndx},datestr(jabts(matchndx)));
    end
    E(matchndx,ndx) = 1;
  end
end

order = graphtopoorder(sparse(E));

data = JLabelData();
for ndx = order(:)'
  
  data.openJabFileNoExps(jabfiles{ndx},false);
  [success,msg] = data.AddExpDir(expdir,false);
  if ~success,
    error(msg);
  end
  if ~forcecompute
    sfn = data.GetFile('scores',1);
    if exist(sfn,'file'),
      Q = load(sfn);
      if Q.timestamp == jabts(ndx),
        continue;
      end
    end
  end
  data.PredictSaveMovie(1);
end