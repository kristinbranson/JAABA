function classifyST(expdir,jab,varargin)
% classifyST(expdir,jab,p1,v1,...)
%
% expdir: (single) experiment dir
% jab: either a jab filename, or a Macguffin
%
% Optional PVs:
% - doforce (default=false): Force reclassification regardless of existing 
%   scorefiles. If false, existing scorefiles with up-to-date timestamps 
%   will not be overwritten.
% - verbose (default=0): Integer specifying verbosity level. Currently
%   expects either 0 or 1
% - jabfilename (default=''). Used only if jab is a Macguffin. Jab filename to
%   store in scorefiles.
% - useGetClassifierFromJabFile (default=false): if true, access classifier
%   using getClassifierFromJabFile (which checks for 'external'
%   classifierfiles and 'fixes' jabfiles etc), rather than reading jabfile
%   directly.
% - numTargets: see classifySTCore.

[doforce,verbose,jabfilename,useGetClassifierFromJabFile,numTargets] = myparse(varargin,...
  'doforce',false,...
  'verbose',0,...
  'jabfilename','',...
  'useGetClassifierFromJabFile',false,...
  'numTargets',1);

assert(ischar(expdir) && exist(expdir,'dir')==7,'Cannot find directory %s.',expdir);

if ischar(jab)
  Q = loadAnonymous(jab);
  if ~isempty(jabfilename)
    warning('classifyST:ignore','jabfilename ''%s'' ignored, using ''%s''.',...
      jabfilename,jab);
  end
  jabfilename = jab;
else
  assert(isa(jab,'Macguffin'),...
    'Input argument ''jab'' must either be a filename or a Macguffin obj');
  Q = jab;
  assert(ischar(jabfilename),'jabfilename must be a string.');
end

Q.modernize(true);
assert(isfield(Q.featureLexicon,'st'),'Jab must represent an ST project.');

scorefilenames = Q.file.scorefilename;
if ischar(scorefilenames)
  scorefilenames = {scorefilenames};
end
behs = Labels.verifyBehaviorNames(Q.behaviors.names); 
if useGetClassifierFromJabFile
  assert(ischar(jab),'Expected input arg ''jab'' to specify filename, since useGetClassifierFromJabFile=true.');
  classifiers = getClassifierFromJabFile(jab); 
  % note: depending on user interaction, getClassifierFromJabFile may
  % update classifier in jab. That is okay though. We have already loaded
  % the jab above (before the update), but we will not use the potentially 
  % out-of-date classifier in the already-loaded jab.
else
  classifiers = {Q.classifierStuff.params};
end
ppParams = {Q.classifierStuff.postProcessParams};
clsTS = getstructarrayfield(Q.classifierStuff,'timeStamp','numericscalar',true);
scoreNorm = getstructarrayfield(Q.classifierStuff,'scoreNorm','numericscalar',true);
assert(iscellstr(scorefilenames) && iscellstr(behs) && iscell(classifiers) && iscell(ppParams));
assert(isequal(numel(scorefilenames),numel(behs),numel(classifiers),...
  numel(ppParams),numel(clsTS),numel(scoreNorm)));
Ncls = numel(classifiers);

usePastOnly = Q.extra.usePastOnly;

for i = 1:Ncls
  scorefn = scorefilenames{i};
  
  if ~doforce 
    % skip this classifier if scorefile exists and is up-to-date
    sfn = fullfile(expdir,scorefn);
    if numel(sfn) < 4 || ~strcmp(sfn(end-3:end),'.mat'),
      sfn = [sfn '.mat']; %#ok<*AGROW>
    end
    if exist(sfn,'file'),
      tmp = load(sfn);
      if tmp.timestamp == clsTS(i),
        [~,name] = fileparts(expdir);
        if verbose>0
          fprintf('Not classifying experiment %s for scorefile %s, up-to-date predictions already exist\n',name,scorefn);
        end
        continue;
      elseif tmp.timestamp > clsTS(i)
        warningNoTrace('classifyST:scorefileTimestampAhead',...
          'Scorefile %s has timestamp more recent than classifier.',scorefn);
      end
    end
  end
  
  if verbose>0
    fprintf('Classifying experiment %s: %s.\n',expdir,scorefn);
  end
  
  allScores = classifySTCore(expdir,classifiers{i},ppParams{i},...
    scoreNorm(i),'usePastOnly',usePastOnly,'numTargets',numTargets);
  
  sf = ScoreFile;
  sf.allScores = allScores;
  sf.behaviorName = behs{i};
  sf.timestamp = clsTS(i);
  sf.jabFileNameAbs = jabfilename;
  sf.version = Q.version;

  savefilename = fullfile(expdir,scorefn);
  save(sf,savefilename);
end
