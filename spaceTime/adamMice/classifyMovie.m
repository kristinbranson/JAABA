function classifyMovie(expdir,jab,varargin)
% classifyMovie(expdir,jab,p1,v1,...)
%
% expdir: (single) experiment dir
% jab: either a jab filename, or a Macguffin
% Optional PVs:
% - doforce (default=false): Force reclassification regardless of existing 
% scorefiles. If false, existing scorefiles with up-to-date timestamps will 
% not be overwritten.
% - verbose (default=0): Integer specifying verbosity level. Currently
% expects either 0 or 1
% - useGetClassifierFromJabFile (default=false): if true, access classifier
% using getClassifierFromJabFile (which checks for 'external'
% classifierfiles and 'fixes' jabfiles etc), rather than reading jabfile
% directly.

[doforce,verbose,useGetClassifierFromJabFile] = myparse(varargin,...
  'doforce',false,...
  'verbose',0,...
  'useGetClassifierFromJabFile',false);

assert(ischar(expdir) && exist(expdir,'dir')==7,'Cannot find directory %s.',expdir);

if ischar(jab),
  Q = loadAnonymous(jab);
else
  assert(isa(jab,'Macguffin'),...
    'Input argument ''jab'' must either be a filename or a Macguffin obj');
  Q = jab;
end

Q.modernize(true);

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
  classifiers = Q.classifierStuff.params;
end
if ~iscell(classifiers)
  classifiers = {classifiers};
end
clsTS = Q.classifierStuff.timeStamp;
scoreNorm = Q.classifierStuff.scoreNorm;
assert(iscellstr(scorefilenames) && iscellstr(behs) && iscell(classifiers));
assert(isequal(numel(scorefilenames),numel(behs),numel(classifiers),numel(clsTS),numel(scoreNorm)));
Ncls = numel(classifiers);

ppParams = Q.classifierStuff.postProcessParams;
usePastOnly = Q.extra.usePastOnly;
assert(isscalar(ppParams));

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
        warningNoTrace('classifyMovie:scorefileTimestampAhead',...
          'Scorefile %s has timestamp more recent than classifier.',scorefn);
      end
    end
  end
  
  if verbose>0
    fprintf('Classifying experiment %s: %s.\n',expdir,scorefn);
  end
  
  allScores = classifyMovieCore(expdir,classifiers{i},...
    ppParams,scoreNorm(i),'usePastOnly',usePastOnly);
  
  savefilename = fullfile(expdir,scorefn);
  writeScoreFile(savefilename,allScores,behs{i},jab,clsTS(i));  
end


