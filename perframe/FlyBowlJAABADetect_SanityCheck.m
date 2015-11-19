% make sure that running JAABADetect one classifier at a time gives the
% same results as running them in batch

function maxerr = FlyBowlJAABADetect_SanityCheck(expdir,classifierparamsfiles,varargin)

[classifierparams,isrelativepath] = ...
  myparse(varargin,'classifierparams',[],'isrelativepath',true);

if isempty(classifierparams),
  classifierparams = ReadClassifierParamsFile(classifierparamsfiles,'isrelativepath',isrelativepath);
end
nbehaviors = numel(classifierparams);

allScores = cell(1,nbehaviors);
maxerr = zeros(1,nbehaviors);
for i = 1:nbehaviors,
  
  behavior = classifierparams(i).behaviors.names{1};
  
  scorefile = fullfile(expdir,classifierparams(i).file.scorefilename);
  if ~exist(scorefile,'file'),
    warning('Score file %s does not exist',scorefile);
    continue;
  end

  allScores(i) = JAABADetect(expdir,'classifierfiles',{classifierparams(i).classifierfile},...
    'configfiles',{classifierparams(i).configfile},'debug',true,'forcecompute',true);
  
  scoresload = load(scorefile);
  nflies = numel(allScores{i});
  nflies2 = numel(scoresload.allScores.scores);
  if nflies ~= nflies2,
    warning('Number of flies does not match for %s',behavior);
    maxerr(i) = inf;
  end
  nflies = min(nflies,nflies2);
  for fly = 1:nflies,
    if numel(scoresload.allScores.scores{fly}) ~= numel(allScores{i}{fly}),
      warning('Number of frames does not match for %s fly %d\n',behavior,fly);
      maxerr(i) = inf;
    end
    realidx = ~isnan(scoresload.allScores.scores{fly});
    if ~all(~realidx == isnan(allScores{i}{fly})),
      fprintf('%s fly %d: nans do not match\n',behavior,fly);
      maxerr(i) = inf;
      continue;
    end
    maxdiff = max(abs(scoresload.allScores.scores{fly}(realidx) - allScores{i}{fly}(realidx)));
    if maxdiff > 0,
      fprintf('%s fly %d: scores do not match, maxdiff = %f\n',behavior,fly,maxdiff);
      maxerr(i) = max(maxerr(i),maxdiff);
    end
  end  
  
end

