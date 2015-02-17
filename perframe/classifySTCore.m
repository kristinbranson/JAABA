function allScores = classifySTCore(expdir,classifier,postProcessParams,scoreNorm,varargin)
% allScores = classifySTCore(expdir,classifier,postProcessParams,scoreNorm,varargin)
%
% expdir: (single) experiment directory 
% classifier: (single) classifier
% postProcessParams: postprocessing parameters struct
% scoreNorm: scalar value for this classifier
%
% Optional PVs: 
%   - usePastOnly. Defaults to false.
%   - numTargets. Positive integer, number of targets to save in
%   scorefiles. Defaults to 1.

assert(ischar(expdir)&&exist(expdir,'dir')==7);
assert(isstruct(classifier));
assert(isstruct(postProcessParams)&&isscalar(postProcessParams)...
  &&isfield(postProcessParams,'method'));
validateattributes(scoreNorm,{'numeric'},{'scalar'});

[usePastOnly,featureWindowRadius,chunkSize,numTargets] = myparse(varargin,...
  'usePastOnly',false,...
  'featureWindowRadius',2,...
  'chunkSize',200,...
  'numTargets',1);

assert(exist(fullfile(expdir,'features.mat'),'file')==2,...
  'Experiment %s: st features file does not exist.');

S = load(fullfile(expdir,'features.mat'));

scores = myBoostClassify(S.curFeatures,classifier);  

if usePastOnly,
  if numel(scores) <= featureWindowRadius+1,
    scores(:) = 0;
  else
    scores = [zeros(featureWindowRadius+1,1);scores(1:end-featureWindowRadius-1)];
  end  
end

scores = scores(:)';

posts = PostProcessor.PostProcess(scores,postProcessParams,scoreNorm);
[t0s,t1s] = get_interval_ends(posts);

allScores = ScoreFile.allScrs(numTargets);
allScores.scores(:) = {scores};
allScores.tStart(:) = 1;
allScores.tEnd(:) = numel(scores);
allScores.postprocessed(:) = {posts};
allScores.postprocessparams = postProcessParams;
allScores.t0s(:) = {t0s};
allScores.t1s(:) = {t1s};
allScores.scoreNorm = 10;
