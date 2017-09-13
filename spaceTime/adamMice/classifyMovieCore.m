function allScores = classifyMovieCore(expdir,classifier,postProcessParams,scoreNorm,varargin)
% allScores = classifyMovieCore(expdir,classifier,postProcessParams,scoreNorm,varargin)
%
% expdir: (single) experiment directory 
% classifier: (single) classifier
% postProcessParams: postprocessing parameters struct
% scoreNorm: scalar value for this classifier
%
% Optional PVs: 
%   - usePastOnly. Defaults to false.

assert(ischar(expdir)&&exist(expdir,'dir')==7);
assert(isstruct(classifier));
assert(isstruct(postProcessParams)&&isscalar(postProcessParams)...
  &&isfield(postProcessParams,'method'));
validateattributes(scoreNorm,{'numeric'},{'scalar'});

[usePastOnly,featureWindowRadius,chunkSize] = myparse(varargin,...
  'usePastOnly',false,...
  'featureWindowRadius',2,...
  'chunkSize',200);

persistent dogenfeatures;

if ~exist(fullfile(expdir,'features.mat'),'file')
  
  assert(false,'XXXAL multiclassifier unsupported codepath');

  % give a warning that features do not exist
  if ~isdeployed && isempty(dogenfeatures),
    res = questdlg(sprintf('features.mat file has not been generated for experiment %s. Perhaps setUpDir has not been run. Attempt to generate features now?',expdir),'Generate features?');
    if strcmpi(res,'yes'),
      dogenfeatures = true;
    else
      return;
    end
  else
    fprintf('features.mat file has not been generated for experiment %s. Attemptting to generate features now.\n',expdir);
  end
    
  if isFrontSideExpDir(expdir)
    error('Experiment %s is a front/side experiment. features.mat file must already be generated.',expdir);
  end
  
  vidfile = fullfile(expdir,Q.x.file.moviefilename);
    
  [readframe,nframes,fid,headerinfo] = get_readframe_fcn(vidfile);
  
  readFrameFcns.readframe = readframe;
  readFrameFcns.nframes = nframes;
  readFrameFcns.headerinfo = headerinfo;
  
  T = load(fullfile(expdir,Q.x.file.trxfilename));
  ips(1,:) = T.trx(1).arena.food; % hardcoded; needs to match genFeatures.m
  ips(2,:) = T.trx(1).arena.mouth;
  ips(3,:) = T.trx(1).arena.perch;
  ips = round(ips);
  face = T.trx(1).arena.face;
  
  scores = zeros(1,nframes);
  nchunks = ceil( (nframes-1)/chunkSize);
  scorecell = cell(1,nchunks);
  for ndx = 1:nchunks
    t0 = (ndx-1)*chunkSize+1;
    t1 = min(ndx*chunkSize,nframes);
    scorecell{ndx} = zeros(1,t1-t0+1);
  end
  
  for ndx = 1:nchunks
    t0 = (ndx-1)*chunkSize+1;
    t1 = min(ndx*chunkSize,nframes);
    [F,H] = genFeatures('ts',t0:t1,'readframeFcns',readFrameFcns,...
      'ips',ips,'face',face);
    F = reshape(F,[],t1-t0+1)';
    H = reshape(H,[],t1-t0+1)';
    curFeatures = [F H];
    curS = myBoostClassify(curFeatures,classifier);
    scorecell{ndx} = curS';
    
    for cndx = 1:ndx
      t0 = (cndx-1)*chunkSize+1;
      t1 = min(cndx*chunkSize,nframes);
      scores(t0:t1) = scorecell{cndx};
    end
    
    P = struct('allScores',[],....
      'classifierfilename',jabfile,...
      'timestamp',0);
    
    [t0s t1s] = get_interval_ends(scores>0);
    
    P.allScores.scores{1} = scores;
    P.allScores.tStart(1) = 1;
    P.allScores.tEnd(1) = nframes;
    P.allScores.t0s{1} = t0s;
    P.allScores.t1s{1} = t1s;
    P.allScores.postprocessed{1} = scores > 0;
    P.allScores.scoreNorm = 10; %#ok<STRNU>
  end
  
  for ndx = 1:nchunks
    t0 = (ndx-1)*chunkSize+1;
    t1 = min(ndx*chunkSize,nframes);
    scores(t0:t1) = scorecell{ndx};
  end
  
  if fid>0,
    fclose(fid);
  end

else
  
  S = load(fullfile(expdir,'features.mat'));
  scores = myBoostClassify(S.curFeatures,classifier);  
end

if usePastOnly,
  if numel(scores) <= featureWindowRadius+1,
    scores(:) = 0;
  else
    scores = [zeros(featureWindowRadius+1,1);scores(1:end-featureWindowRadius-1)];
  end  
end

posts = ApplyPostprocessing.PostProcess(scores,postProcessParams,scoreNorm);
[t0s,t1s] = get_interval_ends(posts);

allScores = struct();
for ndx = 1:3 % hardcoded, 3 "flies"
  allScores.scores{ndx} = scores;
  allScores.tStart(ndx) = 1;
  allScores.tEnd(ndx) = numel(scores);
  allScores.t0s{ndx} = t0s;
  allScores.t1s{ndx} = t1s;
  allScores.postprocessed{ndx} = posts;
end
allScores.scoreNorm = 10; 

