function [targets,expdirs_chosen,framestarts,nframesperseg] = ...
  ChooseExampleIntervals(expdirs,classifierparamsfile,varargin)

if ~iscell(expdirs)
  expdirs = {expdirs};
end
nexps = numel(expdirs);

score_norm_prctile = 80;

% hysteresis for counting bouts
min_score_close = -.1;
max_dist_close = 5;
max_dist_open = 5;
max_score_open = .25;
nframesperseg = 500;
behaviorsuse = [];

% for weighing number of bouts of different behaviors
weightlogoff = .1;

% for weighting the same fly, same sex
downweight_samefly = 5;
downweight_samesex = 2;

% scheme for weighting
weightingscheme = 'sumlogsumscore';

score_thresh = .25;
[nintervals,score_norm_prctile,...
  min_score_close,max_dist_close,...
  max_dist_open,max_score_open,...
  nframesperseg,...
  downweight_samefly,...
  downweight_samesex,...
  weightlogoff,...
  behaviorsuse,...
  weightingscheme] = ...
  myparse(varargin,'nintervals',1,...
  'score_norm_prctile',score_norm_prctile,...
  'min_score_close',min_score_close,...
  'max_dist_close',max_dist_close,...
  'max_dist_open',max_dist_open,...
  'max_score_open',max_score_open,...
  'nframesperseg',nframesperseg,...
  'downweight_samefly',downweight_samefly,...
  'downweight_samesex',downweight_samesex,...
  'weightlogoff',weightlogoff,...
  'behaviorsuse',behaviorsuse,...
  'weightingscheme',weightingscheme);

classifierparams = ReadClassifierParamsFile(classifierparamsfile);
nbehaviors = numel(classifierparams);
behaviors = cell(1,nbehaviors);
scorefns = cell(1,nbehaviors);
for i = 1:nbehaviors,
  if iscell(classifierparams(i).behaviors.names),
    behaviors{i} = sprintf('%s_',classifierparams(i).behaviors.names{:});
    behaviors{i} = behaviors{i}(1:end-1);
  else
    behaviors{i} = classifierparams(i).behaviors.names;
  end
  behaviors{i} = regexprep(behaviors{i},'(_|^)([a-z])','${upper($2)}');
  scorefn = classifierparams(i).file.scorefilename;
  scorefns{i} = regexprep(scorefn,'\.mat$','');
end

%% create trx structure

%global STORED_TRX;
%if isempty(STORED_TRX),

trx = Trx('trxfilestr',classifierparams(1).file.trxfilename,...
  'moviefilestr',classifierparams(1).file.moviefilename,...
  'perframedir',classifierparams(1).file.perframedir);
if isfield(classifierparams,'perframe') && isfield(classifierparams(1).perframe,'params'),
  trx.SetPerFrameParams(classifierparams(1).perframe.params);
end  
if isfield(classifierparams,'perframe') && isfield(classifierparams(1).perframe,'landmark_params'),
  trx.SetLandmarkParams(classifierparams(1).perframe.landmark_params);
end  

for i = 1:nexps,
  trx.AddExpDir(expdirs{i},'openmovie',false);
end

% STORED_TRX = trx;
% else
%   trx = STORED_TRX;
% end

%% process scores

scores = cell(1,trx.nflies);
for fly = 1:trx.nflies,
  scores{fly} = nan(nbehaviors,trx(fly).nframes);
  for behi = 1:nbehaviors,
    scores{fly}(behi,:) = trx(fly).(scorefns{behi});
  end
end

% normalize
score_norms = nan(nbehaviors,1);
for behi = 1:nbehaviors,
  tmp = cell2mat(cellfun(@(x) x(behi,:),scores,'UniformOutput',false));
  score_norms(behi) = prctile(abs(tmp),score_norm_prctile);
end
for fly = 1:trx.nflies,
  for behi = 1:nbehaviors,
    scores{fly}(behi,:) = max(-1,min(1,scores{fly}(behi,:)/score_norms(behi,:)));
  end
end

labels = cell(1,trx.nflies);
for fly = 1:trx.nflies,
  
  labels{fly} = false(nbehaviors,trx(fly).nframes);
  for behi = 1:nbehaviors,
    label0 = scores{fly}(behi,:) > score_thresh;

    % close small holes where score >= min_score_close
    label_lenient = scores{fly}(behi,:) >= min_score_close;
    [t0s,t1s] = get_interval_ends(label0); t1s = t1s-1;
    gapdist = t0s(2:end)-t1s(1:end-1)-1;
    issmallgap = gapdist <= max_dist_close;
    label1 = label0;
    for i = find(issmallgap),
      doclose = all(label_lenient(t1s(i)+1:t0s(i+1)-1));
      if doclose,
        label1(t1s(i)+1:t0s(i+1)-1) = true;
      end
    end
    
    % remove short bouts where max score <= max_score_open
    [t0s,t1s] = get_interval_ends(label1); t1s = t1s-1;
    label_strict = scores{fly}(behi,:) > max_score_open;
    gapdist = t1s-t0s+1;
    issmallgap = gapdist <= max_dist_open;
    label2 = label1;
    for i = find(issmallgap),
      doopen = ~any(label_strict(t0s(i):t1s(i)));
      if doopen,
        label2(t0s(i):t1s(i)) = false;
      end
    end
    labels{fly}(behi,:) = label2;

  end
end

%%

% count the number of bout starts and ends in each interval
fil = ones(1,nframesperseg);
nbouts = cell(1,trx.nflies);
w0 = cell(1,trx.nflies);
for fly = 1:trx.nflies,
  
  isstart = double(~labels{fly}(:,1:end-1) & labels{fly}(:,2:end));
  isend = double(labels{fly}(:,1:end-1) & ~labels{fly}(:,2:end));
  nbouts{fly} = .5*imfilter(isstart,fil,0,'same') + ...
    .5*imfilter(isend,fil,0,'same');
  
  switch weightingscheme,
    case 'sumlogsumscore',
      logsumscores = nan(size(labels{fly}));
      minv = log(weightlogoff);
      logsumscores(:) = minv;
      for behi = 1:nbehaviors,
        [t0s,t1s] = get_interval_ends(labels{fly}(behi,:)); t1s=t1s-1;
        % csscores(i) is the sum of scores from 1:i-1
        csscores = [0,cumsum(scores{fly}(behi,:))];
        sumscores = csscores(t1s+1)-csscores(t0s);
        lsscurr = log(max(0,sumscores)+weightlogoff)./(t1s-t0s+1);
        for i = 1:numel(t0s),
          logsumscores(behi,t0s(i):t1s(i)) = lsscurr(i);
        end
      end
      logsumscores = logsumscores-minv;
      sumlogsumscore = imfilter(logsumscores,fil,0,'same');
      w0{fly} = sumlogsumscore;
    case 'lognbouts',  
      w0{fly} = log(nbouts{fly}+weightlogoff);
    case 'isbout',
      w0{fly} = double(nbouts{fly} >= 1);
    case 'none',
      w0{fly} = ones(size(nbouts{fly}));
    case 'random',
      w0{fly} = rand(size(nbouts{fly}));
  end
end

% make sure weights are non-negative
if ~strcmpi(weightingscheme,'none'),
  tmp = min(min(cell2mat(w0),[],1),[],2);
  w0 = cellfun(@(x) x-tmp,w0,'UniformOutput',false);
end
% 
% % normalize per-behavior
% Z = zeros(nbehaviors,1);
% for behi = 1:nbehaviors,
%   Z(behi) = sum(cell2mat(cellfun(@(x) x(behi,:),w0,'UniformOutput',false)));
% end
% Z = Z / max(Z);
% w0 = cellfun(@(x) bsxfun(@rdivide,x,Z),w0,'UniformOutput',false);

% make sure sex classification is uniform for each bout
ismale = cellfun(@(x) strcmpi(x,'M'),trx.sex,'UniformOutput',false);
isallowed = cell(1,trx.nflies);
for fly = 1:trx.nflies,
  nmale = imfilter(double(ismale{fly}),fil,nan,'same');
  isallowed{fly} = nmale == 0 | nmale == nframesperseg;
end

w0mat = cell2mat(w0);
w0mat(:,~cell2mat(isallowed)) = nan;
nboutsmat = cell2mat(nbouts);
ismalemat = cell2mat(ismale);
flyidx = zeros(1,numel(ismalemat));
flyoffs = cumsum([0,trx.nfliespermovie(1:end-1)]);
flystarts = cumsum([1,trx(1:trx.nflies-1).nframes]);
flyidx(flystarts) = 1;
flyidx = cumsum(flyidx);

framestarts = nan(1,nintervals);
targets = nan(1,nintervals);
expis = nan(1,nintervals);

% only care about some behaviors
if ~isempty(behaviorsuse)
  behaviorsuseidx = ismember(behaviors,behaviorsuse);
else
  behaviorsuseidx = true(1,nbehaviors);
end
for segi = 1:nintervals,
  
  % choose the highest weight interval
  [~,i] = max(sum(w0mat(behaviorsuseidx,:),1));
  ismalecurr = ismalemat(i);
  fly = flyidx(i);
  i0 = i-flystarts(fly)+1-floor(nframesperseg/2);
  t0 = i0-trx(fly).off;
  framestarts(segi) = t0;
  expi = trx.fly2exp(fly);
  flyi = fly-flyoffs(expi);
  targets(segi) = flyi;
  expis(segi) = expi;
  
  % downweight this fly
  idx = flyidx == fly;
  w0mat(:,idx) = w0mat(:,idx) / downweight_samefly;
  
  % downweight this fly's sex
  idx = ismalemat == ismalecurr;
  w0mat(:,idx) = w0mat(:,idx) / downweight_samesex;
  
  % downweight behaviors selected
  w0mat = bsxfun(@rdivide,w0mat,(nboutsmat(:,i)+1));
  
  % don't allow overlapping bouts
  i00 = max(i0-nframesperseg,1)+flystarts(fly)-1;
  i10 = min(i0+2*nframesperseg,trx(fly).nframes)+flystarts(fly)-1;
  w0mat(:,i00:i10) = nan;
  
end

expdirs_chosen = expdirs(expis);