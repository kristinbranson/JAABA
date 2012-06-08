function [scoreidx,labelidx] = LoadScoresFromFile(trx,scorefilestr,expis)

if nargin < 2,
  expis = 1:trx.nexpdirs;
end

labelidx = cell(1,trx.nflies);
scoreidx = cell(1,trx.nflies);
for expii = 1:numel(expis),
  expi = expis(expii);
  flies = trx.exp2flies{expi};
  scores_curr = load(fullfile(trx.expdirs{expi},scorefilestr));
  for flyi = 1:numel(flies),
    fly = flies(flyi);

    T0 = trx.firstframes(fly);
    T1 = trx.endframes(fly);

    n = T1-T0+1;
    off = 1 - T0;
    labelidx{fly} = false(1,n);
    scoreidx{fly} = nan(1,n);

    if fly > numel(scores_curr.allScores.scores),
      continue;
    end
    if ~isfield(scores_curr.allScores,'t0s'),
      labelidx{fly} = scores_curr.allScores.scores{fly} > 0;
    else
      for j = 1:numel(scores_curr.allScores.t0s{flyi}),
        t0 = scores_curr.allScores.t0s{flyi}(j);
        t1 = scores_curr.allScores.t1s{flyi}(j);
        if t0>T1 || t1<T0,
          warning('Labels out of bounds for exp %s, fly %d, label file %s',trx.expdirs{expi},fly,scorefilestr);
          continue;
        end
        t0 = max(T0,t0);
        t1 = min(T1,t1);
        labelidx{fly}(t0+off:t1-1+off) = true;
      end
    end
    t0 = min(T1,max(T0,scores_curr.allScores.tStart(fly)));
    t1 = min(T1,max(T0,scores_curr.allScores.tEnd(fly)));
    scoreidx{fly}(t0+off:t1+off) = scores_curr.allScores.scores{fly}(t0:t1);
  end
end

