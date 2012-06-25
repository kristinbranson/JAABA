function [scoreidx,labelidx] = LoadScoresFromFile(trx,scorefilestr,expi)

labelidx = cell(1,trx.nfliespermovie(expi));
scoreidx = cell(1,trx.nfliespermovie(expi));

flies = trx.exp2flies{expi};
scores_curr = load(fullfile(trx.expdirs{expi},scorefilestr));
for flyi = 1:numel(flies),
  fly = flies(flyi);

  T0 = trx.firstframes(fly);
  T1 = trx.endframes(fly);

  n = T1-T0+1;
  off = 1 - T0;
  labelidx{flyi} = false(1,n);
  scoreidx{flyi} = nan(1,n);

  if flyi > numel(scores_curr.allScores.scores),
    continue;
  end
  if ~isfield(scores_curr.allScores,'t0s'),
    labelidx{flyi} = scores_curr.allScores.scores{flyi} > 0;
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
      labelidx{flyi}(t0+off:t1-1+off) = true;
    end
  end
  t0 = min(T1,max(T0,scores_curr.allScores.tStart(flyi)));
  t1 = min(T1,max(T0,scores_curr.allScores.tEnd(flyi)));
  scoreidx{flyi}(t0+off:t1+off) = scores_curr.allScores.scores{flyi}(t0:t1);
end

