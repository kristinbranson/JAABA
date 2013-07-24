function [scoreidx,labelidx] = LoadScoresFromFile(trx,scorefilestr,expi)

labelidx = cell(1,trx.nfliespermovie(expi));
scoreidx = cell(1,trx.nfliespermovie(expi));

flies = trx.exp2flies{expi};
scores_curr = load(fullfile(trx.expdirs{expi},scorefilestr));
for flyi = 1:numel(flies),
  fly = flies(flyi);

  T0 = trx.firstframes(fly);
  T1 = trx.endframes(fly);

%   n = T1-T0+1;
%   labelidx{flyi} = false(1,n);
%   scoreidx{flyi} = nan(1,n);

  if flyi > numel(scores_curr.allScores.scores),
    continue;
  end
  if isfield(scores_curr.allScores,'postprocessed'),
    labelidx{flyi} = scores_curr.allScores.postprocessed{flyi}(T0:T1);
  else
    labelidx{flyi} = scores_curr.allScores.scores{flyi}(T0:T1) > 0;
  end
  scoreidx{flyi} = scores_curr.allScores.scores{flyi}(T0:T1);
end