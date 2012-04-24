function labelidx = LoadLabelsFromFile(trx,labelfilestr,expis)

if nargin < 2,
  expis = 1:trx.nexpdirs;
end

labelidx = cell(1,trx.nflies);
for expii = 1:numel(expis),
  expi = expis(expii);
  flies = trx.exp2flies{expi};
  labels_curr = load(fullfile(trx.expdirs{expi},labelfilestr));
  for flyi = 1:numel(flies),
    fly = flies(flyi);

    T0 = trx.firstframes(fly);
    T1 = trx.endframes(fly);

    n = T1-T0+1;
    off = 1 - T0;
    labelidx{fly} = false(1,n);

    if fly > numel(labels_curr.t0s),
      continue;
    end
    for j = 1:numel(labels_curr.t0s{flyi}),
      t0 = labels_curr.t0s{labels_curr.flies(flyi)}(j)-labels_curr.off(flyi);
      t1 = labels_curr.t1s{labels_curr.flies(flyi)}(j)-labels_curr.off(flyi);
      if t0>T1 || t1<T0,
        warning('Labels out of bounds for exp %s, fly %d, label file %s',trx.expdirs{expi},fly,labelfilestr);
        continue;
      end
      t0 = max(T0,t0);
      t1 = min(T1,t1);
      labelidx{fly}(t0+off:t1-1+off) = true;
    end
  end
end

