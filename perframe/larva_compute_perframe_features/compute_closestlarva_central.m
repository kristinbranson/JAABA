% closest larva, based on dcentral
function [data,units,mind] = compute_closestlarva_central(trx,n,dosave_d)

if nargin < 3,
  dosave_d = true;
end

larvae = trx.exp2flies{n};
nlarvae = numel(larvae);
closestlarva = cell(1,nlarvae);
mind = cell(1,nlarvae);

for i1 = 1:nlarvae,
  larva1 = larvae(i1);
  larvae2 = larvae;
  dcentral = inf(numel(larvae2),trx(larva1).nframes);
  for i2 = 1:numel(larvae2),
    larva2 = larvae2(i2);
    if larva1 == larva2,
      continue;
    end
    dcentral(i2,:) = dcentral_pair(trx,larva1,larva2);
  end
  [mind{i1},closesti] = min(dcentral,[],1);
  closestlarva{i1} = larvae2(closesti);
  closestlarva{i1}(isnan(mind{i1})) = nan;
end

% so that we don't compute dcentral twice
if dosave_d,
  data = mind; %#ok<NASGU>
  units = parseunits('mm'); %#ok<NASGU>
  filename = trx.GetPerFrameFile('dcentral',n);
  save(filename,'data','units');
end

data = closestlarva;
units = parseunits('unit');
