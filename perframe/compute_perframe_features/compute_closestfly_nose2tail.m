% closest fly, based on dnose2tail
function [data,units,mind] = compute_closestfly_nose2tail(trx,n,dosave_d)

if nargin < 3,
  dosave_d = true;
end

flies = trx.exp2flies{n};
nflies = numel(flies);
closestfly = cell(1,nflies);
mind = cell(1,nflies);

for i1 = 1:nflies,
  fly1 = flies(i1);
  fprintf('fly1 = %d\n',fly1);
  flies2 = flies(trx.roi(fly1)==trx.roi(flies));
  d = nan(numel(flies2),trx(fly1).nframes);  
  
  for i2 = 1:numel(flies2),
    fly2 = flies2(i2);
    if fly1 == fly2,
      continue;
    end
    d(i2,:) = dnose2tail_pair(trx,fly1,fly2);
  end
  [mind{i1},closesti] = min(d,[],1);
  closestfly{i1} = flies2(closesti);
  closestfly{i1}(isnan(mind{i1})) = nan;
end

% so that we don't compute dcenter twice
if dosave_d,
  data = mind; %#ok<NASGU>
  units = parseunits('mm'); %#ok<NASGU>
  filename = trx.GetPerFrameFile('dnose2tail',n);
  save(filename,'data','units');
end

data = closestfly;
units = parseunits('unit');
