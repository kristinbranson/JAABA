% closest fly, based on dcenter
function [data,units,mind] = compute_closestfly_center(trx,n,dosave_d)

if nargin < 3,
  dosave_d = true;
end

flies = trx.exp2flies{n};
nflies = numel(flies);
closestfly = cell(1,nflies);
mind = cell(1,nflies);

for i1 = 1:nflies,
  fly1 = flies(i1);
  dcenter = nan(nflies,trx(fly1).nframes);
  for i2 = 1:nflies,
    fly2 = flies(i2);
    if i1 == i2,
      continue;
    end
    dcenter(i2,:) = dcenter_pair(trx,fly1,fly2);
  end
  [mind{i1},closesti] = min(dcenter,[],1);
  closestfly{i1} = flies(closesti);
  closestfly{i1}(isnan(mind{i1})) = nan;
end

% so that we don't compute dcenter twice
if dosave_d,
  data = mind; %#ok<NASGU>
  units = parseunits('mm'); %#ok<NASGU>
  filename = trx.GetPerFrameFile('dcenter',n);
  save(filename,'data','units');
end

data = closestfly;
units = parseunits('unit');
