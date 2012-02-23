% closest fly, based on dell2nose
function [data,units,mind] = compute_closestfly_ell2nose(trx,n,dosave_d)

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
  d = nan(nflies,trx(fly1).nframes);
  
  % use dcenter2nose and major axis length to compute upper and lower bounds on
  % dell2nose
  [mindupper,dlower] = dell2nose_bounds(trx,fly1,flies);
  
  
  for i2 = 1:nflies,
    fly2 = flies(i2);
    if i1 == i2,
      continue;
    end
    % only try for frames where lower bound is smaller than min upper bound
    idx1try = find(mindupper >= dlower(fly2,:));
    d(i2,:) = dell2nose_pair(trx,fly1,fly2,idx1try);
  end
  [mind{i1},closesti] = min(d,[],1);
  closestfly{i1} = flies(closesti);
  closestfly{i1}(isnan(mind{i1})) = nan;
end

% so that we don't compute dcenter twice
if dosave_d,
  data = mind; %#ok<NASGU>
  units = parseunits('mm'); %#ok<NASGU>
  filename = trx.GetPerFrameFile('dell2nose',n);
  save(filename,'data','units');
end

data = closestfly;
units = parseunits('unit');
