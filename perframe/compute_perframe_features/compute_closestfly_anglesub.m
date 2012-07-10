% closest fly, based on anglesub
function [data,units,maxanglesub] = compute_closestfly_anglesub(trx,n,dosave_d)

if nargin < 3,
  dosave_d = true;
end

flies = trx.exp2flies{n};
nflies = numel(flies);
closestfly = cell(1,nflies);
maxanglesub = cell(1,nflies);

for i1 = 1:nflies,
  fly1 = flies(i1);
  fprintf('fly1 = %d\n',fly1);
  anglesub = nan(nflies,trx(fly1).nframes);
  for i2 = 1:nflies,
    fly2 = flies(i2);
    if i1 == i2,
      continue;
    end
    anglesub(i2,:) = anglesub_pair(trx,fly1,fly2);
  end
  [maxanglesub{i1},closesti] = max(anglesub,[],1);
  closestfly{i1} = flies(closesti);
  closestfly{i1}(isnan(maxanglesub{i1})) = nan;
end

% so that we don't compute dcenter twice
if dosave_d,
  data = maxanglesub; %#ok<NASGU>
  units = parseunits('rad'); %#ok<NASGU>
  filename = trx.GetPerFrameFile('anglesub',n);
  try
    save(filename,'data','units');
  catch ME
    warning('Could not save anglesub data to %s: %s',filename,getReport(ME));
  end         
end

data = closestfly;
units = parseunits('unit');
