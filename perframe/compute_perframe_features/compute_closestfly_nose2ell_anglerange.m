function [data,units,mind] = compute_closestfly_nose2ell_anglerange(trx,n,anglerange,dosave_d)

if nargin < 4,
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
  for i2 = 1:nflies,
    if i1 == i2,
      continue;
    end
    fly2 = flies(i2);
    d(i2,:) = dnose2ell_anglerange_pair(trx,fly1,fly2,anglerange);
  end
  [mind{i1},closesti] = min(d,[],1);
  closestfly{i1} = flies(closesti);
  closestfly{i1}(isnan(mind{i1})) = nan;
  mind{i1}(isnan(mind{i1})) = inf;
end

anglerange_deg = round(anglerange*180/pi);

% so that we don't compute dcenter twice
if dosave_d,
  data = mind; %#ok<NASGU>
  units = parseunits('mm'); %#ok<NASGU>
  tmp = sprintf('dnose2ell_angle_%dto%d',anglerange_deg(1),anglerange_deg(2));
  tmp = strrep(tmp,'-','min');
  filename = trx.GetPerFrameFile(tmp,n);
  save(filename,'data','units');
end

data = closestfly;
units = parseunits('unit');