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

  % compute upper and lower bounds on distance
  
  % initialize
  mindupper = inf(1,trx(fly).nframes);
  dlower = inf(nflies,trx(fly).nframes);
  dnose2center = nan(nflies,trx(fly).nframes);
  
  for i2 = 1:nflies,
    if i1 == i2,
      continue;
    end
    fly2 = flies(i2);
    
    % compute dnose2center, angle from nose to center
    [dcurr,off10,off11,off20,off21,t0,t1,anglefrom1to2] = dnose2center_pair(trx,fly1,fly2);
        
    if off10 > off11,
      continue;
    end
    
    % frames within the angle range
    anglefrom1to2 = modrange(anglefrom1to2,anglerange(1),anglerange(1)+2*pi);
    idx = find(anglefrom1to2 >= anglerange(1) & anglefrom1to2 <= anglerange(2));
    idx1 = idx+off10-1;
    idx2 = idx+off20-1;
    
    mindupper(idx1) = min(mindupper(idx1),dcurr(idx1) + 2*trx(fly2).a_mm(idx2));
    dlower(i2,idx1) = dcurr(idx1) - 2*trx(fly2).a_mm(idx2);
    dnose2center(i2,idx1) = dcurr(idx1) - 2*trx(fly2).a_mm(idx2);
  end
  
  d = nan(nflies,trx(fly1).nframes);
  for i2 = 1:nflies,
    if i1 == i2,
      continue;
    end
    fly2 = flies(i2);
    % this will only be true for flies in angle range
    istry = find(dlower(i2,:) <= mindupper);
    d(i2,:) = dnose2ell_anglerange_pair(trx,fly1,fly2,anglerange,istry);
  end
  [mind{i1},closesti] = min(d,[],1);
  closestfly{i1} = flies(closesti);
  closestfly{i1}(isnan(mind{i1})) = nan;
  mind{i1}(isnan(mind{i1})) = 10000000;
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