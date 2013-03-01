% closest fly, based on dnose2ell
function [data,units,mind,angle] = compute_closestfly_nose2ell(trx,n,dosave_d)

if nargin < 3,
  dosave_d = true;
end

flies = trx.exp2flies{n};
nflies = numel(flies);
closestfly = cell(1,nflies);
mind = cell(1,nflies);
angle = cell(1,nflies);

for i1 = 1:nflies,
  fly1 = flies(i1);
  fprintf('target 1 = %d\n',fly1);

  flies2 = flies(trx.roi(fly1)==trx.roi(flies));
  
  % use dnose2center and major axis length to compute upper and lower bounds on
  % dnose2ell
  % fprintf('CHANGE THIS: THIS IS JUST FOR DEBUGGING NOSE2ELL\n');
  [mindupper,dlower] = dnose2ell_bounds(trx,fly1,flies2);
  %mindupper = zeros(1,trx(fly1).nframes);
  %dlower = zeros(nflies,trx(fly1).nframes);
  
  mind{i1} = inf(1,trx(fly1).nframes);
  closesti = ones(1,trx(fly1).nframes);
  
  
  %d = nan(nflies,trx(fly1).nframes);
  for i2 = 1:numel(flies2),
    fly2 = flies2(i2);
    if fly1 == fly2,
      continue;
    end
    % only try for frames where lower bound is smaller than min upper bound
    idx1try = find(mindupper >= dlower(i2,:));
    [dcurr,anglecurr] = dnose2ell_pair(trx,fly1,fly2,idx1try);
    idx = dcurr < mind{i1};
    mind{i1}(idx) = dcurr(idx);
    closesti(idx) = i2;
    angle{i1}(idx) = anglecurr(idx);
  end
  closestfly{i1} = flies2(closesti);
  closestfly{i1}(isnan(mind{i1})|isinf(mind{i1})) = nan;
end

% so that we don't compute dcenter twice
if dosave_d,
  data = mind; %#ok<NASGU>
  units = parseunits('mm'); %#ok<NASGU>
  filename = trx.GetPerFrameFile('dnose2ell',n);
  try
    save(filename,'data','units');
  catch ME,
    warning('Could not save file %s:\n%s',filename,getReport(ME));
  end
  data = angle; %#ok<NASGU>
  units = parseunits('rad'); %#ok<NASGU>
  filename = trx.GetPerFrameFile('angleonclosestfly',n);
  if exist(filename,'file'),
    delete(filename);
  else
    if isunix,
      [res,link] = unix(sprintf('readlink %s',filename));
      if ~res && ~isempty(link),
        warning('Deleting broken soft link from %s to %s.\n',filename,link);
        unix(sprintf('rm %s',filename));
      end
    end
  end
  save(filename,'data','units');
end

data = closestfly;
units = parseunits('unit');
