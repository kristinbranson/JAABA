% distance to closest region of interest
function [data,units,mind] = compute_closestroi2(trx,n,dosave_d)

if nargin < 3,
  dosave_d = true;
end
 

flies = trx.exp2flies{n};
nflies = numel(flies);
closestfly = cell(1,nflies);
mind = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);
  ROIdata=trx.roi2{n}.data{fly};
  mind{i} = inf(1,trx(fly).nframes);
  closestfly{i} = nan(1,trx(fly).nframes);
  for j = 1:size(ROIdata,1),
    fn = sprintf('dist2roi2_%d',j);
    dcurr = trx(fly).(fn);
    idx = dcurr < mind{i};
    mind{i}(idx) = dcurr(idx);
    closestfly{i}(idx) = j;
  end
end 
  
% so that we don't compute dcenter twice
if dosave_d,
  data = mind; %#ok<NASGU>
  units = parseunits('mm'); %#ok<NASGU>
  filename = trx.GetPerFrameFile('mindist2roi2',n);
  save(filename,'data','units');
end

data = closestfly;
units = parseunits('unit');
