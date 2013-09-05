% distance to region of interest
function [data,units] = compute_angle2roi2(trx,n,idx)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);

  x = trx(fly).x_mm;
  y = trx(fly).y_mm;
  ROIdata=trx.roi2{n}.data{fly}(idx,:);
  xc = ROIdata(1);
  yc = ROIdata(2);
  direction = atan2((yc-y),(xc-x));

  data{i} = modrange(direction - trx(fly).theta_mm,-pi,pi);
end

units = parseunits('rad');