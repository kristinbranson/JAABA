% distance to region of interest
function [data,units] = compute_dist2roi2(trx,n,idx)

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
  distance = sqrt((x-xc).^2+(y-yc).^2);
  switch isnan(ROIdata(4))
    case 1
      radius = ROIdata(3);
      distance = distance - radius;
    case 0
      width = ROIdata(3);
      height = ROIdata(4);
      direction = atan((y-yc)/(x-xc));
      if direction < atan(height/width)
        distance = distance - height/2/cos(direction);
      else
        distance = distance - width/2/cos(direction);
      end
  end

  data{i} = distance;
end

units = parseunits('mm');