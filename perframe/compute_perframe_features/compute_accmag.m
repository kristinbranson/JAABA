% magnitude of acceleration of center
function [data,units] = compute_accmag(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);
  
  % change in center position
  dx = diff(trx(fly).x_mm,1,2);
  dy = diff(trx(fly).y_mm,1,2);
  
  % acceleration magnitude
  if trx(fly).nframes < 2,
    data{i} = [];
  elseif trx(fly).nframes == 2,
    data{i} = 0;
  else
    % speed from frame 1 to 2 minus speed from 2 to 3 / time from 2 to 3
    tmp = sqrt(diff(dx./trx(fly).dt,1,2).^2 + diff(dy./trx(fly).dt,1,2).^2)./trx(fly).dt(2:end);
    data{i} = [0,tmp];
  end

end
units = parseunits('mm/s/s');

