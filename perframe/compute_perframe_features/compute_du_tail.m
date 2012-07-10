% forward motion of tail
function [data,units] = compute_du_tail(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);
  
  if trx(fly).nframes < 2,
    data{i} = [];
  else

    % location of tail
    tailx = trx(fly).x_mm + 2*cos(-trx(fly).theta).*trx(fly).a_mm;
    taily = trx(fly).y_mm + 2*sin(-trx(fly).theta).*trx(fly).a_mm;

    % change in tail location
    dx = diff(tailx,1,2);
    dy = diff(taily,1,2);

    % forward motion of tail
    data{i} = (dx.*cos(trx(fly).theta_mm(1:end-1)) + dy.*sin(trx(fly).theta_mm(1:end-1)))./trx(fly).dt;

  end

end
units = parseunits('mm/s');

