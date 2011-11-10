% velocity direction
function [data,units] = compute_phi(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);
  
  if trx(fly).nframes < 2,
    data{i} = trx(fly).theta_mm;
  else
    % change in center position
    dy1 = [trx(fly).y_mm(2)-trx(fly).y_mm(1),(trx(fly).y_mm(3:end)-trx(fly).y_mm(1:end-2))/2,trx(fly).y_mm(end)-trx(fly).y_mm(end-1)];
    dx1 = [trx(fly).x_mm(2)-trx(fly).x_mm(1),(trx(fly).x_mm(3:end)-trx(fly).x_mm(1:end-2))/2,trx(fly).x_mm(end)-trx(fly).x_mm(end-1)];
    badidx = dy1 == 0 & dx1 == 0;
    data{i} = atan2(dy1,dx1);
    data{i}(badidx) = trx(fly).theta_mm(badidx);
  end
end
units = parseunits('rad');

