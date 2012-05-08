% rotation of nose around mean tail location
function [data,units] = compute_dtheta_tail(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);  
  
  % location of tail
  tailx = trx(fly).x_mm + 2*cos(-trx(fly).theta).*trx(fly).a_mm;
  taily = trx(fly).y_mm + 2*sin(-trx(fly).theta).*trx(fly).a_mm;
  % location of nose
  nosex = trx(fly).x_mm + 2*cos(trx(fly).theta).*trx(fly).a_mm;
  nosey = trx(fly).y_mm + 2*sin(trx(fly).theta).*trx(fly).a_mm;

  if trx(fly).nframes < 2,
    trx(fly).dtheta_tail = [];
  else
    meantailx = (tailx(1:end-1)+tailx(2:end))/2;
    meantaily = (taily(1:end-1)+taily(2:end))/2;
    anglenose1 = atan2(nosey(1:end-1)-meantaily,nosex(1:end-1)-meantailx);
    anglenose2 = atan2(nosey(2:end)-meantaily,nosex(2:end)-meantailx);
    data{i} = modrange(anglenose2-anglenose1,-pi,pi)./trx(fly).dt;
  end

end
units = parseunits('rad/s');

