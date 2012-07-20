function [data,units] = compute_dwing_angle_imbalance(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);

  imbalancer = trx(fly).wing_angler+trx(fly).wing_anglel;
  imbalancel = -imbalancer;
  dimbalancer = diff(imbalancer);
  dimbalancel = diff(imbalancel);
  data{i} = dimbalancel;
  idx = imbalancer(1:end-1) > imbalancel(1:end-1);
  data{i}(idx) = dimbalancer(idx);
  
  data{i} = data{i} ./ trx(fly).dt;
  
end
units = parseunits('rad/s');

