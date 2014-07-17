% angle to the closest roi
function [data,units] = compute_angle2closestroi2(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);

for i = 1:nflies,
  fly = flies(i);
  
  % roi closest to fly
  closestroi = trx(fly).closestroi2;
  
  % position of fly
  xnose_mm1 = trx(fly).xnose_mm;
  ynose_mm1 = trx(fly).ynose_mm;
  theta_mm1 = trx(fly).theta_mm;

  ROIdata=trx.roi2{n}.data{fly};
  
  % loop over all rois
  for j = 1:size(ROIdata,1),
    
    % frames where this roi is closest
    idx = find(closestroi == j);
    if isempty(idx), continue; end
    
    % angle to roi
    dx2 = ROIdata(j,1)-xnose_mm1(idx);
    dy2 = ROIdata(j,2)-ynose_mm1(idx);
    theta2 = atan2(dy2,dx2);
    
    % angle relative to fly's orientation
    data{i}(idx) = modrange(theta2 - theta_mm1(idx),-pi,pi);

  end
end

units = parseunits('rad');