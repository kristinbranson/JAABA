% forward motion of center of rotation
function [data,units] = compute_du_cor(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);
  
  if trx(fly).nframes < 2,
    data{i} = [];
  else
    % center of rotation
    [x_cor_curr,y_cor_curr,x_cor_next,y_cor_next] = rfrac2center(trx,fly,[trx(fly).corfrac_maj;trx(fly).corfrac_min]);
    % change in center of rotation
    dx_cor = x_cor_next - x_cor_curr;
    dy_cor = y_cor_next - y_cor_curr;      
    % forward motion of center of rotation
    data{i} = (dx_cor.*cos(trx(fly).theta_mm(1:end-1)) + dy_cor.*sin(trx(fly).theta_mm(1:end-1)))./trx(fly).dt;

  end

end
units = parseunits('mm/s');

