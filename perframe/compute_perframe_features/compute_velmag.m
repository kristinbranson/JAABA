% magnitude of velocity of center of rotation
function [data,units] = compute_velmag(trx,n)

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

    % magnitude of velocity vector
    data{i} = sqrt(dx_cor.^2 + dy_cor.^2)./trx(fly).dt;
    badidx = isnan(dx_cor);
    data{i}(badidx) = trx(fly).velmag_ctr(badidx);

  end

end
units = parseunits('mm/s');

