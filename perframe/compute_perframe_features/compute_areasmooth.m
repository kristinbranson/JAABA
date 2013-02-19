% area of ellipse, with outliers removed
function [data,units] = compute_areasmooth(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);
  data{i} = SmoothAreaOutliers(trx(fly).area,...
    trx.perframe_params.areasmooth_filterorder,...
    trx.perframe_params.areasmooth_maxfreq,...
    trx.perframe_params.areasmooth_maxerrx);
end
units = parseunits('mm^2');

