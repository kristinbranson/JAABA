% minimum distance from nose of fly to ellipse of any other fly within 30
% degrees 
function [data,units] = compute_dnose2ell_anglerange(trx,n,anglerange,varargin)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);

for i1 = 1:nflies,
  fly1 = flies(i1);
  % access closestfly to ensure that dnose2ell_anglerange is computed
  compute_closestfly_nose2ell_anglerange(trx,n,anglerange,varargin{:});
  data{i1} = trx(fly1).dnose2ell_angle_min30to30;
end
units = parseunits('mm');
