% minimum distance from nose of fly to ellipse of any other fly within 30
% degrees 
function [data,units] = compute_dnose2ell_anglerange(trx,n,anglerange,varargin)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);

% access closestfly to ensure that dnose2ell_anglerange is computed
compute_closestfly_nose2ell_anglerange(trx,n,anglerange,varargin{:});

anglerange_deg = round(anglerange*180/pi);
fn = sprintf('dnose2ell_angle_%dto%d',anglerange_deg(1),anglerange_deg(2));
fn = strrep(fn,'-','min');
for i1 = 1:nflies,
  fly1 = flies(i1);
  data{i1} = trx(fly1).(fn);
end
units = parseunits('mm');
