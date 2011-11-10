% absolute azimuthal angle to from fly to closest fly according to type
function [data,units] = compute_absanglefrom1to2(trx,n,type)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);

for i1 = 1:nflies,
  fly1 = flies(i1);
  data{i1} = abs(trx(fly1).(['anglefrom1to2_',type]));
end

units = parseunits('rad');