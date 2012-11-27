%compute dtailheadang
function [data,units]=compute_dtailheadang(trx,n)
larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
dtailheadang=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    % KB: this angle distance didn't make sense to me, changed to one that
    % made sense to me
    dtailheadang{1,i}=modrange(diff(trx(larva).tailheadang),-pi,pi)./trx(larva).dt;
    % dtailheadang{1,i}=(mod1(trx(larva).tailheadang(2:end)-trx(larva).tailheadang(1:end-1),pi)-pi)./trx(larva).dt;
end

units=parseunits('rad/s');
data=dtailheadang;
