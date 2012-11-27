%compute dtailinflang
function [data,units]=compute_dtailinflang(trx,n)


larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
%absdtailheadang=cell(1,numlarvae);
dtailinflang=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    dtailinflang{1,i}=modrange(diff(trx(larva).tailinflang),-pi,pi)./trx(larva).dt;
    % KB: this angle distance didn't make sense to me, changed to one that
    % made sense to me
    % dtailinflang{1,i}=(mod1(trx(larva).tailinflang(2:end)-trx(larva).tailinflang(1:end-1),pi)-pi)./trx(larva).dt;
end

units=parseunits('rad/s');
data=dtailinflang;
