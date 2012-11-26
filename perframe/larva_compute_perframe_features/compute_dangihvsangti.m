%compute dangihvsangti
function [data,units]=compute_dangihvsangti(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
%absdtailheadang=cell(1,numlarvae);
dangihvsangti=cell(1,numlarvae);
for i=1:numlarvae
  larva=larvae(i);
  
  dangihvsangti{1,i}=modrange(diff(trx(larva).angihvsangti),-pi,pi)./trx(larva).dt;
  % KB: this angle distance didn't make sense to me, changed to one that
  % made sense to me
  % dangihvsangti{1,i}=(mod1(trx(larva).angihvsangti(2:end)-trx(larva).angihvsangti(1:end-1),pi)-pi)./trx(larva).dt;
end

units=parseunits('rad/s');
data=dangihvsangti;
