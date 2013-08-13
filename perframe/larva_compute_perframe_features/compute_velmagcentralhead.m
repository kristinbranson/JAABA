%compute velmag_centralhead
function [data,units]=compute_velmagcentralhead(trx,n)
larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
velmagcentralhead=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    velmagcentralhead{1,i}=bsxfun(@hypot,trx(larva).dxcentralhead_mm,trx(larva).dycentralhead_mm);
end

units=parseunits('rad');
data=velmagcentralhead;