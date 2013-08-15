%compute velang_centralhead
function [data,units]=compute_velangcentralhead(trx,n)
larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
velangcentralhead=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    velangcentralhead{1,i}=bsxfun(@atan2,trx(larva).dycentralhead_mm,trx(larva).dxcentralhead_mm);
end

units=parseunits('rad');
data=velangcentralhead;