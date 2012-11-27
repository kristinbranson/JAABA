%compute velang_tail
function [data,units]=compute_velangtail(trx,n)
larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
velangtail=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    velangtail{1,i}=bsxfun(@atan2,trx(larva).dytail_mm,trx(larva).dxtail_mm);
end

units=parseunits('rad');
data=velangtail;