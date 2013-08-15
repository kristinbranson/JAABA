%compute velang_tailhead
function [data,units]=compute_velangtailhead(trx,n)
larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
velangtailhead=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    velangtailhead{1,i}=bsxfun(@atan2,trx(larva).dytailhead_mm,trx(larva).dxtailhead_mm);
end

units=parseunits('rad');
data=velangtailhead;