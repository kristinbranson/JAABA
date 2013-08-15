%compute velang_tailcentral
function [data,units]=compute_velangtailcentral(trx,n)
larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
velangtailcentral=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    velangtailcentral{1,i}=bsxfun(@atan2,trx(larva).dytailcentral_mm,trx(larva).dxtailcentral_mm);
end

units=parseunits('rad');
data=velangtailcentral;