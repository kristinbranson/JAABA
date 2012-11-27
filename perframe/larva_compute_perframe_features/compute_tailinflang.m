%compute tail-inflection direction and correct trx if necessary
function [data,units]=compute_tailinflang(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
tailinflang=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
     tailinflang{1,i}=bsxfun(@atan2,trx(larva).yinflection_mm-trx(larva).ytail_mm,trx(larva).xinflection_mm-trx(larva).xtail_mm);
   
end
units=parseunits('rad');
data=tailinflang;
