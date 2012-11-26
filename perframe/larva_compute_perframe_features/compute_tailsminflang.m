%compute tailsm-inflection direction and correct trx if necessary
function [data,units]=compute_tailsminflang(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
tailsminflang=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
     tailsminflang{1,i}=bsxfun(@atan2,trx(larva).yinflection_mm-trx(larva).ytailsm_mm,trx(larva).xinflection_mm-trx(larva).xtailsm_mm);
   
end
units=parseunits('rad');
data=tailsminflang;