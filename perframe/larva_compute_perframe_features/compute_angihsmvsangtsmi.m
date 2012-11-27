%compute ang between vector tailsminfl and vector inflheadsm

function [data,units]=compute_angihsmvsangtsmi(trx,n)


larvae = trx.exp2flies{n};
numlarvae = numel(larvae);
angihsmvsangtsmi=cell(1,numlarvae);
absangihsmvsangtsmi=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    absangihsmvsangtsmi{1,i}=real(acos(cos(trx(larva).tailsminflang).*cos(trx(larva).inflheadsmang)+sin(trx(larva).tailsminflang).*sin(trx(larva).inflheadsmang)));
    temp=trx(larva).tailsminflang-pi/2;
    cosperp=sign(cos(temp).*cos(trx(larva).inflheadsmang)+sin(temp).*sin(trx(larva).inflheadsmang));
    angihsmvsangtsmi{1,i}=bsxfun(@times,absangihsmvsangtsmi{1,i},cosperp);
end

units=parseunits('rad');
data=angihsmvsangtsmi;