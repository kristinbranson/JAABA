%compute ang between vector tailinfl and vector inflhead

function [data,units]=compute_angihvsangti(trx,n)


larvae = trx.exp2flies{n};
numlarvae = numel(larvae);
angihvsangti=cell(1,numlarvae);
absangihvsangti=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    absangihvsangti{1,i}=real(acos(cos(trx(larva).tailinflang).*cos(trx(larva).inflheadang)+sin(trx(larva).tailinflang).*sin(trx(larva).inflheadang)));
    temp=trx(larva).tailinflang-pi/2;
    cosperp=sign(cos(temp).*cos(trx(larva).inflheadang)+sin(temp).*sin(trx(larva).inflheadang));
    angihvsangti{1,i}=bsxfun(@times,absangihvsangti{1,i},cosperp);
end

units=parseunits('rad');
data=angihvsangti;
