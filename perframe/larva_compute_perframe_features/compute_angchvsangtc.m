%compute ang between vector tailcentral and vector centralhead

function [data,units]=compute_angchvsangtc(trx,n)
larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
angchvsangtc=cell(1,numlarvae);
absangchvsangtc=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    absangchvsangtc{1,i}=real(acos(cos(trx(larva).tailcentralang).*cos(trx(larva).centralheadang)+sin(trx(larva).tailcentralang).*sin(trx(larva).centralheadang)));
    temp=trx(larva).tailcentralang-pi/2;
    cosperp=sign(cos(temp).*cos(trx(larva).centralheadang)+sin(temp).*sin(trx(larva).centralheadang));
    angchvsangtc{1,i}=bsxfun(@times,absangchvsangtc{1,i},cosperp);
end

units=parseunits('rad');
data=angchvsangtc;
