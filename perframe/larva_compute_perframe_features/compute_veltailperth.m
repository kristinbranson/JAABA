%compute tail speed in the direction perpendicular to tail-head direction
function [data,units]=compute_veltailperth(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
veltailperth=cell(1,numlarvae);
%tailheadangperp=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    %tailheadangperp{1,i}=trx(larva).tailheadang-pi/2;
    veltailperth{1,i}=trx(larva).velmagtail.*(cos(trx(larva).velangtail).*cos(trx(larva).tailheadang(1,1:end-1)+pi/2)+sin(trx(larva).velangtail).*sin(trx(larva).tailheadang(1,1:end-1)+pi/2));
end

units=parseunits('mm/s');
data=veltailperth;