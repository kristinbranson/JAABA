%compute tail speed in the direction perpendicular to the central-head direction
function [data,units]=compute_veltailperch(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
veltailperch=cell(1,numlarvae);
%centralheadangperp=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    %centralheadangperp{1,i}=trx(larva).centralheadang-pi/2;
    veltailperch{1,i}=trx(larva).velmagtail.*(cos(trx(larva).velangtail).*cos(trx(larva).centralheadang(1,1:end-1)+pi/2)+sin(trx(larva).velangtail).*sin(trx(larva).centralheadang(1,1:end-1)+pi/2));
end

units=parseunits('mm/s');
data=veltailperch;