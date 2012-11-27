%compute tail peed in the direction perpendicular to tail-central direction
function [data,units]=compute_veltailpertc(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
veltailpertc=cell(1,numlarvae);
%tailcentralangperp=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    %tailcentralangperp{1,i}=trx(larva).tailcentralang-pi/2;
    veltailpertc{1,i}=trx(larva).velmagtail.*(cos(trx(larva).velangtail).*cos(trx(larva).tailcentralang(1,1:end-1)+pi/2)+sin(trx(larva).velangtail).*sin(trx(larva).tailcentralang(1,1:end-1)+pi/2));
end

units=parseunits('mm/s');
data=veltailpertc;