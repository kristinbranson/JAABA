%compute tail-head direction and correct trx if necessary
function [data,units]=compute_tailheadang(trx,n)
larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
headtailang=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    nspinepts = size(trx(larva).xspine_mm,1);
    xspine=trx(larva).xspine_mm([1,nspinepts],:);
    yspine=trx(larva).yspine_mm([1,nspinepts],:);
    headtailang{1,i}=atan2(yspine(1,:)-yspine(2,:),xspine(1,:)-xspine(2,:));
end
units=parseunits('rad');
data=headtailang;
