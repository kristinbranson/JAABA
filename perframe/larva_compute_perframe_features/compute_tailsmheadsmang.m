%compute tailsm-headsm direction and correct trx if necessary
function [data,units]=compute_tailsmheadsmang(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
headsmtailsmang=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    xspine=trx(larva).xspine_mm([1,11],:);
    yspine=trx(larva).yspine_mm([1,11],:);
    headsmtailsmang{1,i}=atan2(yspine(1,:)-yspine(2,:),xspine(1,:)-xspine(2,:));

end
units=parseunits('rad');
data=headsmtailsmang;
