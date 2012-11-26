%compute dcentralheadang
function [data,units]=compute_dcentralheadang(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
%absdtailheadang=cell(1,numlarvae);
dcentralheadang=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
%     absdtailheadang{1,i}=real(acos(cos(tailheadang{1,i}(1:end-1)).*cos(tailheadang{1,i}(2:end))+sin(tailheadang{1,i}(1:end-1)).*sin(tailheadang{1,i}(2:end))))./trx(i).dt;
%     temp=tailheadang{1,i}-pi/2;
%     cosperp=sign(cos(temp(1:end-1)).*cos(tailheadang{1,i}(2:end))+sin(temp(1:end-1)).*sin(tailheadang{1,i}(2:end)));
%     temp2=bsxfun(@times,absdtailheadang{1,i},cosperp);
    dcentralheadang{1,i}=modrange(diff(trx(larva).centralheadang),-pi,pi)./trx(larva).dt;
    % KB: this angle distance didn't make sense to me, changed to one that
    % made sense to me
    %dcentralheadang{1,i}=(mod1(trx(larva).centralheadang(2:end)-trx(larva).centralheadang(1:end-1),pi)-pi)./trx(larva).dt;
end

units=parseunits('rad/s');
data=dcentralheadang;
