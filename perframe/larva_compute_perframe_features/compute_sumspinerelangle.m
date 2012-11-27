%compute sumspinerelangle

function [data,units]=compute_sumspinerelangle(trx,n)

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
spineangles=cell(1,numlarvae);
spinerelangle=cell(1,numlarvae);
sumspinerelangle=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    spineangles{1,i}=bsxfun(@atan2,trx(larva).yspine_mm(1:end-1,:)-trx(larva).yspine_mm(2:end,:),trx(larva).xspine_mm(1:end-1,:)-trx(larva).xspine_mm(2:end,:));
    % spinerelangle{1,i}=bsxfun(@atan2,sin(spineangles{1,i}(2:end,:))-sin(spineangles{1,i}(1:end-1,:)),cos(spineangles{1,i}(2:end,:))-cos(spineangles{1,i}(1:end-1,:)));
    % KB: this method of measuring the difference between angles didn't
    % make sense to me. 
    spinerelangle{1,i} = modrange(-diff(spineangles{i},1,1),-pi,pi);
    sumspinerelangle{1,i}=sum(spinerelangle{1,i},1);
end

units=parseunits('rad');
data=sumspinerelangle;