%compute dcentralheadmag
function [data,units]=compute_dcentralheadmag(trx,n)

% inputfilename=[outputfolder,'centralheadmag.mat'];
% if ~exist(inputfilename,'file')
%     [trx]=compute_centralheadmag(trx,outputfolder);
% end
% load([outputfolder,'centralheadmag.mat'], 'data')
% centralheadmag=data;
% 

larvae=trx.exp2flies{n};
numlarvae=numel(larvae);
dcentralheadmag=cell(1,numlarvae);
for i=1:numlarvae
    larva=larvae(i);
    dcentralheadmag{1,i}=(trx(larva).centralheadmag(2:end)-trx(larva).centralheadmag(1:end-1))./trx(larva).dt;
end

units=parseunits('mm/s');
data=dcentralheadmag;
% filename=[outputfolder, 'dcentralheadmag.mat'];
% save(filename, 'data', 'units')