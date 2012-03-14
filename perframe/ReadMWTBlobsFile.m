function [blobs,timestamps] = ReadMWTBlobsFile(filename,varargin)

[dotransposeimage] = myparse(varargin,'dotransposeimage',false);

if ~exist(filename,'file'),
  error('Blobs file %s does not exist',filename);
end
fid = fopen(filename,'r');
blobs = [];

while true,
  s = fgetl(fid);
  if ~ischar(s),
    break;
  end
  s = strtrim(s);
  if isempty(s),
    continue;
  end
  if s(1) ~= '%',
    error('Start of new blob should be a %, but was %s',s(1));
  end
  id = str2double(s(3:end));
  if isempty(id) || isnan(id),
    warning('Error reading blob id from string >%s<, skipping',s);
    continue;
  end
  blobcurr = ReadMWTBlobFile(fid,'dotransposeimage',dotransposeimage);
  blobcurr.id = id;
  blobs = structappend(blobs,blobcurr);
  
end

fclose(fid);

% make all the timestamps agree with each other
fnsmerge = {'x','y','a','b','theta','area','width','length'};
fnscopy = {};
[blobs,timestamps] = MergeMWTData(blobs,fnsmerge,fnscopy);