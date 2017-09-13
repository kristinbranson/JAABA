function chunk = readchunk(fid)
%READCHUNK read riff file chunk
%   CHUNK = READCHUNK(FID) reads a four character chunk ID and a 32-bit
%   integer chunk size into the 'ckid' and 'cksize' fields of CHUNK, from
%   the RIFF file associated with FID.

%   Copyright 1984-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.5 $  $Date: 2011/05/13 17:34:38 $

[id, count] = fread(fid,4,'uchar');
chunk.ckid = [char(id)]';
if (count ~= 4)
  error(message('MATLAB:audiovideo:readchunk:badChunkRead'));
end

[chunk.cksize, count] = fread(fid,1,'uint32');
if (count ~= 1)
  error(message('MATLAB:audiovideo:readchunk:badChunkRead'));
end
return;
