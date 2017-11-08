function [msg msgID] = skipchunk(fid,chunk)
%SKIPCHUNK skip chunk in AVI
%   [MSG MSGID] = SKIPCHUNK(FID,CHUNK) skips CHUNK.cksize bytes in the AVI file
%   FID.  MSG contains an error message string if the skip fails, otherwise
%   it is an empty string.

%   Copyright 1984-2007 The MathWorks, Inc.
%   $Revision: 1.1.6.6 $  $Date: 2011/07/20 00:00:24 $

msg = '';
msgID = '';
% Determine if pad byte is necessary
if ( rem(chunk.cksize,2) ) % If the chunk size is odd, there is a pad byte
  pad = 1;
else 
  pad = 0;
end

% Skip the chunk
status = fseek(fid,chunk.cksize + pad,0);
if ( status == -1 )
  msgID = 'MATLAB:audiovideo:skipChunk:incorrectChunkSize';
  msg = getString(message(msgID));
end
return;
