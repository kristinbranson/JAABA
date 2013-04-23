function [success,msg,fps] = ReadFPS_MWT(varargin) 

success = false;
msg = ''; %#ok<NASGU>

MINDT = .001;

[fps,blobsfile,spinefile,datfiles,matfiles] = ...
  myparse(varargin,...
  'fps',30,...
  'blobsfile','','spinefile','','datfiles',{},'matfiles',{});  %#ok<NASGU,ASGLU>

%% check that trx exist
if isempty(blobsfile),
  msg = 'Input blobs file not yet set';
  return;
end
if ischar(blobsfile),
  blobsfile = {blobsfile};
end
if ~exist(blobsfile{1},'file'),
  msg = sprintf('Blobs file %s does not exist',blobsfile{1});
  return;
end

%% read the trx file
trx = ReadMWTBlobFile(blobsfile{1});

%% compute fps

timestamps = unique([trx.timestamps]);
dts = diff(timestamps);
% sanity check to look for rounding errors
if any(dts < MINDT),
  msg = sprintf('Rounding error check failed: there are timestamps < %f apart',MINDT);
  return;
end

fps = 1/nanmedian(dts);

msg = sprintf('Read fps from blobs file %s',blobsfile{1});
success = true;
