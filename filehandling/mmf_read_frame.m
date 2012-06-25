function [im,header,timestamp,bb,mu] = mmf_read_frame(header,framei)

fp = header.fid;

% read from the current location if no frame input
if nargin < 2,
  loc = ftell(fp);
  framei = find(header.frame2file==loc,1);
  if isempty(framei),
    error('Could not find file location %d in frame index',loc);
  end
else
  fseek(fp,header.frame2file(framei),'bof');
end

[mu,header] = mmf_read_mean(header,'framei',framei);

im = mu;
%im_idcode = fread(fp,1,'ulong'); %#ok<NASGU>
%im_headersize = fread(fp,1,'int');
%im_depth = fread(fp,1,'int');
%im_nchannels = fread(fp,1,'int');
%fseek(fp,8,'cof');
%im_numblocks = fread(fp,1,'int');
im_numblocks = header.numblocks(framei);
    
% seek to the end of the header
%fseek(fp,header.frame2file(framei)+im_headersize,'bof');
fseek(fp,header.frame2file(framei),'bof');

bb = nan(im_numblocks,4);
isbaddata = false(1,im_numblocks);
for j = 1:im_numblocks,
  loc = ftell(fp);
  block_rect = fread(fp,4,'int');
  if any(block_rect([1,2]) < 0) || ...
      any(block_rect([1,2]) >= [header.height;header.width]) || ...
      any(block_rect([3,4]) <= 0) || ...
      any(block_rect([1,2])+block_rect([3,4]) >= [header.height;header.width])
    isbaddata(j) = true;
    warning('Bad block location %s read for frame %d, block %d, file loc %d, skipping',mat2str(block_rect),framei,j,loc);
    continue;
  end
      
  bb(j,[2,1,4,3]) = block_rect;
  % then interlaced row ordered image data
  block_data = fread(fp,block_rect(3)*block_rect(4),'*uint8');
  block_data = reshape(block_data,[block_rect(3),block_rect(4)]);
  im(block_rect(1)+1:block_rect(1)+block_rect(3),block_rect(2)+1:block_rect(2)+block_rect(4)) = block_data;
end
bb(isbaddata,:) = [];

timestamp = nan;