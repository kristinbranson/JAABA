
function mywriteUFMFFrame(fig,im,timestamp,fid,i,bb)

FRAME_CHUNK = 1;

gui = get(fig,'UserData');

ncc = size(bb,1);

% get location of this frame
loc = ftell(fid);
% store in index
gui.bg.index(i).frame.loc(end+1) = loc;
% also store timestamp
gui.bg.index(i).frame.timestamp(end+1) = timestamp;

% write chunk type: 1
fwrite(fid,FRAME_CHUNK,'uchar');
% write timestamp: 8
fwrite(fid,timestamp,'double');
% write number of points: 2
fwrite(fid,ncc,'uint32');

dtype = class(im);
for j = 1:ncc,
  % images are sideways: swap x and y, width and height
  % bb(j,1) = xmin
  % bb(j,2) = ymin
  % bb(j,3) = width
  % bb(j,4) = height
  fwrite(fid,[bb(j,[2,1]),bb(j,[4,3])],'ushort');
  tmp = im( (bb(j,1)+1):(bb(j,1)+bb(j,3)),(bb(j,2)+1):(bb(j,2)+bb(j,4)),:);
  fwrite(fid,permute(tmp,[3,2,1]),dtype);
end

set(fig,'UserData',gui);

