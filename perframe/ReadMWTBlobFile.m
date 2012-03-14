function blob = ReadMWTBlobFile(fid,varargin)

dotransposeimage = myparse(varargin,'dotransposeimage',false);

didopen = false;
if ischar(fid),
  didopen = true;
  fid = fopen(fid,'r');
end

blob = struct('timestamp',[],'x',[],'y',[],...
  'area',[],'a',[],'b',[],'theta',[],'length',[],'width',[]);
frame = [];
while true,
  lastloc = ftell(fid);
  s = fgetl(fid);
  if ~ischar(s),
    break;
  end
  if s(1) == '%',
    fseek(fid,lastloc,'bof');
    break;
  end
  [datacurr,n,errmsg] = sscanf(s,'%f',10);
  if n < 10,
    warning('Error reading line >%s<: not enough fields: %s',s,errmsg);
  end
  frame(end+1) = datacurr(1); %#ok<AGROW>
  blob.timestamp(end+1) = datacurr(2);
  if dotransposeimage,
    blob.x(end+1) = datacurr(3)+1;
    blob.y(end+1) = datacurr(4)+1;
  else
    blob.x(end+1) = datacurr(4)+1;
    blob.y(end+1) = datacurr(3)+1;
  end
  blob.area(end+1) = datacurr(5);
  if dotransposeimage,
    u = datacurr(6);
    v = datacurr(7);
  else
    u = datacurr(7);
    v = datacurr(6);
  end
  a = sqrt(u.^2 + v.^2);
  theta = atan2(v,u);
  blob.a(end+1) = a;
  blob.theta(end+1) = theta;
  blob.b(end+1) = datacurr(8);
  blob.length(end+1) = datacurr(9);
  blob.width(end+1) = datacurr(10);
end
blob.firstframe = frame(1);
blob.endframe = frame(end);
blob.off = 1-blob.firstframe;
blob.nframes = blob.endframe-blob.firstframe+1;

% % make sure data is contiguous
% fns = {'x','y','area','a','b','theta','length','width','timestamp'};
% idx = frame + blob.off;
% for i = 1:numel(fns),
%   fn = fns{i};
%   tmp = nan(1,blob.nframes);
%   tmp(idx) = blob.(fn);
%   blob.(fn) = tmp;
% end

if didopen,
  fclose(fid);
end