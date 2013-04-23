function [data,datatypes,success,msg] = ReadChorePerLarvaDatFile(filename)

data = [];
msg = '';
success = false;
if ~exist(filename,'file'),
  msg = sprintf('Dat file %s does not exist',filename);
  return;
end
fid = fopen(filename,'r');
if fid < 0,
  msg = sprintf('Failed to open %s',filename);
  return;
end

% parsing info
datatypes = {'timestamp','s'
  'frame','frame'
  'objectid','unit'
  'objectnumber','unit'
  'goodnumber','unit'
  'persistence','s'
  'area_mm','mm^2'
  'speed','mm/s'
  'angularspeed','deg/s'
  'skeletonlength','mm'
  'width','mm'
  'midlinelength','mm'
  'kink','deg'
  'bias','unit'
  'pathlength_mwt','mm'
  'curvature','deg'
  'changeindirection','unit'
  'x_mm','mm'
  'y_mm','mm'
  'velx','mm/s'
  'vely','mm/s'
  'orientation','deg'
  'crab','mm/s'
  'flux','unit'
  'headangle','deg'
  'headvecx','unit'
  'headvecy','unit'
  'tailvecx','unit'
  'tailvecy','unit'
  'tailx','px'
  'taily','px'
  'hcdist','mm'
  'headx','px'
  'heady','px'};
IDIDX = find(strcmp(datatypes(:,1),'objectid'));
FRAMEIDX = find(strcmp(datatypes(:,1),'frame'));
ndatatypes = size(datatypes,1);
% don't store these in 
fnsskip = {'frame','objectid'};

[~,idxskip] = ismember(fnsskip,datatypes(:,1)');

datamat = nan(0,34);

while true,
  
  s = fgetl(fid);
  if ~ischar(s),
    break;
  end
  s = strtrim(s);
  if isempty(s),
    continue;
  end
  ss = regexp(s,'\s','split');
  d = str2double(ss);
  datamat(end+1,:) = d; %#ok<AGROW>
  
end

fclose(fid);

ids = datamat(:,IDIDX); %#ok<FNDSB>
[uniqueids,~,ididx] = unique(ids);
for i = 1:numel(uniqueids),
  datacurr = struct;
  datacurr.id = uniqueids(i);
  idx = find(ididx == i);
  % find index for each frame
  minf = min(datamat(idx,FRAMEIDX));
  maxf = max(datamat(idx,FRAMEIDX));
  [ism,fidx] = ismember(minf:maxf,datamat(idx,FRAMEIDX));
  datacurr.firstframe = minf+1;
  datacurr.endframe = maxf+1;
  datacurr.nframes = maxf-minf+1;
  datacurr.off = 1-minf;
  for j = 1:ndatatypes,
    fn = datatypes{j,1};
    if any(j == idxskip),
      continue;
    end
    datacurr.(fn) = nan(1,datacurr.nframes);
    datacurr.(fn)(ism) = datamat(idx(fidx),j);
  end  
  data = structappend(data,datacurr);
end

success = true;