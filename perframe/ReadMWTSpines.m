function [trx,timestamps] = ReadMWTSpines(spinenames,ids,varargin)

[dotransposeimage] = myparse(varargin,'dotransposeimage',false);

nspines = numel(spinenames);
% read everthing in from the spine files
trx = [];
isfirst = true;
npts = 0;
for i = 1:nspines,

  trk = struct('id',ids(i),'xspine',nan(npts,0),'yspine',nan(npts,0));
  timestamps = [];
  fid = fopen(spinenames{i},'r');
  j = 1;
  while true,
    s = fgetl(fid);
    if ~ischar(s),
      break;
    end
    d = sscanf(s,'%f');
    if numel(d) < 3,
      warning('Spine line %s length < 3, skipping this line',s);
      continue;
    end
    if mod(numel(d),2) ~= 1,
      warning('Spine line %s length is not odd, skipping this line',s);
      continue;
    end
    d = d';
    timestamp = d(1);
    x = d(2:2:end-1);
    y = d(3:2:end);
    if isfirst,
      npts = numel(x);
      trk.xspine = nan(npts,0);
      trk.yspine = nan(npts,0);
      isfirst = false;
    else
      if numel(x) ~= npts,
        warning('Spine at line >%s<, id %d must have %d points',s,i,npts);
      end
    end
    if dotransposeimage,
      trk.xspine(:,j) = x+1;
      trk.yspine(:,j) = y+1;
    else
      trk.xspine(:,j) = y+1;
      trk.yspine(:,j) = x+1;
    end
    timestamps(j) = timestamp; %#ok<*AGROW>
    j = j + 1;
  end
  fclose(fid);
  trk.timestamp = timestamps;
  trx = structappend(trx,trk);
end

% make all the timestamps agree with each other
fnsmerge = {'xspine','yspine'};
fnscopy = {};
[trx,timestamps] = MergeMWTData(trx,fnsmerge,fnscopy);


