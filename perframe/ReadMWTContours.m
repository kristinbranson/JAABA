function [trx,timestamps] = ReadMWTContours(contournames,ids,varargin)

[dotransposeimage] = myparse(varargin,'dotransposeimage',false);

ncontours = numel(contournames);
% read everthing in from the contour files
trx = [];
for i = 1:ncontours,

  trk = struct('id',ids(i),'xcontour',{{}},'ycontour',{{}});
  timestamps = [];
  fid = fopen(contournames{i},'r');
  j = 1;
  while true,
    s = fgetl(fid);
    if ~ischar(s),
      break;
    end
    d = sscanf(s,'%f');
    if numel(d) < 3,
      warning('Contour line %s length < 3, skipping this line',s);
      continue;
    end
    if mod(numel(d),2) ~= 1,
      warning('Contour line %s length is not odd, skipping this line',s);
      continue;
    end
    d = d';
    timestamp = d(1);
    x = d(2:2:end-1);
    y = d(3:2:end);
    if dotransposeimage,
      trk.xcontour{j} = x+1;
      trk.ycontour{j} = y+1;
    else
      trk.xcontour{j} = y+1;
      trk.ycontour{j} = x+1;
    end
    timestamps(j) = timestamp; %#ok<*AGROW>
    j = j + 1;
  end
  fclose(fid);
  trk.timestamp = timestamps;
  trx = structappend(trx,trk);
end

% make all the timestamps agree with each other
fnsmerge = {'xcontour','ycontour'};
fnscopy = {};
[trx,timestamps] = MergeMWTData(trx,fnsmerge,fnscopy);


