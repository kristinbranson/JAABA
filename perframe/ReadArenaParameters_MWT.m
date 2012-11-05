function [success,msg,pxpermm] = ReadArenaParameters_MWT(varargin) 

success = false;
msg = ''; %#ok<NASGU>

datext = 'dat';


[pxpermm,blobsfile,spinefile,datfiles] = ...
  myparse(varargin,...
  'pxpermm',1,...
  'blobsfile','','spinefile','','datfiles',{});  %#ok<ASGLU>

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

%% look for x_mm or y_mm
xi = [];
yi = [];
for i = 1:numel(datfiles),
  [~,name] = myfileparts(datfiles{i});
  match = regexp(name,['\.(.+)\.',datext,'$'],'tokens','once');
  if isempty(match),
    continue;
  end
  fn = match{1};
  if strcmp(fn,'x'),
    xi = i;
  elseif strcmp(fn,'y'),
    yi = i;
  end
end

if isempty(xi) && isempty(yi),
  msg = 'None of the input data files corresponds to the x- or y-positions';
  return;
end

%% read the trx file
trx = ReadMWTBlobFile(blobsfile);
removeid = false(1,numel(trx));
trxids = [trx.id];

%% read data file(s)

for i = [xi,yi],
  [~,name] = myfileparts(datfiles{i});
  match = regexp(name,['\.(.+)\.',datext,'$'],'tokens','once');
  fn = match{1};
  perframefn = [fn,'_mm'];
  
  % read in data
  datacurr = ReadChoreDatFile(datfiles{i});
  
  % match up with identities
  dataids = [datacurr.id];
  newremoveid = ~ismember(trxids,dataids);
  removeid = removeid | newremoveid;
  [oldid,idx] = ismember(dataids,trxids);
  if any(~oldid),
    idx(~oldid) = [];
  end
  for jj = 1:numel(idx),
    j = idx(jj);
    [~,f0] = min(abs(trx(j).timestamps-datacurr(jj).timestamps(1)));
    [~,f1] = min(abs(trx(j).timestamps-datacurr(jj).timestamps(end)));
    if f1-f0+1 ~= numel(datacurr(jj).value),
      removeid(j) = true;
      continue;
%       msg = sprintf('timestamps don''t match up for datfile %s',datfiles{i});
%       return;
    end
    trx(j).(perframefn) = nan(1,trx(j).nframes);
    trx(j).(perframefn)(f0:f1) = datacurr(jj).value;
  end
  
end

if any(removeid),
  trx(removeid) = [];
end

%% 

if ~isempty(xi),
  idx = [trx.x_mm] > 0;
  tmp = [trx.x]./[trx.x_mm];
  s = nanstd(tmp);
  if s > 1,
    msg = sprintf('Standard deviation of pxpermm readings = %.1f too high, something is wrong',s);
    return;
  end
  pxpermm_x = nanmean(tmp(idx));
end
if ~isempty(yi),
  idx = [trx.y_mm] > 0;
  tmp = [trx.y]./[trx.y_mm];
  s = nanstd(tmp);
  if s > 1,
    msg = sprintf('Standard deviation of pxpermm readings = %.1f too high, something is wrong',s);
    return;
  end
  pxpermm_y = nanmean(tmp(idx));
end

if ~isempty(xi) && ~isempty(yi),
  pxpermm = mean([pxpermm_x,pxpermm_y]);
  msg = sprintf('Computed pxpermm from x and y positions in %s, %s, and %s',blobsfile{1},datfiles{xi},datfiles{yi});
elseif ~isempty(xi),
  pxpermm = pxpermm_x;
  msg = sprintf('Computed pxpermm from x positions in %s and %s',blobsfile{1},datfiles{xi});
else
  pxpermm = pxpermm_y;
  msg = sprintf('Computed pxpermm from y positions in %s and %s',blobsfile{1},datfiles{yi});
end

success = true;
