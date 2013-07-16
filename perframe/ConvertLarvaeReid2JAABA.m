function [success,msg] = ConvertLarvaeReid2JAABA(varargin)

success = false;
msg = {};

MINDT = .001;

trxfns = {'x','y','a','b','theta','x_mm','y_mm','a_mm','b_mm','theta_mm','dt','timestamps'};
trxunits = {
  parseunits('px'),...
  parseunits('px'),...
  parseunits('px'),...
  parseunits('px'),...
  parseunits('rad'),...
  parseunits('mm'),...
  parseunits('mm'),...
  parseunits('mm'),...
  parseunits('mm'),...
  parseunits('rad'),...
  parseunits('sec'),...
  parseunits('sec')};

perframeunits = {...
  'timestamp','s'
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
  'heady','px'
  'heading','deg'
  'bodyangle','deg'
  'headspeed','mm/s'
  'tailspeed','mm/s'
  'centerspeed','mm/s'
  'bodyanglespeed','deg/s'
  'concheadspeed','unit/s'
  'dfromsource','mm'
  'source','unit'
  'Chead','unit'
  'Ctail','unit'
  'sourceangle','deg'
  'curvingrate','deg/s'
  'bearingangforcurvingrate','deg'
  'cdotc','unit/s'};

[inmoviefile,...
  blobsfile,datfiles,kinmatfile,eventmatfile,...
  expdir,moviefilestr,trxfilestr,perframedirstr,...
  arenatype,arenacenterx,arenacentery,...
  arenaradius,arenawidth,arenaheight,...
  pxpermm,fps,overridefps,overridearena,...
  dosoftlink] = myparse(varargin,...
  'inmoviefile','',...
  'blobsfile','','datfiles',{},'kinmatfile','','eventmatfile','',...
  'expdir','','moviefilestr','movie.ufmf','trxfilestr','trx.mat','perframedirstr','perframe',...
  'arenatype','None','arenacenterx',0,'arenacentery',0,...
  'arenaradius',123,'arenawidth',123,'arenaheight',123,...
  'pxpermm',1,'fps',30,...
  'overridefps',false,'overridearena',false,...
  'dosoftlink',false);

% check that required inputs are given
if isempty(blobsfile),
  msg = 'Input blobs file is empty';
  return;
end
if ~iscell(blobsfile),
  blobsfile = {blobsfile};
end
for i = 1:numel(blobsfile),
  if ~exist(blobsfile{i},'file'),
    msg = sprintf('Input blobs file %s does not exist',blobsfile{i});
    return;
  end
end

for i = 1:numel(datfiles),
  if ~exist(datfiles{i},'file'),
    msg = sprintf('Input dat file %s does not exist',datfiles{i});
    return;
  end
end

if ~exist(kinmatfile,'file'),
  msg = sprintf('Input mat file %s does not exist',kinmatfile);
  return;
end

if ~isempty(eventmatfile) && ~exist(eventmatfile,'file'),
  msg = sprintf('Input mat file %s does not exist',eventmatfile);
  return;
end

% output file locations
moviefile = fullfile(expdir,moviefilestr);
trxfile = fullfile(expdir,trxfilestr);
perframedir = fullfile(expdir,perframedirstr);

%% read in the trx from the blobs file

trx = ReadMWTBlobFile(blobsfile);
nflies = numel(trx);

% make sure all the timestamps agree with each other
timestamps = unique([trx.timestamps]);
dts = diff(timestamps);
% sanity check to look for rounding errors
if any(dts < MINDT),
  msg = sprintf('Rounding error check failed: there are timestamps < %f apart',MINDT);
  return;
end

% add dt, fps
if ~overridefps,
  fps = 1/nanmedian(dts);
end
for i = 1:nflies,
  trx(i).dt = diff(trx(i).timestamps);
  trx(i).fps = fps;
end

trxids = [trx.id];

%% read in kin mat file

try
  load(kinmatfile,'kinData');
catch %#ok<CTCH>
  msg = sprintf('Error loading kinData from kin mat file %s',kinmatfile);
  return;
end

dataall = [];
for i = 1:numel(kinData), %#ok<USENS>
  dataall = structappend(dataall,kinData{i});
end
fns = setdiff(fieldnames(dataall),{'frame'});
for i = 1:numel(dataall),
  id = unique(dataall(i).objID);
  if numel(id) ~= 1,
    msg = sprintf('More than one objID found in kinData{%d}',i);
    return;
  end
  dataall(i).id = id; %#ok<AGROW>
  dataall(i).firstframe = min(dataall(i).frame)+1; %#ok<AGROW>
  dataall(i).endframe = max(dataall(i).frame)+1; %#ok<AGROW>
  dataall(i).nframes = dataall(i).endframe-dataall(i).firstframe+1; %#ok<AGROW>
  [ism,fidx] = ismember(dataall(i).firstframe-1:dataall(i).endframe-1,dataall(i).frame);
  if any(~ism),
    fprintf('%d frames with no data for target %d\n',nnz(~ism),id);
  end
  
  n = dataall(i).nframes;
  for j = 1:numel(fns),
    fn = fns{j};
    tmp = nan(1,n);
    tmp(ism) = dataall(i).(fn)(fidx);
    dataall(i).(fn) = tmp; %#ok<AGROW>
  end
  
end

%% read in the eventData

try
  load(eventmatfile,'eventData');
catch %#ok<CTCH>
  msg = sprintf('Error loading eventData from kin mat file %s',eventmatfile);
  return;
end

if numel(eventData) > numel(kinData), %#ok<USENS>
  msg = sprintf('Less trajectories for kinData file %s than eventData file %s',kinmatfile,eventmatfile);
  return;
end

% add eventData to dataall
fns = fieldnames(eventData{1});
for i = 1:numel(eventData),
  datacurr = eventData{i};
  n = dataall(i).nframes;
  for j = 1:numel(fns),
    fn = fns{j};
    if max(datacurr.(fn)) > n,
      msg = sprintf('Maximum index for %s bigger than number of frames in trajectory %d',fn,i);
      return;
    end
    tmp = zeros(1,n);
    tmp(datacurr.(fn)) = 1;
    dataall(i).(['event_',fn]) = tmp; %#ok<AGROW>
  end
end



%% match up ids

% find starts of videos

% which sequence does the kindata correspond to
dataids = [dataall.id];
isstart = [true,dataids(1:end-1)>dataids(2:end)];
idxstart = find(isstart);
idxend = [idxstart(2:end)-1,numel(dataids)];

maxmatches = -1;
for i = 1:numel(idxstart),
  idxcurr = idxstart(i):idxend(i);
  dataids_curr = [dataall(idxcurr).id];
  [ism,order] = ismember(dataids_curr,trxids);
  ismatch = [dataall(idxcurr(ism)).nframes] == [trx(order(ism)).nframes];
  fracmatches = nnz(ismatch) / nnz(ism);
  if fracmatches > maxmatches,
    bestmatch = i;
    maxmatches = fracmatches;
  end
end

msg{end+1} = sprintf('Choosing to match kindata sequence %d with %s, fraction of ids for which nframes matches = %f',bestmatch,blobsfile{1},maxmatches);
dataall = dataall(idxstart(bestmatch):idxend(bestmatch));
dataids = [dataall.id];

if numel(unique(dataids)) < numel(dataids),
  msg = 'Some repeat ids found in dataids';
  return;
end
badtrxids = ~ismember(trxids,dataids);
baddataids = ~ismember(dataids,trxids);
if any(badtrxids) || any(baddataids),
  msg{end+1} = sprintf('%d ids in blobs not in kindata, %d ids in kindata not in blobs',nnz(badtrxids),nnz(baddataids));
end
dataall(baddataids) = [];
dataids(baddataids) = [];
trx(badtrxids) = [];
trxids(badtrxids) = [];
[~,idx] = ismember(dataids,trxids);
dataall = dataall(idx);

badids = false(1,numel(trx));
for i = 1:numel(trx),
  if trx(i).firstframe ~= dataall(i).firstframe || ...
      trx(i).endframe ~= dataall(i).endframe,
    badids(i) = true;
  end
end
if any(badids),
  msg{end+1} = sprintf('%d id nframes do not match between kindata and blobs files for some ids, removing these',nnz(badids));
  dataall(badids) = [];
  trx(badids) = [];
end

if isempty(trx),
  msg{end+1} = 'After removing mismatches between kindata and blobs, no trajectories left.';
  return;
end

%% add trajectory features that are in kindata files

for i = 1:numel(trx),
  if isfield(dataall,'headx'),
    trx(i).headx = dataall(i).headx;
  end
  if isfield(dataall,'heady'),
    trx(i).heady = dataall(i).heady;
  end
  if isfield(dataall,'tailx'),
    trx(i).tailx = dataall(i).tailx;
  end
  if isfield(dataall,'taily'),
    trx(i).tailx = dataall(i).taily;
  end
  if isfield(dataall,'cx'),
    trx(i).x_mm = dataall(i).cx;
  end
  if isfield(dataall,'cy'),
    trx(i).y_mm = dataall(i).cy;
  end
end

%% convert to mm

if strcmpi(arenatype,'None'),
  arenacenterx = 0;
  arenacentery = 0;
end

if isfield(trx,'x_mm'),
  if ~overridearena,
    idx = [trx.x_mm] > 0;
    tmp = [trx.x]./[trx.x_mm];
    s = nanstd(tmp);
    if s < 1,
      tmp = nanmean(tmp(idx));
      if ~isnan(tmp),
        pxpermm = tmp;
      end
    end
  end
  dx = [trx.x_mm] - ([trx.x]-arenacenterx)/pxpermm;
  if nanstd(dx) > 1,
    msg = sprintf('pxpermm does not match the relationship between x and x_mm');
    return;
  end
  arenacenterx_mm = nanmean(dx);
else
  arenacenterx_mm = 0;
end

if isfield(trx,'y_mm'),
  dy = [trx.y_mm] - ([trx.y]-arenacentery)/pxpermm;
  if nanstd(dy) > 1,
    msg = sprintf('pxpermm does not match the relationship between y and y_mm');
    return;
  end
  arenacentery_mm = nanmean(dy);
else
  arenacentery_mm = 0;
end


for i = 1:numel(trx),
  trx(i).pxpermm = pxpermm;
  if ~isfield(trx,'x_mm'),
    trx(i).x_mm = (trx(i).x - arenacenterx) / pxpermm;
  end
  if ~isfield(trx,'y_mm'),
    trx(i).y_mm = (trx(i).y - arenacentery) / pxpermm;
  end
  trx(i).a_mm = trx(i).a / pxpermm;
  trx(i).b_mm = trx(i).b / pxpermm;
  trx(i).theta_mm = trx(i).theta;
  if isfield(trx,'xspine_mm'),
    trx(i).xspine = trx(i).xspine_mm * pxpermm;
  end
  if isfield(trx,'yspine_mm'),
    trx(i).yspine = trx(i).yspine_mm * pxpermm;
  end
end

%% over-ride timestamps

if overridefps,
  timestamps = timestamps(1)+(0:numel(timestamps)-1)/fps;
  for i = 1:numel(trx),
    trx(i).timestamps = timestamps(trx(i).firstframe:trx(i).endframe);
    trx(i).dt(:) = 1/fps;
  end
end

%% set landmark parameters

arenaradius_mm = arenaradius / pxpermm;
arenawidth_mm = arenawidth / pxpermm;
arenaheight_mm = arenaheight / pxpermm;
trx = SetLandmarkParameters(trx,arenatype,arenacenterx_mm,arenacentery_mm,...
  arenaradius_mm,arenawidth_mm,arenaheight_mm); 

%% create directory

% create the experiment directory
if ~exist(expdir,'dir'),
  [success1,msg] = mkdir(expdir);
  if ~success1,
    return;
  end
end

%% create per-frame directory

if ~exist(perframedir,'dir'),
  [success1,msg] = mkdir(perframedir);
  if ~success1,
    return;
  end
end

%% save the per-frame data

perframefns = setdiff(fieldnames(dataall),{'frame','objID','objN','goodN','persists','id','firstframe','endframe','nframes','times'});

for i = 1:numel(perframefns),
  fn = perframefns{i};
  data = {dataall.(fn)}; %#ok<NASGU>
  j = find(strcmp(perframeunits(:,1),fn),1);
  if ~isempty(j),
    units = parseunits(perframeunits{j,2}); %#ok<NASGU>
  else
    units = parseunits('unit'); %#ok<NASGU>
  end
  outfile = fullfile(perframedir,[fn,'.mat']);
  try
    save(outfile,'data','units');
  catch ME,
    msg = getReport(ME);
    return;
  end
end

for i = 1:numel(trxfns),
  fn = trxfns{i};
  if isfield(trx,fn),
    data = {trx.(fn)}; %#ok<NASGU>
    units = trxunits{i}; %#ok<NASGU>
    %units = parseunits('unit'); %#ok<NASGU>
    outfile = fullfile(perframedir,[fn,'.mat']);
    try
      save(outfile,'data','units');
    catch ME,
      msg = getReport(ME);
      return;
    end
  end
end

%% create the trx file

try
  save(trxfile,'trx','timestamps');
catch ME,
  msg = getReport(ME);
  return;
end

if ~exist(trxfile,'file'),
  msg = sprintf('Failed to save trx to file %s',trxfile);
  return;
end

%% copy/soft-link movie

if ~isempty(inmoviefile),
  
  if strcmp(fullfile(inmoviefile),fullfile(moviefile)),
    fprintf('Input and out movie files are the same, not copying/linking.\n');
  else
    
    if dosoftlink,
      if isunix,
        cmd = sprintf('ln -s %s %s',inmoviefile,moviefile);
        unix(cmd);
        % test to make sure it worked
        [status,result] = unix(sprintf('readlink %s',moviefile));
        result = strtrim(result);
        if status ~= 0 || ~strcmp(result,inmoviefile),
          warndlg(sprintf('Failed to make soft link, copying %s to %s instead',inmoviefile,moviefile));
          dosoftlink = false;
        end
      elseif ispc,
        cmd = sprintf('mkshortcut.vbs /target:"%s" /shortcut:"%s"',inmoviefile,moviefile);
        fprintf('Making a Windows shortcut file at "%s" with target "%s"\n',inmoviefile,moviefile);
        system(cmd);
        % test to make sure that worked
        [equalmoviefile,didfind] = GetPCShortcutFileActualPath(moviefile);
        if ~didfind || ~strcmp(equalmoviefile,inmoviefile),
          warndlg(sprintf('Failed to make shortcut, copying %s to %s instead',inmoviefile,moviefile));
          dosoftlink = false;
        end
      else
        warndlg(sprintf('Unknown OS, not soft-linking movie file %s',inmoviefile));
        dosoftlink = false;
      end
    end
    
    if ~dosoftlink,
      [success1,msg1] = copyfile(inmoviefile,moviefile);
      if ~success1,
        msg = msg1;
        success = false;
        return;
      end
    end
    
  end
end

success = true;