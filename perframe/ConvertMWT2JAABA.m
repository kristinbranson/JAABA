function [success,msg] = ConvertMWT2JAABA(varargin)

success = false;
msg = {};

MINDT = .001;
MAXTIMESTAMPERR = .01;
datext = 'dat';
perframedict = {'x','x_mm'
  'y','y_mm'
  'area','area_mm'};
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

[inmoviefile,...
  blobsfile,spinefile,datfiles,...
  expdir,moviefilestr,trxfilestr,perframedirstr,...
  arenatype,arenacenterx,arenacentery,...
  arenaradius,arenawidth,arenaheight,...
  pxpermm,fps,overridefps,overridearena,...
  dosoftlink] = myparse(varargin,...
  'inmoviefile','',...
  'blobsfile','','spinefile','','datfiles',{},...
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

if ~isempty(spinefile) && ~exist(spinefile,'file'),
  msg = sprintf('Input spine file %s does not exist',spinefile);
  return;
end

for i = 1:numel(datfiles),
  if ~exist(datfiles{i},'file'),
    msg = sprintf('Input dat file %s does not exist',datfiles{i});
    return;
  end
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

%% read in spine file

haveremoved = false;
removeid = false(1,numel(trx));
trxids = [trx.id];

if ~isempty(spinefile),
  
  % read in data
  datacurr = ReadChoreSpineFile(spinefile);
  
  % match up with identities
  dataids = [datacurr.id];
  newremoveid = ~ismember(trxids,dataids);
  if haveremoved && any(newremoveid ~= removeid),
    warning('removeid and newremoveid do not match');
  end
  removeid = removeid | newremoveid;
  haveremoved = true;
  [oldid,idx] = ismember(dataids,trxids);
  if any(~oldid),
    warning('ids found in spine file %s that are not in trx',spinefile);
    idx(~oldid) = [];
  end
  for jj = 1:numel(idx),
    j = idx(jj);
    [~,f0] = min(abs(trx(j).timestamps-datacurr(jj).timestamps(1)));
    [~,f1] = min(abs(trx(j).timestamps-datacurr(jj).timestamps(end)));
    if f1-f0+1 ~= size(datacurr(jj).xspine,2),
      msg = sprintf('timestamps don''t match up for spine file %s',spinefile);
      return;
    end
    maxerr = max(abs(trx(j).timestamps(f0:f1)-datacurr(jj).timestamps));
    if maxerr > MAXTIMESTAMPERR,
      msg = sprintf('timestamps don''t match up for spine file %s',spinefile);
      return
    end
    trx(j).xspine_mm = nan(size(datacurr(jj).xspine,1),trx(j).nframes);
    trx(j).yspine_mm = nan(size(datacurr(jj).xspine,1),trx(j).nframes);
    trx(j).xspine_mm(:,f0:f1) = datacurr(jj).xspine;
    trx(j).yspine_mm(:,f0:f1) = datacurr(jj).yspine;
  end

end

%% read dat files
perframefns = cell(numel(datfiles),2);

for i = 1:numel(datfiles),
  [~,name] = myfileparts(datfiles{i});
  match = regexp(name,['\.(.+)\.',datext,'$'],'tokens','once');
  if isempty(match),
    msg = sprintf('Could not parse dat file type for file %s',name);
    return;
  end
  fn = match{1};
  j = find(strcmp(fn,perframedict(:,1)),1);
  if ~isempty(j),
    fn = perframedict{j,2};
  end
  if ismember(fn,trxfns),
    perframefn = fn;
  else
    perframefn = ['perframe_',fn];
  end
  perframefns(i,:) = {fn,perframefn};
  
  % read in data
  datacurr = ReadChoreDatFile(datfiles{i});
  
  % match up with identities
  dataids = [datacurr.id];
  newremoveid = ~ismember(trxids,dataids);
  if haveremoved && any(newremoveid ~= removeid),
    warning('removeid and newremoveid do not match');
  end
  removeid = removeid | newremoveid;
  haveremoved = true;
  [oldid,idx] = ismember(dataids,trxids);
  if any(~oldid),
    warning('ids found in dat file %s that are not in trx',datfiles{i});
    idx(~oldid) = [];
  end
  for jj = 1:numel(idx),
    j = idx(jj);
    [~,f0] = min(abs(trx(j).timestamps-datacurr(jj).timestamps(1)));
    [~,f1] = min(abs(trx(j).timestamps-datacurr(jj).timestamps(end)));
    if f1-f0+1 ~= numel(datacurr(jj).value),
      msg = sprintf('timestamps don''t match up for datfile %s',datfiles{i});
      return;
    end
    maxerr = max(abs(trx(j).timestamps(f0:f1)-datacurr(jj).timestamps));
    if maxerr > MAXTIMESTAMPERR,
      msg = sprintf('timestamps don''t match up for datfile %s',datfiles{i});
      return;
    end
    trx(j).(perframefn) = nan(1,trx(j).nframes);
    trx(j).(perframefn)(f0:f1) = datacurr(jj).value;
  end
  
end

if any(removeid),
  trx(removeid) = [];
end


% replace first frames with nans
for i = 1:size(perframefns,1),
  fn = perframefns{i,2};
  off = inf;
  for j = 1:numel(trx),
    if isempty(trx(j).(fn)),
      continue;
    end
    if trx(j).(fn)(1) ~= 0,
      off = 0;
      break;
    end
    if all(trx(j).(fn) == 0),
      continue;
    end
    k = find(trx(j).(fn) ~= 0,1);
    off = min(off,k-1);
  end
  if off > 0,
    fprintf('First %d frames of data for %s are always 0, setting to nan\n',off,fn);
  end
  for j = 1:numel(trx),
    trx(j).(fn)(1:min(numel(trx(j).(fn)),off)) = nan;
  end
end

% replace last frames with nans
for i = 1:size(perframefns,1),
  fn = perframefns{i,2};
  off = inf;
  for j = 1:numel(trx),
    if isempty(trx(j).(fn)),
      continue;
    end
    if trx(j).(fn)(end) ~= 0,
      off = 0;
      break;
    end
    if all(trx(j).(fn) == 0),
      continue;
    end
    k = find(trx(j).(fn) ~= 0,1,'last');
    off = min(off,numel(trx(j).(fn))-k);
  end
  if off > 0,
    fprintf('Last %d frames of data for %s are always 0, setting to nan\n',off,fn);
  end
  for j = 1:numel(trx),
    trx(j).(fn)(max(1,numel(trx(j).(fn))-off+1):end) = nan;
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

for i = 1:size(perframefns,1),
  fn = perframefns{i,1};
  perframefn = perframefns{i,2};
  if ~ismember(fn,trxfns),
    data = {trx.(perframefn)}; %#ok<NASGU>
    % units not set yet
    units = parseunits('unit'); %#ok<NASGU>
    outfile = fullfile(perframedir,[fn,'.mat']);
    try
      save(outfile,'data','units');
    catch ME,
      msg = getReport(ME);
      return;
    end
    trx = rmfield(trx,perframefn);
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