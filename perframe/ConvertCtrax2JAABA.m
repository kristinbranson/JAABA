function [success,msg] = ConvertCtrax2JAABA(varargin)

success = false;
msg = {};

[inmoviefile,intrxfile,annfile,...
  expdir,moviefilestr,trxfilestr,perframedirstr,...
  arenatype,arenacenterx,arenacentery,...
  arenaradius,arenawidth,arenaheight,...
  pxpermm,fps,overridefps,overridearena,...
  dosoftlink,...
  dofliplr,doflipud,dotransposeimage,...
  inperframedir] = myparse(varargin,...
  'inmoviefile','','intrxfile','','annfile','',...
  'expdir','','moviefilestr','movie.ufmf','trxfilestr','trx.mat','perframedirstr','perframe',...
  'arenatype','None','arenacenterx',0,'arenacentery',0,...
  'arenaradius',123,'arenawidth',123,'arenaheight',123,...
  'pxpermm',1,'fps',30,...
  'overridefps',false,'overridearena',false,...
  'dosoftlink',false,...
  'fliplr',false,...
  'flipud',false,...
  'dotransposeimage',false,...
  'inperframedir','');

wingunits = struct(...
  'nwingsdetected',parseunits('unit'),...
  'wing_areal',parseunits('px^2'),...
  'wing_arear',parseunits('px^2'),...
  'wing_trough_angle',parseunits('rad'));

%% check that required inputs are given
if isempty(inmoviefile),
  msg = 'Input movie file is empty';
  return;
end
if isempty(intrxfile),
  msg = 'Input trx mat file is empty';
  return;
end
% if isempty(annfile),
%   msg = 'Input ann file is empty';
%   return;
% end

% output file locations
moviefile = fullfile(expdir,moviefilestr);
trxfile = fullfile(expdir,trxfilestr);
perframedir = fullfile(expdir,perframedirstr);

%% load in the trx
try
  [trx,~,success1,timestamps] = load_tracks(intrxfile,moviefile,...
    'dosave',false,'annname',annfile,'verbose',false);
  if ~success1,
    msg = sprintf('Could not load trx from file %s',intrxfile);
    return;
  end
catch ME,
  msg = sprintf('Could not load trx from file %s: %s',intrxfile,getReport(ME));
  return;
end

iswings = isfield(trx,'xwingl');

if iswings,
  if isempty(inperframedir),
    msg = 'Per-frame directory containing wing data is not set';
    return;
  end
  fns = fieldnames(wingunits);
  for i = 1:numel(fns),
    fn = fns{i};
    filename = fullfile(inperframedir,[fn,'.mat']);
    if ~exist(filename,'file'),
      msg = sprintf('Wing feature file %s does not exist.',filename);
      return;
    end
  end
end

if dofliplr || doflipud,
  try
    % get image size
    [readframe,~,fid] = get_readframe_fcn(inmoviefile);
    im = readframe(1);
    [imheight,imwidth,~] = size(im);
    if fid > 1,
      fclose(fid);
    end
  catch ME,
    msg = getReport(ME);
    return;
  end
end

% flip if necessary
if dofliplr,
  for i = 1:numel(trx),
    trx(i).x = imwidth-trx(i).x+1;
    trx(i).theta = modrange(pi - trx(i).theta,-pi,pi);
    if iswings,
      trx(i).xwingl = imwidth-trx(i).xwingl+1;
      trx(i).xwingr = imwidth-trx(i).xwingr+1;
    end
  end
  msg{end+1} = 'Flipped the trajectories left-right.';
end
if doflipud,
  for i = 1:numel(trx),
    trx(i).y = imheight-trx(i).y+1;
    trx(i).theta = modrange(-trx(i).theta,-pi,pi);
    if iswings,
      trx(i).ywingl = imheight-trx(i).xwingl+1;
      trx(i).ywingr = imheight-trx(i).xwingr+1;
    end
  end
  msg{end+1} = 'Flipped the trajectories up-down.';
end
if dotransposeimage,
  for i = 1:numel(tx),
    tmp = trx(i).x;
    trx(i).x = trx(i).y;
    trx(i).y = tmp;
    c = cos(trx(i).theta);
    s = sin(trx(i).theta);
    trx(i).theta = atan2(c,s);
    if iswings,
      tmp = trx(i).xwingl;
      trx(i).xwingl = trx(i).ywingl;
      trx(i).ywingl = tmp;
      tmp = trx(i).xwingr;
      trx(i).xwingr = trx(i).ywingr;
      trx(i).ywingr = tmp;
    end
  end
  msg{end+1} = 'Switched x and y in trajectories.';
end

% compute frame rate if not overriding and possible
if ~overridefps && ~isempty(timestamps),
  fps = 1/nanmedian(diff(timestamps));
  msg{end+1} = sprintf('Computed frame rate = %f fps from timestamps.',fps);
else
  msg{end+1} = sprintf('Using input frame rate %f fps to set timestamps',fps);
end

if ~overridearena && isfield(trx,'pxpermm'),
  tmp = nanmean([trx.pxpermm]);
  if ~isnan(tmp),
    pxpermm = tmp;
    msg{end+1} = sprintf('Read pxpermm = %f from trxfile.',pxpermm);
  end
else
  msg{end+1} = sprintf('Using input pixels-per-mm scaling = %f.',pxpermm);
end

if ~overridearena && isfield(trx,'arena') && all(isfield(trx(1).arena,{'x','y','r'})) && ...
  ~isempty(trx(1).arena.r) && ~isempty(trx(1).arena.x) && ~isempty(trx(1).arena.y),
  arenaradius = trx(1).arena.r;
  arenacenterx = trx(1).arena.x;
  arenacentery = trx(1).arena.y;
  arenatype = 'Circle';
  msg{end+1} = sprintf('Read circle arena parameters from trxfile. Center = (%.1f, %.1f) px, radius = %.1f px.',arenacenterx,arenacentery,arenaradius);
elseif strcmpi(arenatype,'None'),
  arenacenterx = 0;
  arenacentery = 0;
end

% apply to trajectories
for fly = 1:length(trx),
  
  % scale and translate position
  trx(fly).x_mm = (trx(fly).x - arenacenterx) / pxpermm;
  trx(fly).y_mm = (trx(fly).y - arenacentery) / pxpermm;
  trx(fly).a_mm = trx(fly).a / pxpermm;
  trx(fly).b_mm = trx(fly).b / pxpermm;
  trx(fly).theta_mm = trx(fly).theta;
  
  % add in dt
  if overridefps || isempty(timestamps),
    trx(fly).dt = repmat(1/fps,[1,trx(fly).nframes-1]);
  else
    trx(fly).dt = diff(timestamps(trx(fly).firstframe:trx(fly).endframe));
    trx(fly).dt(isnan(trx(fly).dt)) = 1/fps;
  end

  trx(fly).fps = 1/median(trx(fly).dt);
  trx(fly).pxpermm = pxpermm;
  
end

msg{end+1} = sprintf('Read trx for %d flies, frame range [%d,%d]',numel(trx),min([trx.firstframe]),max([trx.endframe]));

if isempty(timestamps) || overridefps,
  timestamps = (0:max([trx.endframe]))*1/fps; %#ok<NASGU>
end

% set landmark parameters
arenaradius_mm = arenaradius / pxpermm;
arenawidth_mm = arenawidth / pxpermm;
arenaheight_mm = arenaheight / pxpermm;
arenacenterx_mm = 0;
arenacentery_mm = 0;
trx = SetLandmarkParameters(trx,arenatype,arenacenterx_mm,arenacentery_mm,...
  arenaradius_mm,arenawidth_mm,arenaheight_mm); 

if strcmpi(arenatype,'Circle')
  msg{end+1} = sprintf('Set circular arena parameters. Center = (0,0) mm, radius = %.1f mm.',arenaradius_mm);
elseif strcmpi(arenatype,'Rectangle'),
  msg{end+1} = sprintf('Set rectangular arena parameters. Center = (0,0) mm, width= %.1f mm, height = %.1f mm',arenawidth_mm,arenaheight_mm);  
end
msg{end+1} = sprintf('x_mm ranges within [%.1f,%.1f], y_mm ranges within [%.1f,%.1f]',...
  min([trx.x_mm]),max([trx.x_mm]),min([trx.y_mm]),max([trx.y_mm]));

% create the experiment directory
if ~exist(expdir,'dir'),
  [success1,msg1] = mkdir(expdir);
  if ~success1,
    msg = msg1;
    return;
  end
  msg{end+1} = sprintf('Created experiment directory %s.',expdir);
end

% save the trx file
try
  save(trxfile,'trx','timestamps');
catch ME,
  msg = sprintf('Could not save to file %s: %s',trxfile,getReport(ME));
  return;
end
if ~exist(trxfile,'file'),
  msg = sprintf('Failed to save trx to file %s',trxfile);
  return;
end
msg{end+1} = sprintf('Saved trx to file %s.',trxfile);

% copy/soft-link movie
if dosoftlink,
  if exist(moviefile,'file'),
    delete(moviefile);
  end
  if isunix,
    cmd = sprintf('ln -s %s %s',inmoviefile,moviefile);
    unix(cmd);
    % test to make sure it worked
    [status,result] = unix(sprintf('readlink %s',moviefile));
    result = strtrim(result);
    if status ~= 0 || ~strcmp(result,inmoviefile),
      res = questdlg(sprintf('Failed to make soft link. Copy %s to %s instead?',inmoviefile,moviefile));
      if ~strcmpi(res,'Yes'),
        msg = sprintf('Failed to make soft link from %s to %s.',inmoviefile,moviefile);
        return;
      end
      dosoftlink = false;
    end
  elseif ispc,
    if exist([moviefile,'.lnk'],'file'),
      delete([moviefile,'.lnk']);
    end
    cmd = sprintf('mkshortcut.vbs /target:"%s" /shortcut:"%s"',inmoviefile,moviefile);
    fprintf('Making a Windows shortcut file at "%s" with target "%s"\n',inmoviefile,moviefile);
    system(cmd);
    % test to make sure that worked
    [equalmoviefile,didfind] = GetPCShortcutFileActualPath(moviefile);
    if ~didfind || ~strcmp(equalmoviefile,inmoviefile),
      res = questdlg(sprintf('Failed to make shortcut. Copy %s to %s instead?',inmoviefile,moviefile));
      if ~strcmpi(res,'Yes'),
        msg = sprintf('Failed to make shortcut from %s to %s.',inmoviefile,moviefile);
        return;
      end
      dosoftlink = false;
    end
  else
    res = questdlg(sprintf('Unknown OS, cannot soft-link movie file %s. Copy instead?',inmoviefile));
    if ~strcmpi(res,'Yes'),
      msg = sprintf('Failed to make softlink from %s to %s.',inmoviefile,moviefile);
      return;
    end
    dosoftlink = false;
  end  
  if dosoftlink,
    msg{end+1} = sprintf('Made a link to movie file %s at %s',inmoviefile,moviefile);
  end
end
  
if ~dosoftlink,
  if ispc,
    if exist([moviefile,'.lnk'],'file'),
      delete([moviefile,'.lnk']);
    end
  end
  if exist(moviefile,'file'),
    delete(moviefile);
  end
  [success1,msg1] = copyfile(inmoviefile,moviefile);
  if ~success1,
    msg = msg1;
    success = false;
    return;
  end
  msg{end+1} = sprintf('Copied movie file %s to %s',inmoviefile,moviefile);
end

% make per-frame directory
if ~exist(perframedir,'dir'),
  [success1,msg1] = mkdir(perframedir);
  if ~success1,
    msg = msg1;
    return;
  end
end

% copy over wing features
if iswings,
  fns = fieldnames(wingunits);
  for i = 1:numel(fns),
    fn = fns{i};
    filenamein = fullfile(inperframedir,[fn,'.mat']);
    filenameout = fullfile(perframedir,[fn,'.mat']);
    if exist(filenameout,'file'),
      delete(filenameout);
    end
    [success1,msg1] = copyfile(filenamein,filenameout);
    if ~success1,
      msg = msg1;
      return;
    end
  end
end

success = true;