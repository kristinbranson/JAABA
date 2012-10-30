function [success,msg] = ConvertCtrax2JAABA(varargin)

success = false;
msg = '';

[inmoviefile,intrxfile,annfile,...
  expdir,moviefilestr,trxfilestr,perframedirstr,...
  arenatype,arenacenterx,arenacentery,...
  arenaradius,arenawidth,arenaheight,...
  pxpermm,fps,overridefps,...
  dosoftlink] = myparse(varargin,...
  'inmoviefile','','intrxfile','','annfile','',...
  'expdir','','moviefilestr','movie.ufmf','trxfilestr','trx.mat','perframedirstr','perframe',...
  'arenatype','None','arenacenterx',0,'arenacentery',0,...
  'arenaradius',123,'arenawidth',123,'arenaheight',123,...
  'pxpermm',1,'fps',30,...
  'overridefps',false,...
  'dosoftlink',false);

% check that required inputs are given
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

% load in the trx
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

% compute frame rate if not overriding and possible
if ~overridefps && ~isempty(timestamps),
  fps = 1/median(diff(timestamps));
end

if strcmpi(arenatype,'None'),
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

if isempty(timestamps),
  timestamps = (0:max([trx.endframe]))*1/fps; %#ok<NASGU>
end

% set landmark parameters
arenaradius_mm = arenaradius / pxpermm;
arenawidth_mm = arenawidth / pxpermm;
arenaheight_mm = arenaheight / pxpermm;
arenacenterx_mm = 0;
arenacentery_mm = 0;
trx = SetLandmarkParameters(trx,arenatype,arenacenterx_mm,arenacentery_mm,...
  arenaradius_mm,arenawidth_mm,arenaheight_mm); %#ok<NASGU>

% create the experiment directory
if ~exist(expdir,'dir'),
  [success1,msg1] = mkdir(expdir);
  if ~success1,
    msg = msg1;
    return;
  end
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

% copy/soft-link movie
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

% make per-frame directory
if ~exist(perframedir,'dir'),
  [success1,msg1] = mkdir(perframedir);
  if ~success1,
    msg = msg1;
    return;
  end
end

success = true;