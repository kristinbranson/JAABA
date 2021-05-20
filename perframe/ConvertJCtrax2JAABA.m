function [success,msg] = ConvertJCtrax2JAABA(varargin)

success = false;
msg = {};

[inmoviefile,intrxfile,infofile,inctraxfile,...
  expdir,moviefilestr,trxfilestr,perframedirstr,...
  dosoftlink,...
  doflipud,dofliplr,dotransposeimage,...
  sex,movieheight,moviewidth] = myparse(varargin,...
  'inmoviefile','','intrxfile','',...
  'infofile','',...
  'inctraxfile','',...
  'expdir','','moviefilestr','movie.avi','trxfilestr','trx.mat','perframedirstr','perframe',...
  'dosoftlink',false,...
  'doflipud',false,'dofliplr',false,'dotransposeimage',false,...
  'sex',{},...
  'movieheight',[],'moviewidth',[]);

pfunits = struct(...
  'nwingsdetected',parseunits('unit'),...
  'wing_areal',parseunits('mm^2'),...
  'wing_arear',parseunits('mm^2'),...
  'wing_trough_angle',parseunits('rad'),...
  'fgarea',parseunits('mm^2'),...
  'imgcontrast',parseunits('unit'),...
  'minfgdist',parseunits('unit'),...
  'area_mm',parseunits('mm^2'));

%% check that required inputs are given
if isempty(inmoviefile),
  msg = 'Input movie file is empty';
  return;
end
if isempty(intrxfile),
  msg = 'Input trx mat file is empty';
  return;
end
if isempty(infofile) && isempty(inctraxfile),
  msg = 'Info mat file is empty';
  return;
end

% output file locations
moviefile = fullfile(expdir,moviefilestr);
trxfile = fullfile(expdir,trxfilestr);
perframedir = fullfile(expdir,perframedirstr);

%% load in data

td = load(intrxfile);
[nflies,nframes,nfeats] = size(td.trk.data);
if doflipud || dofliplr && (isempty(movieheight) || isempty(moviewidth)),

  % read frame size, moviename is either the name of the movie or an mmreader object
  [readframe,~,fid] = get_readframe_fcn(moviefile);
  im = readframe(1);
  if isempty(movieheight),
    movieheight = size(im,1);
  end
  if isempty(moviewidth),
    moviewidth = size(im,2);
  end
  if fid > 1,
    fclose(fid);
  end
end

msg{end+1} = sprintf('Loaded trajectories from file %s',trxfile);

%% read arena parameters, pxpermm, and fps

if ~isempty(infofile),

  infodata = load(infofile);
  
  fps = infodata.hinfo.FPS;
  pxpermm = infodata.hinfo.PPM;
  r = any(infodata.hinfo.mask,2);
  arenar1 = find(r,1);
  arenar2 = find(r,1,'last');
  c = any(infodata.hinfo.mask,1);
  arenac1 = find(c,1);
  arenac2 = find(c,1,'last');
  arenacenterx = (arenac1+arenac2)/2;
  arenacentery = (arenar1+arenar2)/2;
  arenaradius = ((arenar2-arenar1+1)+(arenac2-arenac1+1))/2;
  arenatype = 'Circle';

else
  
  ctraxdata = load(inctraxfile);
  fps = median([ctraxdata.trx.fps]);
  pxpermm = median([ctraxdata.trx.pxpermm]);
  arenacenterx = median(cellfun(@(x) x.x, {ctraxdata.trx.arena}));
  arenacentery = median(cellfun(@(x) x.y, {ctraxdata.trx.arena}));
  arenaradius = median(cellfun(@(x) x.r, {ctraxdata.trx.arena}));
  arenatype = 'Circle';
  
end

%% convert

idxnums = struct;
idxnums.x = find(strcmp(td.trk.names,'pos x'));
idxnums.y = find(strcmp(td.trk.names,'pos y'));
idxnums.theta = find(strcmp(td.trk.names,'ori'));
idxnums.a = find(strcmp(td.trk.names,'major axis len'));
idxnums.b = find(strcmp(td.trk.names,'minor axis len'));
idxnums.xwingl = find(strcmp(td.trk.names,'wing l x'));
idxnums.ywingl = find(strcmp(td.trk.names,'wing l y'));
idxnums.xwingr = find(strcmp(td.trk.names,'wing r x'));
idxnums.ywingr = find(strcmp(td.trk.names,'wing r y'));
idxnums.area = find(strcmp(td.trk.names,'body area'));
idxnums.wing_anglel = find(strcmp(td.trk.names,'wing l ang'));
idxnums.wing_angler = find(strcmp(td.trk.names,'wing r ang'));
idxnums.wing_lengthl = find(strcmp(td.trk.names,'wing l len'));
idxnums.wing_lengthr = find(strcmp(td.trk.names,'wing r len'));
idxnums.fgarea = find(strcmp(td.trk.names,'fg area'));
idxnums.imgcontrast = find(strcmp(td.trk.names,'img contrast'));
idxnums.minfgdist = find(strcmp(td.trk.names,'min fg dist'));

if isfield(td.trk,'frame_ids'),
  timestamps = td.trk.frame_ids/fps;
elseif ~isempty(inctraxfile),
  timestamps = ctraxdata.timestamps;
  td.trk.frame_ids = (0:size(td.trk.data,2)-1);
else
  timestamps = (0:size(td.trk.data,2)-1)/fps;
  td.trk.frame_ids = (0:size(td.trk.data,2)-1);
end
  

trx = struct;
pd = struct;
for fly = 1:nflies,

  trx(fly).id = fly;
  
  i0 = find(~isnan(td.trk.data(fly,:,1)),1);
  i1 = find(~isnan(td.trk.data(fly,:,1)),1,'last');
  trx(fly).firstframe = i0+td.trk.frame_ids(1);
  trx(fly).endframe = i1+td.trk.frame_ids(1);
  trx(fly).nframes = i1-i0+1;
  trx(fly).off = 1-trx(fly).firstframe;
  trx(fly).timestamps = timestamps(i0:i1);
  trx(fly).dt = diff(trx(fly).timestamps);

  if ~iscell(sex),
    trx(fly).sex = sex;
  elseif numel(sex) >= fly,
    trx(fly).sex = sex{fly};
  else
    trx(fly).sex = '?';
  end
  trx(fly).pxpermm = pxpermm;
  trx(fly).fps = fps;
  
  % pixel-based stuff
  trx(fly).x = td.trk.data(fly,i0:i1,idxnums.x);
  trx(fly).y = td.trk.data(fly,i0:i1,idxnums.y);
  trx(fly).theta = -td.trk.data(fly,i0:i1,idxnums.theta);
  trx(fly).a = td.trk.data(fly,i0:i1,idxnums.a)/4;
  trx(fly).b = td.trk.data(fly,i0:i1,idxnums.b)/4;
  
  trx(fly).xwingl = td.trk.data(fly,i0:i1,idxnums.xwingl);
  trx(fly).ywingl = td.trk.data(fly,i0:i1,idxnums.ywingl);
  trx(fly).xwingr = td.trk.data(fly,i0:i1,idxnums.xwingr);
  trx(fly).ywingr = td.trk.data(fly,i0:i1,idxnums.ywingr);
    
  % wing angles
  trx(fly).wing_anglel = -td.trk.data(fly,i0:i1,idxnums.wing_anglel);
  trx(fly).wing_angler = td.trk.data(fly,i0:i1,idxnums.wing_angler);
  % remove nans
  islwing = ~isnan(td.trk.data(fly,i0:i1,idxnums.wing_anglel));
  isrwing = ~isnan(td.trk.data(fly,i0:i1,idxnums.wing_angler));
  trx(fly).wing_anglel(~islwing) = 0;
  trx(fly).wing_angler(~isrwing) = 0;
  
  % convert to milimeters
  trx(fly).x_mm = trx(fly).x / pxpermm;
  trx(fly).y_mm = trx(fly).y / pxpermm;
  trx(fly).theta_mm = trx(fly).theta;
  trx(fly).a_mm = trx(fly).a / pxpermm;
  trx(fly).b_mm = trx(fly).b / pxpermm;

  % area
  pd.area_mm{fly} = td.trk.data(fly,i0:i1,idxnums.area) / pxpermm^2;
  
  % always set nwingsdetected to 2
  pd.nwingsdetected{fly} = islwing+isrwing;
  
  % this is really the wing length
  pd.wing_areal{fly} = td.trk.data(fly,i0:i1,idxnums.wing_lengthl)/pxpermm;
  pd.wing_arear{fly} = td.trk.data(fly,i0:i1,idxnums.wing_lengthr)/pxpermm;
  pd.wing_areal{fly}(~islwing) = 0;
  pd.wing_arear{fly}(~isrwing) = 0;
  
  % set trough angle just to be between the two wings
  pd.wing_trough_angle{fly} = (trx(fly).wing_angler - trx(fly).wing_anglel)/2;

  % extra stuff
  pd.fgarea{fly} = td.trk.data(fly,i0:i1,idxnums.fgarea) / pxpermm^2;
  pd.imgcontrast{fly} = td.trk.data(fly,i0:i1,idxnums.imgcontrast);
  pd.minfgdist{fly} = td.trk.data(fly,i0:i1,idxnums.minfgdist);
    
end

msg{end+1} = sprintf('Read trx for %d flies, frame range [%d,%d]',numel(trx),min([trx.firstframe]),max([trx.endframe]));

%% set landmark parameters
arenaradius_mm = arenaradius / pxpermm;
arenacenterx_mm = arenacenterx / pxpermm;
arenacentery_mm = arenacentery / pxpermm;
trx = SetLandmarkParameters(trx,arenatype,arenacenterx_mm,arenacentery_mm,...
  arenaradius_mm,nan,nan); 

if strcmpi(arenatype,'Circle')
  msg{end+1} = sprintf('Set circular arena parameters. Center = (%.1f,%.1f) mm, radius = %.1f mm.',arenacenterx_mm,arenacentery_mm,arenaradius_mm);
elseif strcmpi(arenatype,'Rectangle'),
  msg{end+1} = sprintf('Set rectangular arena parameters. Center = (%.1f,%.1f) mm, width= %.1f mm, height = %.1f mm',arenacenterx_mm,arenacentery_mm,arenawidth_mm,arenaheight_mm);  
end
msg{end+1} = sprintf('x_mm ranges within [%.1f,%.1f], y_mm ranges within [%.1f,%.1f]',...
  min([trx.x_mm]),max([trx.x_mm]),min([trx.y_mm]),max([trx.y_mm]));

% flipud if necessary
if doflipud,
  for fly = 1:nflies,
    trx(fly).y = movieheight - trx(fly).y;
    trx(fly).ywingl = movieheight - trx(fly).ywingl;
    trx(fly).ywingr = movieheight - trx(fly).ywingr;
    trx(fly).theta = -trx(fly).theta;
  end
end

% fliplr if necessary
if dofliplr,
  for fly = 1:nflies,
    trx(fly).x = moviewidth - trx(fly).x;
    trx(fly).xwingl = moviewidth - trx(fly).xwingl;
    trx(fly).xwingr = moviewidth - trx(fly).xwingr;
    trx(fly).theta = modrange(pi-trx(fly).theta,-pi,pi);
  end
end

% transpose if necessary
if dotransposeimage,
  for fly = 1:nflies,
    tmp = trx(fly).x;
    trx(fly).x = trx(fly).y;
    trx(fly).y = tmp;
    c = cos(trx(fly).theta);
    s = sin(trx(fly).theta);
    trx(fly).theta = atan2(c,s);
    tmp = trx(fly).xwingl;
    trx(fly).xwingl = trx(fly).ywingl;
    trx(fly).ywingl = tmp;
    tmp = trx(fly).xwingr;
    trx(fly).xwingr = trx(fly).ywingr;
    trx(fly).ywingr = tmp;
  end
end

if ~isempty(moviefilestr),
  for fly = 1:nflies,
    trx(fly).moviename = moviefile;
  end
end
for fly = 1:nflies,
  trx(fly).matname = trxfile;
  trx(fly).originaltrackfile = intrxfile;
end

%% create the experiment directory
if ~exist(expdir,'dir'),
  [success1,msg1] = mkdir(expdir);
  if ~success1,
    msg = msg1;
    return;
  end
  msg{end+1} = sprintf('Created experiment directory %s.',expdir);
end

%% save the trx file
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

%% copy/soft-link movie

if ~isempty(inmoviefile),
  
  if strcmp(fullfile(inmoviefile),fullfile(moviefile)),
    msg{end+1} = 'Input and out movie files are the same, not copying/linking.';
  else
    
    
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
        inmoviefile = fullfile(inmoviefile);
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
    
  end
  
end

%% make per-frame directory

if ~exist(perframedir,'dir'),
  [success1,msg1] = mkdir(perframedir);
  if ~success1,
    msg = msg1;
    return;
  end
end

%% save per-frame features

fns = fieldnames(pd);
for i = 1:numel(fns),
  fn = fns{i};
  s = struct('data',{pd.(fn)},'units',pfunits.(fn)); %#ok<NASGU>
  filename = fullfile(perframedir,[fn,'.mat']);
  save(filename,'-struct','s');
end

%% done

success = true;