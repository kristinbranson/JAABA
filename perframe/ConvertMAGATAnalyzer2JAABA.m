function [success,msg] = ConvertMAGATAnalyzer2JAABA(varargin)

SCALE = 10; % centimeters instead of millimeters

success = false;
msg = '';

%% set parameters

[inmoviefile,expfile,...
  expdir,moviefilestr,trxfilestr,perframedirstr,...
  arenatype,arenacenterx,arenacentery,...
  arenaradius,arenawidth,arenaheight,...
  pxpermm,fps,overridefps,...
  dosoftlink,dotransposeimage] = myparse(varargin,...
  'inmoviefile','','expfile','',...
  'expdir','','moviefilestr','movie.ufmf','trxfilestr','trx.mat','perframedirstr','perframe',...
  'arenatype','None','arenacenterx',0,'arenacentery',0,...
  'arenaradius',123,'arenawidth',123,'arenaheight',123,...
  'pxpermm',1,'fps',30,...
  'overridefps',false,...
  'dosoftlink',false,...
  'dotransposeimage',false);

%% check for files

if isempty(expfile),
  msg = 'Experiment file is empty';
  return;
end

if isempty(inmoviefile),
  msg = 'Video file is empty';
  return;
end

if ~exist(expfile,'file'),
  msg = sprintf('Experiment file %s does not exist',expfile);
  return;
end

if ~exist(inmoviefile,'file'),
  msg = sprintf('Video file %s does not exist',inmoviefile);
  return;
end

%% load in experiment data

try
  expdata = load(expfile);
catch ME,
  msg = getReport(ME);
  return;
end

%% compute trajectories


if dotransposeimage,
  XIND = 1;
  YIND = 2;
else
  XIND = 2;
  YIND = 1;
end


nflies = numel(expdata.experiment_1.track);
trx = [];
perframedata = struct;

for i = 1:nflies,
  
  fprintf('Larva %d...\n',i);
  
  trk = struct;
  trk.nframes = expdata.experiment_1.track(i).npts;

  % centroid
  loc = double(cat(2,expdata.experiment_1.track(i).pt.loc));
  trk.x_mm = loc(1,:);
  trk.y_mm = loc(2,:);

  % grab other data about the larva
  mid = double(cat(2,expdata.experiment_1.track(i).pt.mid));
  head = double(cat(2,expdata.experiment_1.track(i).pt.head));
  tail = double(cat(2,expdata.experiment_1.track(i).pt.tail));
  spine = {expdata.experiment_1.track(i).pt.spine};
  contour = {expdata.experiment_1.track(i).pt.contour};
  
  
  % convert to pixels, as everything else will be for plotting purposes
  loc_px = cat(1,expdata.experiment_1.camcalinfo.r2cX(loc(1,:),loc(2,:)),...
    expdata.experiment_1.camcalinfo.r2cY(loc(1,:),loc(2,:)));
  trk.x = loc_px(XIND,:)+1;
  trk.y = loc_px(YIND,:)+1;

  % midpoint
  xmid_cm = mid(XIND,:);
  ymid_cm = mid(YIND,:);
  perframedata.xmid_mm{i} = xmid_cm*SCALE;
  perframedata.ymid_mm{i} = ymid_cm*SCALE;
  mid_px = cat(1,expdata.experiment_1.camcalinfo.r2cX(mid(1,:),mid(2,:)),...
    expdata.experiment_1.camcalinfo.r2cY(mid(1,:),mid(2,:)));
  trk.xmid = mid_px(XIND,:)+1;
  trk.ymid = mid_px(YIND,:)+1;

  % head
  xhead_cm = head(XIND,:);
  yhead_cm = head(YIND,:);
  perframedata.xhead_mm{i} = xhead_cm*SCALE;
  perframedata.yhead_mm{i} = yhead_cm*SCALE;
  head_px = cat(1,expdata.experiment_1.camcalinfo.r2cX(head(1,:),head(2,:)),...
    expdata.experiment_1.camcalinfo.r2cY(head(1,:),head(2,:)));
  trk.xhead = head_px(XIND,:)+1;
  trk.yhead = head_px(YIND,:)+1;
  
  % tail
  xtail_cm = tail(XIND,:);
  ytail_cm = tail(YIND,:);
  perframedata.xtail_mm{i} = xtail_cm*SCALE;
  perframedata.ytail_mm{i} = ytail_cm*SCALE;
  tail_px = cat(1,expdata.experiment_1.camcalinfo.r2cX(tail(1,:),tail(2,:)),...
    expdata.experiment_1.camcalinfo.r2cY(tail(1,:),tail(2,:)));
  trk.xtail = tail_px(XIND,:)+1;
  trk.ytail = tail_px(YIND,:)+1;

  % spine
  nspinepts = max(cellfun(@(x) size(x,2),spine));
  isspine = ~cellfun(@isempty,spine);
  xspine_cm = nan(nspinepts,trk.nframes);
  yspine_cm = nan(nspinepts,trk.nframes);
  xspine_cm(:,isspine) = cell2mat(cellfun(@(x) double(x(1,:)'),spine(isspine),'UniformOutput',false));
  yspine_cm(:,isspine) = cell2mat(cellfun(@(x) double(x(2,:)'),spine(isspine),'UniformOutput',false));
  if dotransposeimage,
    trk.xspine = expdata.experiment_1.camcalinfo.r2cY(xspine_cm,yspine_cm)+1;
    trk.yspine = expdata.experiment_1.camcalinfo.r2cX(xspine_cm,yspine_cm)+1;
    trk.xspine_mm = xspine_cm*SCALE;
    trk.yspine_mm = yspine_cm*SCALE;
  else
    trk.yspine = expdata.experiment_1.camcalinfo.r2cX(xspine_cm,yspine_cm)+1;
    trk.xspine = expdata.experiment_1.camcalinfo.r2cY(xspine_cm,yspine_cm)+1;
    trk.xspine_mm = yspine_cm*SCALE;
    trk.yspine_mm = xspine_cm*SCALE;
  end
  
  % contour
  xcontour_cm = cellfun(@(x) double(x(1,:)),contour,'UniformOutput',false);
  ycontour_cm = cellfun(@(x) double(x(2,:)),contour,'UniformOutput',false);
  if dotransposeimage,
    trk.xcontour = cellfun(@(x,y) expdata.experiment_1.camcalinfo.r2cY(x,y)+1, xcontour_cm,ycontour_cm,'UniformOutput',false);
    trk.ycontour = cellfun(@(x,y) expdata.experiment_1.camcalinfo.r2cX(x,y)+1, xcontour_cm,ycontour_cm,'UniformOutput',false);
    perframedata.xcontour_mm{i} = cellfun(@(x) x*SCALE,xcontour_cm,'UniformOutput',false);
    perframedata.ycontour_mm{i} = cellfun(@(x) x*SCALE,ycontour_cm,'UniformOutput',false);
  else
    trk.ycontour = cellfun(@(x,y) expdata.experiment_1.camcalinfo.r2cX(x,y)+1, xcontour_cm,ycontour_cm,'UniformOutput',false);
    trk.xcontour = cellfun(@(x,y) expdata.experiment_1.camcalinfo.r2cY(x,y)+1, xcontour_cm,ycontour_cm,'UniformOutput',false);
    perframedata.xcontour_mm{i} = cellfun(@(x) x*SCALE,ycontour_cm,'UniformOutput',false);
    perframedata.ycontour_mm{i} = cellfun(@(x) x*SCALE,xcontour_cm,'UniformOutput',false);
  end

  % covariance matrix is already in pixels; use it to get major, minor,
  % orientation
  S = cellfun(@(x) [x(1),x(2);x(2),x(3)], {expdata.experiment_1.track(i).pt.cov},'UniformOutput',false);
  S = cat(3,S{:});
  if ~dotransposeimage,
    S = [S(2,2,:),S(2,1,:);S(1,2,:),S(1,1,:)];
  end
  [a,b,theta] = cov2ell(S);
  % note that we don't worry about the zero-indexing correction because the
  % difference will cancel out
  thetahead = atan2(trk.yhead-trk.ytail,trk.xhead-trk.xtail);
  thetaflip = modrange(thetahead+modrange(theta-thetahead,-pi/2,pi/2),-pi,pi);
  trk.a = a/2;
  trk.b = b/2;
  trk.theta = thetaflip;
  
  % convert a, b, theta to mm
  if dotransposeimage,
    x1 = trk.x-1 + a.*cos(thetaflip);
    x2 = trk.x-1 - a.*cos(thetaflip);
    y1 = trk.y-1 + a.*sin(thetaflip);
    y2 = trk.y-1 - a.*sin(thetaflip);
  else
    y1 = trk.x-1 + a.*cos(thetaflip);
    y2 = trk.x-1 - a.*cos(thetaflip);
    x1 = trk.y-1 + a.*sin(thetaflip);
    x2 = trk.y-1 - a.*sin(thetaflip);
  end
  x1_cm = expdata.experiment_1.camcalinfo.c2rX(x1,y1);
  y1_cm = expdata.experiment_1.camcalinfo.c2rY(x1,y1);
  x2_cm = expdata.experiment_1.camcalinfo.c2rX(x2,y2);
  y2_cm = expdata.experiment_1.camcalinfo.c2rY(x2,y2);
  trk.theta_mm = atan2(y2_cm-y1_cm,x2_cm-x1_cm);
  trk.a_mm = sqrt((y2_cm-y1_cm).^2+(x2_cm-x1_cm).^2)/4;
  if dotransposeimage,
    x1 = trk.x-1 + b.*cos(thetaflip+pi/2);
    x2 = trk.x-1 - b.*cos(thetaflip+pi/2);
    y1 = trk.y-1 + b.*sin(thetaflip+pi/2);
    y2 = trk.y-1 - b.*sin(thetaflip+pi/2);
  else
    y1 = trk.x-1 + b.*cos(thetaflip+pi/2);
    y2 = trk.x-1 - b.*cos(thetaflip+pi/2);
    x1 = trk.y-1 + b.*sin(thetaflip+pi/2);
    x2 = trk.y-1 - b.*sin(thetaflip+pi/2);
  end
  x1_cm = expdata.experiment_1.camcalinfo.c2rX(x1,y1);
  y1_cm = expdata.experiment_1.camcalinfo.c2rY(x1,y1);
  x2_cm = expdata.experiment_1.camcalinfo.c2rX(x2,y2);
  y2_cm = expdata.experiment_1.camcalinfo.c2rY(x2,y2);
  trk.b_mm = sqrt((y2_cm-y1_cm).^2+(x2_cm-x1_cm).^2)/4;
    
  % all this was actually in centimeters, so convert to millimeters
  trk.x_mm = trk.x_mm*SCALE;
  trk.y_mm = trk.y_mm*SCALE;
  trk.a_mm = trk.a_mm*SCALE;
  trk.b_mm = trk.b_mm*SCALE;

  % approximate pxpermm
  trk.pxpermm = 1/expdata.experiment_1.camcalinfo.realUnitsPerPixel/SCALE;

  % frame rate
  trk.dt = diff([expdata.experiment_1.track(i).pt.et]);
  
  trk.firstframe = double(expdata.experiment_1.track(i).startFrame+1);
  trk.endframe = trk.firstframe+trk.nframes-1;
  trk.off = 1-trk.firstframe;
  trk.id = expdata.experiment_1.track(i).trackNum;
  trx = structappend(trx,trk);
  
end

%% 

%% create the trx file

outmatname = fullfile(outexpdir,trxfilestr);
save(outmatname,'trx');

%% per-frame data

units = parseunits('unit'); %#ok<NASGU>

perframedir = fullfile(outexpdir,perframedirstr);
if ~exist(perframedir,'dir'),
  mkdir(perframedir);
end

% save the perframedata we computed above
fns = fieldnames(perframedata);
for i = 1:numel(fns),
  fn = fns{i};
  data = perframedata.(fn); %#ok<NASGU>
  matfilename = fullfile(perframedir,[fn,'.mat']);
  fprintf('Saving data to %s...\n',matfilename);
  save(matfilename,'data','units');
end

% for some reason, the derived measurements have different indices
idx = cell(1,nflies);
for i = 1:nflies,
  idx{i} = expdata.experiment_1.track(i).getDerivedQuantity('mapptstointerped',false);
end

fielddict = struct;
fielddict.vel = 'vel_vel';
fielddict.speed = 'vel_speed';
fielddict.vnorm = 'vel_vnorm';
fielddict.theta = 'vel_theta';
fielddict.adjspeed = 'vel_adjspeed';
%fielddict.speed_diff_local = 'speed_diff_local';
fielddict.deltatheta = 'acc_deltatheta';
fielddict.ddtheta = 'acc_ddtheta';
fielddict.curv = 'acc_curv';
fielddict.lrdtheta = 'lrdtheta';
fielddict.pathLength = 'pathLength';
fielddict.covRatio = 'cov_covRatio';
fielddict.covTheta = 'cov_covTheta';
fielddict.covMinor = 'cov_covMinor';
fielddict.covMajor = 'cov_covMajor';
fielddict.dcovRatio = 'dcovRatio';
fielddict.sarea = 'sarea';

fnsin = fieldnames(fielddict);

for j = 1:numel(fnsin),
  fnin = fnsin{j};
  fnout = fielddict.(fnin);
  matfilename = fullfile(perframedir,[fnout,'.mat']);
  fprintf('%s -> %s...\n',fnin,matfilename);

  data = cell(1,nflies);
  
  for i = 1:nflies,
    data{i} = expdata.experiment_1.track(i).getDerivedQuantity(fnin);    
  end
  save(matfilename,'data','units');  
end

% also save unsmoothed area
data = cell(1,nflies);
for i = 1:nflies,
  data{i} = [expdata.experiment_1.track(i).pt.area];
end
matfilename = fullfile(perframedir,['area','.mat']);
save(matfilename,'data','units');

%% copy over the mmf file

outmoviefile = fullfile(outexpdir,moviefilestr);
if isunix && makesoftlink,
  if exist(outmoviefile,'file'),
    delete(outmoviefile);
  end
  cmd = sprintf('ln -s %s %s',inmoviefile,outmoviefile);
  unix(cmd);
else
  [success,msg] = copyfile(inmoviefile,outmoviefile);
  if ~success,
    error('Error copying file %s to %s: %s',inmoviefile,outmoviefile,msg);
  end
end

