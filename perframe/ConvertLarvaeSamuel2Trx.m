function [outexpdir] = ConvertLarvaeSamuel2Trx(inexpdir,rootoutputdir,varargin)

%% set parameters

% default parameters
trxfilestr = 'trx.mat';
moviefilestr = 'movie.mmf';
perframedirstr = 'perframe';
DEBUG = false;
makesoftlink = true;


[trxfilestr,moviefilestr,perframedirstr,DEBUG,makesoftlink] = ...
  myparse(varargin,...
  'trxfilestr',trxfilestr,...
  'moviefilestr',moviefilestr,...
  'perframedirstr',perframedirstr,...
  'debug',DEBUG,...
  'makesoftlink',makesoftlink); %#ok<ASGLU>

%% parse experiment directory name
[rootdatadir,inexperiment_name] = fileparts(inexpdir); %#ok<ASGLU>
timestamp = inexperiment_name;

% find experiment file
files = dir(fullfile(inexpdir,['*experiment_',timestamp,'*.mat']));
if isempty(files),
  error('Could not find experiment file in directory %s',inexpdir);
end
for i = 1:numel(files),
  metadata = regexp(files(i).name,['^(?<line_name>.+)@(?<effector>.+)_experiment_',timestamp,'\.mat$'],'names','once');
  if ~isempty(metadata),
    experimentfile = fullfile(inexpdir,files(i).name);
    break;
  end
end
if isempty(metadata),
  error('Could not find experiment file in directory %s',inexpdir);
end
metadata.exp_datetime = timestamp;

%% find movie file
moviefilepattern = sprintf('%s@%s@%s*.mmf',metadata.exp_datetime,metadata.line_name,metadata.effector);
files = dir(fullfile(inexpdir,moviefilepattern));
if isempty(files),
  error('Could not find movie file in %s',inexpdir);
end
inmoviefile = fullfile(inexpdir,files(1).name);

%% create output directory

% make sure the root output directory exists
if ~exist(rootoutputdir,'dir'),
  mkdir(rootoutputdir);
end

% create the experiment directory
outexperiment_name = sprintf('%s_%s_%s',metadata.line_name,metadata.effector,metadata.exp_datetime);
outexpdir = fullfile(rootoutputdir,outexperiment_name);
if ~exist(outexpdir,'dir'),
  mkdir(outexpdir);
end

%% load in experiment data
expdata = load(experimentfile);

%% compute trajectories
nflies = numel(expdata.experiment_1.track);
trx = [];
perframedata = struct;
SCALE = 10; % centimeters instead of millimeters

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
  
  % approximate pxpermm
  trk.pxpermm = 1/expdata.experiment_1.camcalinfo.realUnitsPerPixel;
  
  % convert to pixels, as everything else will be for plotting purposes
  trk.x = expdata.experiment_1.camcalinfo.r2cX(trk.x_mm,trk.y_mm);
  trk.y = expdata.experiment_1.camcalinfo.r2cY(trk.x_mm,trk.y_mm);

  % midpoint
  xmid_mm = mid(1,:);
  ymid_mm = mid(2,:);
  trk.xmid = expdata.experiment_1.camcalinfo.r2cX(xmid_mm,ymid_mm);
  trk.ymid = expdata.experiment_1.camcalinfo.r2cY(xmid_mm,ymid_mm);

  % head
  xhead_mm = head(1,:);
  yhead_mm = head(2,:);
  trk.xhead = expdata.experiment_1.camcalinfo.r2cX(xhead_mm,yhead_mm);
  trk.yhead = expdata.experiment_1.camcalinfo.r2cY(xhead_mm,yhead_mm);
  
  % tail
  xtail_mm = tail(1,:);
  ytail_mm = tail(2,:);
  trk.xtail = expdata.experiment_1.camcalinfo.r2cX(xtail_mm,ytail_mm);
  trk.ytail = expdata.experiment_1.camcalinfo.r2cY(xtail_mm,ytail_mm);

  % spine
  nspinepts = max(cellfun(@(x) size(x,2),spine));
  isspine = ~cellfun(@isempty,spine);
  xspine_mm = nan(nspinepts,trk.nframes);
  yspine_mm = nan(nspinepts,trk.nframes);
  xspine_mm(:,isspine) = cell2mat(cellfun(@(x) double(x(1,:)'),spine(isspine),'UniformOutput',false));
  yspine_mm(:,isspine) = cell2mat(cellfun(@(x) double(x(2,:)'),spine(isspine),'UniformOutput',false));
  trk.xspine = expdata.experiment_1.camcalinfo.r2cX(xspine_mm,yspine_mm);
  trk.yspine = expdata.experiment_1.camcalinfo.r2cY(xspine_mm,yspine_mm);

  % contour
  xcontour_mm = cellfun(@(x) double(x(1,:)),contour,'UniformOutput',false);
  ycontour_mm = cellfun(@(x) double(x(2,:)),contour,'UniformOutput',false);
  trk.xcontour = cellfun(@(x,y) expdata.experiment_1.camcalinfo.r2cX(x,y), xcontour_mm,ycontour_mm,'UniformOutput',false);
  trk.ycontour = cellfun(@(x,y) expdata.experiment_1.camcalinfo.r2cY(x,y), xcontour_mm,ycontour_mm,'UniformOutput',false);

  % covariance matrix is already in pixels; use it to get major, minor,
  % orientation
  S = cellfun(@(x) [x(1),x(2);x(2),x(3)], {expdata.experiment_1.track(i).pt.cov},'UniformOutput',false);
  S = cat(3,S{:});
  [a,b,theta] = cov2ell(S);
  thetahead = atan2(trk.yhead-trk.ytail,trk.xhead-trk.xtail);
  thetaflip = modrange(thetahead+modrange(theta-thetahead,-pi/2,pi/2),-pi,pi);
  trk.a = a/2;
  trk.b = b/2;
  trk.theta = thetaflip;
 
  % convert a, b, theta to mm
  x1 = trk.x + a.*cos(thetaflip);
  x2 = trk.x - a.*cos(thetaflip);
  y1 = trk.y + a.*sin(thetaflip);
  y2 = trk.y - a.*sin(thetaflip);
  x1_mm = expdata.experiment_1.camcalinfo.c2rX(x1,y1);
  y1_mm = expdata.experiment_1.camcalinfo.c2rY(x1,y1);
  x2_mm = expdata.experiment_1.camcalinfo.c2rX(x1,y1);
  y2_mm = expdata.experiment_1.camcalinfo.c2rY(x2,y2);
  trk.theta_mm = atan2(y2_mm-y1_mm,x2_mm-x1_mm);
  trk.a_mm = sqrt((y2_mm-y1_mm).^2+(x2_mm-x1_mm).^2)/4;
  x1 = trk.x + b.*cos(thetaflip+pi/2);
  x2 = trk.x - b.*cos(thetaflip+pi/2);
  y1 = trk.y + b.*sin(thetaflip+pi/2);
  y2 = trk.y - b.*sin(thetaflip+pi/2);
  x1_mm = expdata.experiment_1.camcalinfo.c2rX(x1,y1);
  y1_mm = expdata.experiment_1.camcalinfo.c2rY(x1,y1);
  x2_mm = expdata.experiment_1.camcalinfo.c2rX(x1,y1);
  y2_mm = expdata.experiment_1.camcalinfo.c2rY(x2,y2);
  trk.b_mm = sqrt((y2_mm-y1_mm).^2+(x2_mm-x1_mm).^2)/4;
    
  % all this was actually in centimeters, so convert to millimeters
  trk.x_mm = trk.x_mm*SCALE;
  trk.y_mm = trk.y_mm*SCALE;
  trk.a_mm = trk.a_mm*SCALE;
  trk.b_mm = trk.b_mm*SCALE;
  trk.pxpermm = trk.pxpermm/SCALE;

  % frame rate
  trk.dt = diff([expdata.experiment_1.track(i).pt.et]);
  
  trk.firstframe = double(expdata.experiment_1.track(i).startFrame+1);
  trk.endframe = trk.firstframe+trk.nframes-1;
  trk.off = 1-trk.firstframe;
  trk.id = expdata.experiment_1.track(i).trackNum;
  trx = structappend(trx,trk);

  % TODO: compute some useful per-frame data
  
  % speed of each body point
  perframedata.headspeed{i} = sqrt(diff(yhead_mm).^2+diff(xhead_mm).^2)*SCALE./trk.dt;
  perframedata.midspeed{i} = sqrt(diff(ymid_mm).^2+diff(xmid_mm).^2)*SCALE./trk.dt;
  perframedata.tailspeed{i} = sqrt(diff(ytail_mm).^2+diff(xtail_mm).^2)*SCALE./trk.dt;

  % head direction
  perframedata.headdir{i} = atan2(yhead_mm-ymid_mm,xhead_mm-xmid_mm);
  perframedata.dheaddir{i} = diff(perframedata.headdir{i})./trk.dt;
  
  % tail to mid direction
  perframedata.taildir{i} = atan2(ymid_mm-ytail_mm,xmid_mm-xtail_mm);
  perframedata.dtaildir{i} = diff(perframedata.taildir{i})./trk.dt;
  
  % bend
  perframedata.bend{i} = modrange(perframedata.headdir{i}-perframedata.taildir{i},-pi,pi);
  
  % change in bend angle
  perframedata.dbend{i} = modrange(diff(perframedata.bend{i}),-pi,pi);
  
  % skeleton length
  perframedata.spinelength{i} = sum(sqrt(diff(yspine_mm).^2+diff(xspine_mm).^2))*SCALE;

  % change in skeleton length
  perframedata.dspinelength{i} = diff(perframedata.spinelength{i}) ./ trk.dt;
  
  % velocity relative to head direction
  perframedata.headvelforward_headdir{i} = (diff(xhead_mm).*cos(perframedata.headdir{i}(1:end-1)) + ...
    diff(yhead_mm).*sin(perframedata.headdir{i}(1:end-1))) * SCALE ./ trk.dt;
  perframedata.headvelsideways_headdir{i} = (diff(xhead_mm).*cos(perframedata.headdir{i}(1:end-1)+pi/2) + ...
    diff(yhead_mm).*sin(perframedata.headdir{i}(1:end-1))+pi/2) * SCALE ./ trk.dt;
  perframedata.midvelforward_headdir{i} = (diff(xmid_mm).*cos(perframedata.headdir{i}(1:end-1)) + ...
    diff(ymid_mm).*sin(perframedata.headdir{i}(1:end-1))) * SCALE ./ trk.dt;
  perframedata.midvelsideways_headdir{i} = (diff(xmid_mm).*cos(perframedata.headdir{i}(1:end-1)+pi/2) + ...
    diff(ymid_mm).*sin(perframedata.headdir{i}(1:end-1))+pi/2) * SCALE ./ trk.dt;
  perframedata.tailvelforward_headdir{i} = (diff(xtail_mm).*cos(perframedata.headdir{i}(1:end-1)) + ...
    diff(ytail_mm).*sin(perframedata.headdir{i}(1:end-1))) * SCALE ./ trk.dt;
  perframedata.tailvelsideways_headdir{i} = (diff(xtail_mm).*cos(perframedata.headdir{i}(1:end-1)+pi/2) + ...
    diff(ytail_mm).*sin(perframedata.headdir{i}(1:end-1))+pi/2) * SCALE ./ trk.dt;
  perframedata.centroidvelforward_headdir{i} = (diff(trk.x_mm).*cos(perframedata.headdir{i}(1:end-1)) + ...
    diff(trk.y_mm).*sin(perframedata.headdir{i}(1:end-1))) ./ trk.dt;
  perframedata.centroidvelsideways_headdir{i} = (diff(trk.x_mm).*cos(perframedata.headdir{i}(1:end-1)+pi/2) + ...
    diff(trk.y_mm).*sin(perframedata.headdir{i}(1:end-1))+pi/2) ./ trk.dt;

  % velocity relative to tail dir
  perframedata.headvelforward_taildir{i} = (diff(xhead_mm).*cos(perframedata.taildir{i}(1:end-1)) + ...
    diff(yhead_mm).*sin(perframedata.taildir{i}(1:end-1))) * SCALE ./ trk.dt;
  perframedata.headvelsideways_taildir{i} = (diff(xhead_mm).*cos(perframedata.taildir{i}(1:end-1)+pi/2) + ...
    diff(yhead_mm).*sin(perframedata.taildir{i}(1:end-1))+pi/2) * SCALE ./ trk.dt;
  perframedata.midvelforward_taildir{i} = (diff(xmid_mm).*cos(perframedata.taildir{i}(1:end-1)) + ...
    diff(ymid_mm).*sin(perframedata.taildir{i}(1:end-1))) * SCALE ./ trk.dt;
  perframedata.midvelsideways_taildir{i} = (diff(xmid_mm).*cos(perframedata.taildir{i}(1:end-1)+pi/2) + ...
    diff(ymid_mm).*sin(perframedata.taildir{i}(1:end-1))+pi/2) * SCALE ./ trk.dt;
  perframedata.tailvelforward_taildir{i} = (diff(xtail_mm).*cos(perframedata.taildir{i}(1:end-1)) + ...
    diff(ytail_mm).*sin(perframedata.taildir{i}(1:end-1))) * SCALE ./ trk.dt;
  perframedata.tailvelsideways_taildir{i} = (diff(xtail_mm).*cos(perframedata.taildir{i}(1:end-1)+pi/2) + ...
    diff(ytail_mm).*sin(perframedata.taildir{i}(1:end-1))+pi/2) * SCALE ./ trk.dt;
  perframedata.centroidvelforward_taildir{i} = (diff(trk.x_mm).*cos(perframedata.taildir{i}(1:end-1)) + ...
    diff(trk.y_mm).*sin(perframedata.taildir{i}(1:end-1))) ./ trk.dt;
  perframedata.centroidvelsideways_taildir{i} = (diff(trk.x_mm).*cos(perframedata.taildir{i}(1:end-1)+pi/2) + ...
    diff(trk.y_mm).*sin(perframedata.taildir{i}(1:end-1))+pi/2) ./ trk.dt;

  % velocity relative to ellipse orientation
  perframedata.headvelforward_bodydir{i} = (diff(xhead_mm).*cos(trk.theta_mm(1:end-1)) + ...
    diff(yhead_mm).*sin(trk.theta_mm(1:end-1))) * SCALE ./ trk.dt;
  perframedata.headvelsideways_bodydir{i} = (diff(xhead_mm).*cos(trk.theta_mm(1:end-1)+pi/2) + ...
    diff(yhead_mm).*sin(trk.theta_mm(1:end-1))+pi/2) * SCALE ./ trk.dt;
  perframedata.midvelforward_bodydir{i} = (diff(xmid_mm).*cos(trk.theta_mm(1:end-1)) + ...
    diff(ymid_mm).*sin(trk.theta_mm(1:end-1))) * SCALE ./ trk.dt;
  perframedata.midvelsideways_bodydir{i} = (diff(xmid_mm).*cos(trk.theta_mm(1:end-1)+pi/2) + ...
    diff(ymid_mm).*sin(trk.theta_mm(1:end-1))+pi/2) * SCALE ./ trk.dt;
  perframedata.tailvelforward_bodydir{i} = (diff(xtail_mm).*cos(trk.theta_mm(1:end-1)) + ...
    diff(ytail_mm).*sin(trk.theta_mm(1:end-1))) * SCALE ./ trk.dt;
  perframedata.tailvelsideways_bodydir{i} = (diff(xtail_mm).*cos(trk.theta_mm(1:end-1)+pi/2) + ...
    diff(ytail_mm).*sin(trk.theta_mm(1:end-1))+pi/2) * SCALE ./ trk.dt;
  perframedata.centroidvelforward_bodydir{i} = (diff(trk.x_mm).*cos(trk.theta_mm(1:end-1)) + ...
    diff(trk.y_mm).*sin(trk.theta_mm(1:end-1))) ./ trk.dt;
  perframedata.centroidvelsideways_bodydir{i} = (diff(trk.x_mm).*cos(trk.theta_mm(1:end-1)+pi/2) + ...
    diff(trk.y_mm).*sin(trk.theta_mm(1:end-1))+pi/2) ./ trk.dt;
  
end

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

