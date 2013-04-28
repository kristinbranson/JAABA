function [success,msg] = ConvertLarvaeLouis2JAABA(varargin)

success = false;
msg = {};

%% constants

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


%% parse parameters

[inmoviefile,indatafile,inkinmatfile,expi,...
  expdir,moviefilestr,trxfilestr,perframedirstr,...
  arenatype,arenacenterx,arenacentery,...
  arenaradius,arenawidth,arenaheight,...
  pxpermm,fps,overridefps,overridearena,...
  dosoftlink,ncontourpts] = myparse(varargin,...
  'inmoviefile','','indatafile','','inkinmatfile','','expi',[],...
  'expdir','','moviefilestr','movie.avi','trxfilestr','trx.mat','perframedirstr','perframe',...
  'arenatype','None','arenacenterx',0,'arenacentery',0,...
  'arenaradius',123,'arenawidth',123,'arenaheight',123,...
  'pxpermm',1,'fps',30,...
  'overridefps',false,'overridearena',false,...
  'dosoftlink',false,...
  'ncontourpts',20); %#ok<ASGLU>

% get experiment_name from full path
[~,experiment_name] = myfileparts(expdir); %#ok<NASGU>

% check that required inputs exist
ismovie = ~isempty(inmoviefile);
if ismovie && ~exist(inmoviefile,'file'),
  msg = sprintf('Input movie file %s does not exist',inmoviefile);
  return;
end
if ~exist(indatafile,'file'),
  msg = sprintf('Input data file %s does not exist',indatafile);
  return;
end
if ~exist(inkinmatfile','file'),
  msg = sprintf('Input kinData file %s does not exist',inkinmatfile);
  return;
end

%% output file names

trxfile = fullfile(expdir,trxfilestr);
moviefile = fullfile(expdir,moviefilestr);
perframedir = fullfile(expdir,perframedirstr);

%% convert trx

%indata = importdata(indataname);
indata.data = [];
fid = fopen(indatafile,'r');
if fid < 0,
  error('Could not open file %s for reading',indatafile);
end
while true,
  s = fgetl(fid);
  if ~ischar(s),
    break;
  end
  s = strtrim(s);
  if isempty(s),
    continue;
  end
  ss = str2double(strsplit(s,','));
  indata.data(:,end+1) = ss;
end
%indata.data = [(0:size(indata.data,1)-1)',indata.data];
fclose(fid);

% % read in relevant parameters from config file
% fid = fopen(inconfigfile,'r');
% if ~ismovie,
%   pxpermm = 1;
% end
% readpxpermm = overridearena || ~ismovie;
% readtickpermm = false;
% while true,
%   s = fgetl(fid);
%   if ~ischar(s),
%     break;
%   end
%   s = strtrim(s);
%   if isempty(s),
%     continue;
%   end
%   if ~readpxpermm,
%     m = regexp(s,'^Camera Calibration \(um per pixel\): (.*)$','tokens','once');
%     if ~isempty(m),
%       pxpermm = 1000/str2double(m{1});
%       readpxpermm = true;
%     end
%   end
%   if ~readtickpermm,
%     m = regexp(s,'^Stage Calibration \(tick per mm\): x = (.*), y= (.*)$','tokens','once');
%     if ~isempty(m),
%       xtickpermm = str2double(m{1});
%       ytickpermm = str2double(m{2});
%       readtickpermm = true;
%     end
%   end
% end
% fclose(fid);
% if ~readtickpermm,
%   msg = sprintf('Could not read stage ticks per mm from %s',inconfigfile);
%   return;
% end
% msg{end+1} = sprintf('pxpermm = %f, x-ticks per mm = %f, y-ticks per mm = %f',pxpermm,xtickpermm,ytickpermm);

timestamps = indata.data(2,:);

nframes = size(indata.data,2);

% load in kinData
try
  load(inkinmatfile,'kinData');
catch %#ok<CTCH>
  msg = sprintf('Could not load kinData from file %s',inkinmatfile);
  return;
end

canmatch = false(1,numel(kinData)); %#ok<NODEF>
for i = 1:numel(kinData),
  if isfield(kinData{i},'centroidposition'),
    canmatch(i) = nframes == size(kinData{i}.centroidposition,1);
  end
end

if ~any(canmatch),
  msg = sprintf('No cells in kinData have the right number of frames (%d)',nframes);
  return;
end  

if isempty(expi),
  
  if nnz(canmatch) == 1,
    expi = find(canmatch);
    msg{end+1} = sprintf('Based on number of frames, setting expi = %d',expi);
  else
    idxmatch = find(canmatch(:));
    [sel,ok] = listdlg('Which experiment is this within the kinData file?',cellstr(num2str(idxmatch)));
    if ~ok,
      msg = 'No experiment index selected';
      return;
    end
    expi = idxmatch(sel);
    msg{end+1} = sprintf('Selected experiment index %d',expi);
  end
end

kinData = kinData{expi};

trx = struct;

if ismovie,
  trx.x = indata.data(10,:);
  trx.y = indata.data(11,:);
  trx.xspine = indata.data([3,5,7],:);
  trx.yspine = indata.data([4,6,8],:);
  trx.a = indata.data(9,:) / 4;
  trx.b = trx.a;
  trx.theta = indata.data(13,:)*pi/180;
else
  trx.x = kinData.centroidposition(:,1)';
  trx.y = kinData.centroidposition(:,2)';
  trx.xspine = [kinData.headposition(:,1)';kinData.centroidposition(:,1)';kinData.tailposition(:,1)'];
  trx.yspine = [kinData.headposition(:,2)';kinData.centroidposition(:,2)';kinData.tailposition(:,2)'];
  trx.a = kinData.skeL' / 4;
  trx.b = trx.a;
  trx.theta = kinData.bodyangle'*pi/180;
end
trx.dt = diff(timestamps);
trx.id = 1; % only one larva
trx.moviename = moviefile;
trx.firstframe = 1;
trx.endframe = numel(trx.x);
trx.nframes = numel(trx.x);
trx.off = 0;
trx.timestamps = timestamps;
trx.dt = diff(timestamps);

stagex = indata.data(14,:);
stagey = indata.data(15,:);

% regress to find scaling from pixels to mm
[xcoeffs,~,resid] = regress(trx.x(:),[indata.data(14,:)',indata.data(10,:)',ones(nframes,1)]);
maxerr = max(abs(resid));
if maxerr > .01,
  msg = sprintf('Could not find a low error regression from pixels to mms for x, maximum residual = %f with x_mm = %f*stagex + %f*x_px + %f',maxerr,xcoeffs(1),xcoeffs(2),xcoeffs(3));
  return;
end
msg{end+1} = sprintf('x-coordinate px to mm regression has error %f: x_mm = %f*stagex + %f*x_px + %f',maxerr,xcoeffs(1),xcoeffs(2),xcoeffs(3));

[ycoeffs,~,resid] = regress(trx.y(:),[indata.data(15,:)',indata.data(11,:)',ones(nframes,1)]);
maxerr = max(abs(resid));
if maxerr > .01,
  msg = sprintf('Could not find a low error regression from pixels to mms for y, maximum residual = %f with y_mm = %f*stagey + %f*y_px + %f',maxerr,ycoeffs(1),ycoeffs(2),ycoeffs(3));
  return;
end
msg{end+1} = sprintf('y-coordinate px to mm regression has error %f: y_mm = %f*stagey + %f*y_px + %f',maxerr,ycoeffs(1),ycoeffs(2),ycoeffs(3));

pxpermm = (xcoeffs(2)+ycoeffs(2))/2;

if ismovie,
  contours_px = getContoursLarvaeLouis(indata.data(19:end,:)',ncontourpts);
else
  contours_px0 = getContoursLarvaeLouis(indata.data(19:end,:)',ncontourpts);
  contours_px = nan(size(contours_px0));
  contours_px(:,:,1) = bsxfun(@plus,indata.data(14,:)'*xcoeffs(1),contours_px0(:,:,1)*xcoeffs(2)) + xcoeffs(3);
  contours_px(:,:,2) = bsxfun(@plus,indata.data(15,:)'*ycoeffs(1),contours_px0(:,:,2)*ycoeffs(2)) + ycoeffs(3);
  
end

trx.x_mm = kinData.centroidposition(:,1)';
trx.y_mm = kinData.centroidposition(:,2)';
trx.a_mm = kinData.skeL';
trx.b_mm = kinData.skeL';
trx.theta_mm = trx.theta;

trx.xcontour = mat2cell(contours_px(:,:,1)',ncontourpts,ones(1,trx.nframes));
trx.ycontour = mat2cell(contours_px(:,:,2)',ncontourpts,ones(1,trx.nframes));

%% convert to mm

if ~strcmpi(arenatype,'None'),
  msg{end+1} = sprintf('Arena type input as %s, but must be None for this type of data. Setting to None.',arenatype);
end

arenacenterx_mm = mean(trx.x_mm);
arenacentery_mm = mean(trx.y_mm);
arenatype = 'None';

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
  [success1,msg1] = mkdir(expdir);
  if ~success1,
    msg = msg1;
    return;
  end
end

%% create per-frame directory

if ~exist(perframedir,'dir'),
  [success1,msg1] = mkdir(perframedir);
  if ~success1,
    msg = msg1;
    return;
  end
end

%% save the per-frame data

perframefns1 = {
  'skeL'
  'led'
  'ledSpeed'
  'ledSpeedRaw'
  'headSpeed'
  'tailSpeed'
  'centroidSpeed'
  'midSpeed'
  'bodyangle'
  'headangle'
  'bodyangleSpeed'
  'headangleSpeed'
  'mode'
  };
perframefns2 = {
  'headposition'
  'tailposition'
  'centroidposition'
  'midposition'
  };

for i = 1:numel(perframefns1),
  fn = perframefns1{i};
  data = {kinData.(fn)(:)'}; %#ok<NASGU>
  units = parseunits('unit'); %#ok<NASGU>
  outfile = fullfile(perframedir,[fn,'.mat']);
  try
    save(outfile,'data','units');
  catch ME,
    msg = getReport(ME);
    return;
  end
end
for i = 1:numel(perframefns2),
  fn = perframefns2{i};
  data1 = kinData.(fn); 
  data = {data1(:,1)'}; %#ok<NASGU>
  units = parseunits('mm'); %#ok<NASGU>
  outfn = [fn,'_x_mm'];
  outfile = fullfile(perframedir,[outfn,'.mat']);
  try
    save(outfile,'data','units');
  catch ME,
    msg = getReport(ME);
    return;
  end
  data = {data1(:,2)'}; %#ok<NASGU>
  units = parseunits('mm'); %#ok<NASGU>
  outfn = [fn,'_y_mm'];
  outfile = fullfile(perframedir,[outfn,'.mat']);
  try
    save(outfile,'data','units');
  catch ME,
    msg = getReport(ME);
    return;
  end
end

data = {stagex}; %#ok<NASGU>
units = parseunits('unit'); %#ok<NASGU>
outfile = fullfile(perframedir,'stagex.mat');
try
  save(outfile,'data','units');
catch ME,
  msg = getReport(ME);
  return;
end
data = {stagey}; %#ok<NASGU>
outfile = fullfile(perframedir,'stagey.mat');
try
  save(outfile,'data','units');
catch ME,
  msg = getReport(ME);
  return;
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

if ismovie,
  
  if strcmp(fullfile(inmoviefile),fullfile(moviefile)),
    fprintf('Input and out movie files are the same, not copying/linking.\n');
  else
    
    if dosoftlink,
      if isunix,
        
        % delete existing files
        if exist(inmoviefile,'file') && strcmp(inmoviefile,moviefile),
        elseif exist(inmoviefile,'file'),
          [status,result] = unix(sprintf('readlink %s',moviefile));
          result = strtrim(result);
          if status == 0,
            msg{end+1} = sprintf('Softlink to %s existed, deleted',result);
            unix(sprintf('rm %s',moviefile));
          else
            res = questdlg(sprintf('File %s exists. Delete and overwrite with softlink to %s?',moviefile,inmoviefile),'Overwrite?','Yes','No','Cancel','Cancel');
            if strcmpi(res,'No'),
              msg = sprintf('File %s already existed',moviefile);
              return;
            elseif strcmpi(res,'Yes'),
              msg{end+1} = sprintf('File %s existed and was deleted',moviefile);
              delete(moviefile);
            end
          end
        end

        if ~exist(moviefile,'file'),
          cmd = sprintf('ln -s %s %s',inmoviefile,moviefile);
          unix(cmd);
          % test to make sure it worked
          [status,result] = unix(sprintf('readlink %s',moviefile));
          result = strtrim(result);
          if status ~= 0 || ~strcmp(result,inmoviefile),
            warndlg(sprintf('Failed to make soft link, copying %s to %s instead',inmoviefile,moviefile));
            dosoftlink = false;
          end
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
