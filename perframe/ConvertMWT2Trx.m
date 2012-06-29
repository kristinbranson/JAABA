function [outexpdir] = ConvertMWT2Trx(expdir,rootoutputdir,varargin)

%% set parameters

% default parameters
blobsext = 'blobs';
datext = 'dat';
spineext = 'spine';
trxfilestr = 'trx.mat';
perframedirstr = 'perframe';
MINDT = .001;
MAXTIMESTAMPERR = .01;
DEBUG = false;
makesoftlink = true;

perframedict = {'x','x_mm'
  'y','y_mm'
  'area','area_mm'};
trxfns = {'x','y','a','b','theta','x_mm','y_mm','a_mm','b_mm','theta_mm'};

[~,experiment_name] = myfileparts(expdir);

[blobsext,datext,spineext,trxfilestr,perframedirstr,...
  DEBUG,makesoftlink,perframedict,trxfns] = ...
  myparse(varargin,...
  'blobsext',blobsext,...
  'datext',datext,...
  'spineext',spineext,...
  'trxfilestr',trxfilestr,...
  'perframedirstr',perframedirstr,...
  'debug',DEBUG,...
  'makesoftlink',makesoftlink,...
  'perframedict',perframedict,...
  'trxfns',trxfns);

%% read in blobs files

% find blobs files
files = dir(fullfile(expdir,[experiment_name,'*.',blobsext]));
blobsnames = cellfun(@(x) fullfile(expdir,x),{files.name},'UniformOutput',false);

if isempty(blobsnames),
  error('Could not find any blobs files in %s',expdir);
end
fprintf('Reading %d blobs files...\n',numel(blobsnames));
trx = ReadMWTBlobFile(blobsnames);
nflies = numel(trx);

% make sure all the timestamps agree with each other
timestamps = unique([trx.timestamp]);
dts = diff(timestamps);
% sanity check to look for rounding errors
if any(dts < MINDT),
  error('Rounding error check failed: there are timestamps < %f apart',MINDT);
end

% add dt, fps
fps = nanmean(dts);
for i = 1:nflies,
  trx(i).dt = diff(trx(i).timestamp);
  trx(i).fps = fps;
end

%% read in spine files

haveremoved = false;
removeid = false(1,numel(trx));
trxids = [trx.id];

% find dat files
files = dir(fullfile(expdir,[experiment_name,'*.',spineext]));
if isempty(files),
  warning('Spine file does not exist in directory %s',expdir);
else
  spinename = fullfile(expdir,files(1).name);

  fprintf('Reading spine file %s...\n',spinename);
  
  % read in data
  datacurr = ReadChoreSpineFile(spinename);
  
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
    warning('ids found in spine file %s that are not in trx',spinename);
    idx(~oldid) = [];
  end
  for jj = 1:numel(idx),
    j = idx(jj);
    [~,f0] = min(abs(trx(j).timestamp-datacurr(jj).timestamp(1)));
    [~,f1] = min(abs(trx(j).timestamp-datacurr(jj).timestamp(end)));
    if f1-f0+1 ~= size(datacurr(jj).xspine,2),
      error('timestamps don''t match up for spine');
    end
    maxerr = max(abs(trx(j).timestamp(f0:f1)-datacurr(jj).timestamp));
    if maxerr > MAXTIMESTAMPERR,
      error('timestamps don''t match up for spine');
    end
    trx(j).xspine_mm = nan(size(datacurr(jj).xspine,1),trx(j).nframes);
    trx(j).yspine_mm = nan(size(datacurr(jj).xspine,1),trx(j).nframes);
    trx(j).xspine_mm(:,f0:f1) = datacurr(jj).xspine;
    trx(j).yspine_mm(:,f0:f1) = datacurr(jj).yspine;
  end
  
end

%% read in dat files

% find dat files
files = dir(fullfile(expdir,[experiment_name,'*.',datext]));
datnames = {};
perframefns = cell(0,2);

for i = 1:numel(files),
  match = regexp(files(i).name,['\.(.+)\.',datext,'$'],'tokens','once');
  if isempty(match),
    continue;
  end
  datnames{end+1} = fullfile(expdir,files(i).name); %#ok<AGROW>
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
  perframefns(end+1,:) = {fn,perframefn}; %#ok<AGROW>
  
  fprintf('Reading in dat file %s...\n',datnames{end});
  
  % read in data
  datacurr = ReadChoreDatFile(datnames{end});
  
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
    warning('ids found in dat file %s that are not in trx',datnames{end});
    idx(~oldid) = [];
  end
  for jj = 1:numel(idx),
    j = idx(jj);
    [~,f0] = min(abs(trx(j).timestamp-datacurr(jj).timestamp(1)));
    [~,f1] = min(abs(trx(j).timestamp-datacurr(jj).timestamp(end)));
    if f1-f0+1 ~= numel(datacurr(jj).value),
      error('timestamps don''t match up for %s',fn);
    end
    maxerr = max(abs(trx(j).timestamp(f0:f1)-datacurr(jj).timestamp));
    if maxerr > MAXTIMESTAMPERR,
      error('timestamps don''t match up for %s',fn);
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

%% compute pxpermm

pxpermm = nanmean([trx.x trx.y] ./ [trx.x_mm trx.y_mm]);
for i = 1:numel(trx),
  trx(i).pxpermm = pxpermm;
  if ~isfield(trx,'x_mm'),
    trx(i).x_mm = trx(i).x / pxpermm;
  end
  if ~isfield(trx,'y_mm'),
    trx(i).y_mm = trx(i).y / pxpermm;
  end
  trx(i).a_mm = trx(i).a / pxpermm;
  trx(i).b_mm = trx(i).b / pxpermm;
  trx(i).theta_mm = trx(i).theta;
  trx(i).xspine = trx(i).xspine_mm * pxpermm;
  trx(i).yspine = trx(i).yspine_mm * pxpermm;
end


%% plot the trajectories

if DEBUG,

ncontours = numel(trx);
nframes = numel(timestamps);
colors = jet(ncontours)*.7;
for t = 1:nframes,
  hold off;
  for i = 1:ncontours,
    if t > trx(i).endframe || t < trx(i).firstframe, 
      continue;
    end
    j = t + trx(i).off;
    plot(trx(i).xcontour{j},trx(i).ycontour{j},'.-','color',colors(i,:));
    hold on;
    drawellipse(trx(i).x(j),trx(i).y(j),trx(i).theta(j),trx(i).a(j)*2,trx(i).b(j)*2,'color','r');
  end
  axis equal;
  drawnow;
end

end

%% create output directory

fprintf('Saving...\n');

% make sure the root output directory exists
if ~exist(rootoutputdir,'dir'),
  mkdir(rootoutputdir);
end

% create the experiment directory
outexpdir = fullfile(rootoutputdir,experiment_name);
if ~exist(outexpdir,'dir'),
  mkdir(outexpdir);
end

%% create per-frame directory

perframedir = fullfile(outexpdir,perframedirstr);
if ~exist(perframedir,'dir'),
  mkdir(perframedir);
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
    save(outfile,'data','units');
    trx = rmfield(trx,perframefn);
  end
end

for i = 1:numel(trxfns),
  fn = trxfns{i};
  if isfield(trx,fn),
    data = {trx.(fn)}; %#ok<NASGU>
    % units not set yet
    units = parseunits('unit'); %#ok<NASGU>
    outfile = fullfile(perframedir,[fn,'.mat']);
    save(outfile,'data','units');    
  end
end
%% create the trx file

outmatname = fullfile(outexpdir,trxfilestr);
save(outmatname,'trx');
