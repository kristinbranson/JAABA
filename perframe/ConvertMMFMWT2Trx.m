function [outexpdir] = ConvertMMFMWT2Trx(inexpdir,rootoutputdir,varargin)

%inexpdir = '/groups/branson/bransonlab/projects/larvae_species/data/rawdata/Amovieshort';
% rootoutputdir = 'C:\Data\larvae_species\data';
[rootdatadir,experiment_name] = fileparts(inexpdir);

%% set parameters

% default parameters
contourext = 'outline';
spineext = 'spine';
blobext = 'blobs';
movieext = 'mmf';
trxfilestr = 'trx.mat';
perframedirstr = 'perframe';
moviefilestr = 'movie.mmf';
DEBUG = false;
makesoftlink = true;
% for now pxpermm is a static input
pxpermm = 1;
dotransposeimage = false;

[contourext,spineext,blobext,movieext,...
  analysisname,...
  trxfilestr,perframedirstr,moviefilestr,...
  pxpermm,dotransposeimage,...
  DEBUG,makesoftlink] = ...
  myparse(varargin,...
  'contourext',contourext,...
  'spineext',spineext,...
  'blobext',blobext,...
  'movieext',movieext,...
  'analysisdir','',...
  'trxfilestr',trxfilestr,...
  'perframedirstr',perframedirstr,...
  'moviefilestr',moviefilestr,...
  'pxpermm',pxpermm,...
  'dotransposeimage',dotransposeimage,...
  'debug',DEBUG,...
  'makesoftlink',makesoftlink);

% if no analysis directory explicitly given, choose the latest
if isempty(analysisname),
  files = dir(inexpdir);
  isanalysisdir = [files.isdir] & ~cellfun(@isempty,regexp({files.name},'^\d{8}_\d{6}$','once'));
  analysisdirs = {files(isanalysisdir).name};
  analysisdirs = sort(analysisdirs);
  analysisname = analysisdirs{end};
end

%% names of files

% find contour files for each trajectory
analysisdir = fullfile(inexpdir,analysisname);
if ~exist(analysisdir,'dir'),
  error('Analysis directory %s does not exist',analysisdir);
end

files = dir(fullfile(analysisdir,[experiment_name,'.*.',contourext]));
contournames = {};
ids = [];
for i = 1:numel(files),
  match = regexp(files(i).name,['\.(\d+)\.',contourext,'$'],'tokens','once');
  if isempty(match),
    continue;
  end
  contournames{end+1} = fullfile(analysisdir,files(i).name); %#ok<AGROW>
  ids(end+1) = str2double(match{1}); %#ok<AGROW>
end

if isempty(contournames),
  error('Could not find any contours matching the pattern %s.*.%s',fullfile(rootdatadir,experiment_name),contourext);
end

% find corresponding spine files
isbaddata = false(1,numel(contournames));
spinenames = cell(1,numel(contournames));
for i = 1:numel(ids),
 spinenames{i} = [contournames{i}(1:end-numel(contourext)),spineext];
 if ~exist(spinenames{i},'file'),
   warning('Spine file %s does not exist for contour file %s, id %d',spinenames{i},contournames{i},ids(i));
   isbaddata(i) = true;
 end
end

spinenames(isbaddata) = [];
contournames(isbaddata) = [];
ids(isbaddata) = [];

% find blobs file
nblobsfiles = ceil(max(ids)/1000);
blobnames = cell(1,nblobsfiles);
for i = 1:nblobsfiles,
  blobnames{i} = fullfile(analysisdir,sprintf('%s_%05dk.%s',experiment_name,i-1,blobext));
  if ~exist(blobnames{i},'file'),
    warning('Blob file %s does not exist',blobnames{i});
  end
end

ncontours = numel(contournames);

%% read blobs files

blobs = [];
for i = 1:numel(blobnames),
  blobs = structappend(blobs,ReadMWTBlobsFile(blobnames{i},'dotransposeimage',dotransposeimage));
end

%% read all contour files

contours = ReadMWTContours(contournames,ids,'dotransposeimage',dotransposeimage);

%% read all spine files

spines = ReadMWTSpines(spinenames,ids,'dotransposeimage',dotransposeimage);

%% merge together

blobfnsmerge = {'x','y','a','b','theta','area','width','length'};
blobfnscopy = {};
contourfnsmerge = {'xcontour','ycontour'};
contourfnscopy = {};
spinefnsmerge = {'xspine','yspine'};
spinefnscopy = {};
[trx,timestamps] = MergeMWTData(blobs,blobfnsmerge,blobfnscopy,...
  contours,contourfnsmerge,contourfnscopy,...
  spines,spinefnsmerge,spinefnscopy);

%% convert to mm
fnsmm = {'x','y','a','b','theta','area','width','length'};
for i = 1:numel(trx),
  for j = 1:numel(fnsmm),
    fnin = fnsmm{j};
    fnout = [fnin,'_mm'];
    trx(i).(fnout) = trx(i).(fnin)/pxpermm;
  end
end

%% add stuff

arena = struct('x',nan,'y',nan,'r',nan);
dts = diff(timestamps);
fps = 1/nanmean(dts);
for i = 1:numel(trx),
  trx(i).arena = arena;
  trx(i).fps = fps;
  trx(i).pxpermm = pxpermm;
end

%% plot the trajectories

if DEBUG,

nframes = numel(timestamps);
colors = jet(ncontours)*.7;
for t = 1:nframes,
  hold off;
  for i = 1:ncontours,
    if t > trx(i).endframe || t < trx(i).firstframe, 
      continue;
    end
    j = t + trx(i).off;
    plot([trx(i).xcontour{j},trx(i).xcontour{j}(1)],[trx(i).ycontour{j},trx(i).ycontour{j}(1)],'-','color',colors(i,:));
    hold on;
    plot(trx(i).xspine(:,j),trx(i).yspine(:,j),'.-','color',colors(i,:));
    drawellipse(trx(i).x(j),trx(i).y(j),trx(i).theta(j),trx(i).a(j)*2,trx(i).b(j)*2,'color','r');
  end
  axis equal;
  drawnow;
end

end

%% create output directory

% make sure the root output directory exists
if ~exist(rootoutputdir,'dir'),
  mkdir(rootoutputdir);
end

% create the experiment directory
outexpdir = fullfile(rootoutputdir,experiment_name);
if ~exist(outexpdir,'dir'),
  mkdir(outexpdir);
end

%% create the trx file

outmatname = fullfile(outexpdir,trxfilestr);
save(outmatname,'trx');

%% copy over the movie
inmoviename = fullfile(inexpdir,[experiment_name,'.',movieext]);
outmoviename = fullfile(outexpdir,moviefilestr);

if ~exist(outmoviename,'file'),
  if ~exist(inmoviename,'file'),
    error('Movie file %s does not exist',inmoviename);
  end
  if isunix && makesoftlink,
    cmd = sprintf('ln -s %s %s',inmoviename,outmoviename);
    unix(cmd);
  else
    [success,msg] = copyfile(inmoviename,outmoviename);
    if ~success,
      error('Error copying file %s to %s: %s',inmoviename,outmoviename,msg);
    end
  end
end

%% save the per-frame features
perframedir = fullfile(outexpdir,perframedirstr);
if ~exist(perframedir,'dir'),
  mkdir(perframedir);
end

% to do: save choreography's per-frame features
