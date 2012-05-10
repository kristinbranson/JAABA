function ConvertMouseTrx2Trx(rootinputdir,rootoutputdir,inexpname,inseqfilestr,pxpermm,sex,varargin)

% names of files within each experiment
outtrxfilestr = 'trx.mat';
outseqfilestr = 'movie.seq';
outseqindexfilestr = 'movie.mat';
makelinks = true;

[outtrxfilestr,outseqfilestr,outseqindexfilestr,makelinks] = ...
  myparse(varargin,'outtrxfilestr',outtrxfilestr,...
  'outseqfilestr',outseqfilestr,...
  'outseqindexfilestr',outseqindexfilestr,...
  'makelinks',makelinks);

%% get paths to input files

% name of input seq file
inseqfile = fullfile(rootinputdir,inseqfilestr);
if ~exist(inseqfile,'file'),
  error('Seq file %s does not exist',inseqfile);
end

% name of input index file, assumed to be <inexpname>.mat
inseqindexfile = fullfile(rootinputdir,[inexpname,'.mat']);
if ~exist(inseqindexfile,'file'),
  error('Seq index file %s does not exist',inseqindexfile);
end

% name of input trx file, assumed to be <inexpname>t.mat
intrxfile = fullfile(rootinputdir,[inexpname,'t.mat']);
if ~exist(intrxfile,'file'),
  error('Trajectory file %s does not exist',intrxfile);
end

%% get paths to output files

% remove bad characters from names
outexpname = regexprep(inexpname,'[^\w.]+','_');
if ~strcmp(inexpname,outexpname),
  fprintf('Converting %s to %s\n',inexpname,outexpname);
end

% create the output experiment directory
outexpdir = fullfile(rootoutputdir,outexpname);
if ~exist(outexpdir,'dir'),
  mkdir(outexpdir);
end

outtrxfile = fullfile(outexpdir,outtrxfilestr);
outseqfile = fullfile(outexpdir,outseqfilestr);
outseqindexfile = fullfile(outexpdir,outseqindexfilestr);

%% link to/copy files

if makelinks,
  if ispc,
    cmd = sprintf('mkshortcut.vbs /target:"%s" /shortcut:"%s"',inseqfile,outseqfile);
    fprintf('Making a Windows shortcut file at "%s" with target "%s"\n',outseqfile,inseqfile);
  else
    cmd = sprintf('ln -s "%s" "%s"',inseqfile,outseqfile);
    fprintf('Soft-linking from "%s" with target "%s"\n',outseqfile,inseqfile);
  end
  system(cmd);

  if ispc,
    cmd = sprintf('mkshortcut.vbs /target:"%s" /shortcut:"%s"',inseqindexfile,outseqindexfile);
    fprintf('Making a Windows shortcut file at "%s" with target "%s"\n',outseqindexfile,inseqindexfile);
  else
    cmd = sprintf('ln -s "%s" "%s"',inseqindexfile,outseqindexfile);
    fprintf('Soft-linking from "%s" with target "%s"\n',outseqindexfile,inseqindexfile);
  end
  system(cmd);
  
else
  fprintf('Copying "%s" to "%s"\n',inseqfile,outseqfile);
  copyfile(inseqfile,outseqfile);

  fprintf('Copying "%s" to "%s"\n',inseqindexfile,outseqindexfile);
  copyfile(inseqindexfile,outseqindexfile);
end

%%

% trx = ConvertMouseTrx2Trx(intrxname,moviename,outtrxname,pxpermm,sex)

% load data
headerinfo = r_readseqinfo(inseqfile);
tmp = load(intrxfile);

% count
nmice = numel(tmp.astrctTrackers);
nframes = numel(tmp.astrctTrackers(1).m_afX);

% create new structure
trx = struct('x',{tmp.astrctTrackers.m_afX},...
  'y',{tmp.astrctTrackers.m_afY},...
  'theta',cellfun(@(x) -x, {tmp.astrctTrackers.m_afTheta},'UniformOutput',false),...
  'a',cellfun(@(x) x/2,{tmp.astrctTrackers.m_afA},'UniformOutput',false),...
  'b',cellfun(@(x) x/2,{tmp.astrctTrackers.m_afB},'UniformOutput',false),...
  'firstframe',num2cell(ones(1,nmice)),...
  'arena',cell(1,nmice),...
  'off',num2cell(zeros(1,nmice)),...
  'nframes',num2cell(nframes(ones(1,nmice))),...
  'endframe',num2cell(nframes(ones(1,nmice))),...
  'timestamps',repmat({headerinfo.m_afTimestamp-headerinfo.m_afTimestamp(1)},[1,nmice]),...
  'moviename',repmat({moviename},[1,nmice]),...
  'annname',repmat({intrxname},[1,nmice]),...
  'matname',repmat({outtrxname},[1,nmice]),...
  'x_mm',cellfun(@(x) x/pxpermm,{tmp.astrctTrackers.m_afX},'UniformOutput',false),...
  'y_mm',cellfun(@(x) x/pxpermm,{tmp.astrctTrackers.m_afY},'UniformOutput',false),...
  'a_mm',cellfun(@(x) x/pxpermm/2,{tmp.astrctTrackers.m_afA},'UniformOutput',false),...
  'b_mm',cellfun(@(x) x/pxpermm/2,{tmp.astrctTrackers.m_afB},'UniformOutput',false),...
  'theta_mm',cellfun(@(x) -x, {tmp.astrctTrackers.m_afTheta},'UniformOutput',false),...
  'dt',repmat({diff(headerinfo.m_afTimestamp)},[1,nmice]),...
  'fps',repmat({1/median(diff(headerinfo.m_afTimestamp))},[1,nmice]),...
  'pxpermm',repmat({pxpermm},[1,nmice]),...
  'sex',sex); %#ok<NASGU>

% save to file
save(outtrxfile,'trx');
