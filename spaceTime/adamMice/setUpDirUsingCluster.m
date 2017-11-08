function [expdirlistfile,expdirs,expdirsCreated,expdirs_linux] = setUpDirUsingCluster(rootdirs,varargin)

if nargin < 1,
  rootdirs = {};
end

expdirs = {};
expdirsCreated = {};
expdirs_linux = {};

% ignore features arguments -- these will be computed on the cluster
[expdirlistfile,doforce,linuxroot,args] = myparse_nocheck(varargin,'expdirlistfile','','doforce',false,...
  'linuxroot','/groups/hantman/hantmanlab');

% file to output directories to process to
if isempty(expdirlistfile),
  expdirlistfile = sprintf('JAABADirsForCluster%s.txt',datestr(now,'yyyymmddTHHMMSS'));
end

argisremove = find(strcmpi(args(1:2:end-1),'features'));
if ~isempty(argisremove),
  args([argisremove,argisremove+1]) = [];
end
args = [args,'features',false,'doforce',doforce];

% select root directories if not given
if isempty(rootdirs),
  p = fileparts(mfilename('fullpath'));
  rcfilename = fullfile(p,'.setUpDir.mat');
  startdir = '.';
  if exist(rcfilename,'file'),
    try %#ok<TRYNC>
      rc = load(rcfilename);
      if exist(rc.startdir,'dir'),
        startdir = rc.startdir;
      end
    end
  end
  rootdirs = uipickfiles('FilterSpec',startdir,'DirsOnly',true,'Output','cell');
  if isnumeric(rootdirs) || isempty(rootdirs),
    fprintf('No root directories selected, returning\n');
    return;
  end
end

% if only one root directory is given
if ~iscell(rootdirs),
  rootdirs = {rootdirs};
end

if ~isempty(rootdirs),
  p = fileparts(mfilename('fullpath'));
  rcfilename = fullfile(p,'.setUpDir.mat');
  p = fileparts(rootdirs{end});
  if exist(p,'dir'),
    startdir = p; %#ok<NASGU>
    try %#ok<TRYNC>
      save(rcfilename,'startdir');
    end
  end
end


% run the annotation and directory creating steps
expdirsCreated = {};
for i = 1:numel(rootdirs),
  expdirs1 = setUpDir(rootdirs{i},args{:});
  if isempty(expdirs1),
    fprintf('%d: No directories created within root directory %s.\n',i,rootdirs{i});
  else
    fprintf('%d: Created the following directories within root directory %s:\n',i,rootdirs{i});
    fprintf('  %s\n',expdirs1{:});
    expdirsCreated = [expdirsCreated,expdirs1(:)']; %#ok<AGROW>
  end
end

% get a list of experiment directories to compute features for
expdirs = {};
for i = 1:numel(rootdirs),
  expdirs1 = genAllFeaturesExpDirList(rootdirs{i},'doforce',doforce,args{:});
  if isempty(expdirs1),
    fprintf('%d: No directories within root directory %s to compute features for.\n',i,rootdirs{i});
  else
    fprintf('%d: The following directories within root directory %s need processing:\n',i,rootdirs{i});
    fprintf('  %s\n',expdirs1{:});
    expdirs = [expdirs,expdirs1]; %#ok<AGROW>
  end
end

if isempty(expdirs),
  fprintf('No directories to process.\n');
  return;
end

expdirs_linux = {};
if ispc,
  for i = 1:numel(expdirs),
    % strip off the <driveletter>:\
    m = regexp(expdirs{i},'^[A-Z]:\\(.*)$','tokens','once');
    if isempty(m),
      fprintf('JAABA directory %s does not start with <driveletter>:\ -- not sure how to make this compatible with the cluster, skipping.\n',expdirs{i});
      continue;
    end
    % change file separator
    m = strrep(m{1},'\','/');
    expdirs_linux{end+1} = [linuxroot,'/',m]; %#ok<AGROW>
  end
elseif ismac,
  % file share is in the Mac path
  [~,fileshare] = fileparts(linuxroot);
  for i = 1:numel(expdirs),
    % strip off the /Volumes part
    m = regexp(expdirs{i},'^/Volumes/([^/]*)/(.*)','tokens','once');
    if isempty(m),
      fprintf('JAABA directory %s could not be parsed as /Volumes/(something)/(rest) -- not sure how to make this compatible with the cluster, skipping.\n',expdirs{i});
      continue;
    end
    if ~strcmp(fileshare,m{1}),
      fprintf('JAABA directory %s does not match specified file share %s -- not sure how to make this compatible with the cluster, skipping.\n',expdirs{i},fileshare);
      continue;
    end
    expdirs_linux{end+1} = fullfile(linuxroot,m{2}); %#ok<AGROW>
  end
else
  expdirs_linux = cell(size(expdirs));
  for i = 1:numel(expdirs),
    expdirs_linux{i} = absolutifyFileName(expdirs{i},pwd);
  end
end

fprintf('Outputting list of JAABA directories to process on the cluster to file %s...\n',expdirlistfile);
fid = fopen(expdirlistfile,'w');
if fid < 0,
  error('Could not write to file %s',expdirlistfile);
end

fprintf(fid,'%s\n',expdirs_linux{:});
fclose(fid);

fprintf('\nOutput %d JAABA experiments to file %s to be processed on the cluster.\n\n',numel(expdirs_linux),expdirlistfile);


if isunix && ~ismac,
  expdirlistfile_abs = absolutifyFileName(expdirlistfile,pwd);
  fprintf('To process on the cluster, do the following:\n');
  fprintf('1. Log into login2 (this is done with the terminal command "ssh login2")\n\n');
  fprintf('2. On login2, run the command:\n');
  fprintf('    /groups/branson/bransonlab/share/qsub_genAllFeatures.pl %s\n\n',expdirlistfile_abs);
  fprintf('To check on the progress of your jobs, you can type "qstat" at the command prompt to get a list of jobs running.\n');
else
  fprintf('To process on the cluster, do the following:\n');
  fprintf('1. Put the file created by this function, %s, on the Janelia file system, e.g. in the hantmanlab file share.\n',expdirlistfile);
  fprintf('2. Log into login2 (this is done with the terminal command "ssh login2")\n\n');
  fprintf('3. On login2, run the command:\n');
  fprintf('    /groups/branson/bransonlab/share/qsub_genAllFeatures.pl <JAABADirectoryListFile>\n');
  [~,n,e] = fileparts(expdirlistfile);
  n = [n,e];
  fprintf('where <JAABADirectoryListFile> is the *full Linux path* to the location of the file %s on the Janelia file server,\n',n);
  fprintf('Example: /groups/branson/bransonlab/share/qsub_GenAllFeatures.pl /groups/hantman/hantmanlab/JAABADirLists/%s\n\n',n);
  fprintf('To check on the progress of your jobs, you can type "qstat" at the command prompt to get a list of jobs running.\n');
end