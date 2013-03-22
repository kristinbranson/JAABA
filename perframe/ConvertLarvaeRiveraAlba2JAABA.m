function [success,msg] = ConvertLarvaeRiveraAlba2JAABA(varargin)

success = false;
msg = '';

[inmoviefile,intrxfile,...
  expdir,moviefilestr,trxfilestr,perframedirstr,...
  dosoftlink] = myparse(varargin,...
  'inmoviefile','','intrxfile','',...
  'expdir','','moviefilestr','movie.ufmf','trxfilestr','trx.mat','perframedirstr','perframe',...
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

% output file locations
moviefile = fullfile(expdir,moviefilestr);
trxfile = fullfile(expdir,trxfilestr);
perframedir = fullfile(expdir,perframedirstr);

% load in the trx
try
  td = load(intrxfile,'-regexp','^(trx|timestamps)$');
  if ~isfield(td,'trx'),
    msg = sprintf('No trx variable in file %s',intrxfile);
    return;
  end
  trx = td.trx;
  if isfield(td,'timestamps'),
    timestamps = td.timestamps;
  else
    timestamps = [];
  end
  clear td;
catch ME,
  msg = sprintf('Could not load trx from file %s: %s',intrxfile,getReport(ME));
  return;
end

% compute timestamps
if isempty(timestamps),
  nframes = max([trx.endframe]);
  dt = zeros(1,nframes-1);
  n = zeros(1,nframes-1);
  for i = 1:numel(trx),
    dt(trx(i).firstframe:trx(i).endframe-1) = ...
      dt(trx(i).firstframe:trx(i).endframe-1) + trx(i).dt;
    n(trx(i).firstframe:trx(i).endframe-1) = ...
      n(trx(i).firstframe:trx(i).endframe-1) + 1;
  end
  dt(n>0) = dt(n>0) / n(n>0);
  mediandt = median(dt(dt>0));
  dt(n==0) = mediandt;
  timestamps = [0,cumsum(dt)]; %#ok<NASGU>
end

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
if strcmp(fullfile(inmoviefile),fullfile(moviefile)),
  fprintf('Input and out movie files are the same, not copying/linking.\n');
else
  
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