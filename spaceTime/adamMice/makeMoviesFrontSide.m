function expdirs = makeMoviesFrontSide(rootdir,varargin)

[doforce,dosoftlink,rotatefrontview,rotatesideview] = myparse(varargin,...
  'doforce',false,...
  'dosoftlink',false,...
  'rotatefrontview',0,'rotatesideview',0);

fcn = @(x)lclCreateJAABADirIfFrontSideMoviesPresent(x,doforce,dosoftlink,...
  rotatefrontview,rotatesideview);
m = depthFirstSearch(rootdir,fcn);
expdirs = m.values;
expdirs = cat(1,expdirs{:});

function out = lclCreateJAABADirIfFrontSideMoviesPresent(edir,doforce,dosoftlink,rotatefrontview,rotatesideview)

out = cell(0,1); % jaaba dirs created

FRONT_MOVIE_NAME = 'movie_frt.avi';
SIDE_MOVIE_NAME = 'movie_sde.avi';
TRX_NAME = 'trx.mat';

dd = dir(fullfile(edir,'*_front_*avi'));
for ndx = 1:numel(dd)
  if dd(ndx).isdir,
    continue;
  end
  frontnameext = dd(ndx).name;
  [~,name,ext] = fileparts(frontnameext);
  fullfrontname = fullfile(edir,frontnameext);
  
  % guid does not need to match
  m = regexp(name,'^(?<pre>.*)_front_(?<mid>.*)(?<guid>_guid_[a-f0-9]+)(?<post>_.*)?$','names','once');
  if ~isempty(m),
    shortname = [m.pre,'_front_',m.mid,m.post];
    re = fullfile(edir,[m.pre,'_side_',m.mid,'_guid_*',m.post,ext]);
    tmp = dir(re);
    if ~isempty(tmp),
      m2 = regexp({tmp.name},['^',regexptranslate('escape',m.pre),'_side_',...
        regexptranslate('escape',m.mid),'(?<guid>_guid_[a-f0-9]+)',...
        regexptranslate('escape',m.post),...
        regexptranslate('escape',ext),'$'],'once');
      ismatch = ~cellfun(@isempty,m2);
      tmp = tmp(ismatch);
    end

    if ~isempty(tmp),
      if numel(tmp) > 1,
        warning('Multiple videos match %s, choosing first one: %s',re,tmp(1).name);
      end
      sidenameext = tmp(1).name;
    else
      warning('Could not find side video matching %s, skipping %s',re,frontnameext);
      continue;
    end
  else
    shortname = name;
    sidenameext = regexprep(frontnameext,'_front_','_side_');
  end
  fullsidename = fullfile(edir,sidenameext);


  % create JAABAdir
  jdir = fullfile(edir,regexprep(shortname,'_front_','_'));
  if ~exist(jdir,'dir')
    mkdir(jdir);
  end
      
  lclCopyOrLinkMov(jdir,fullfrontname,fullfile(jdir,FRONT_MOVIE_NAME),doforce,dosoftlink,rotatefrontview);
  lclCopyOrLinkMov(jdir,fullsidename,fullfile(jdir,SIDE_MOVIE_NAME),doforce,dosoftlink,rotatesideview);
  
  lclCreateTrx(jdir,TRX_NAME,{FRONT_MOVIE_NAME,SIDE_MOVIE_NAME},doforce);
  
  lclCreatePFdir(jdir,TRX_NAME,doforce);
      
  out{end+1,1} = jdir; %#ok<AGROW>
end

function lclCopyOrLinkMov(jaabadir,src,dest,doforce,dosoftlink,rotateview)
    
if doforce
  if exist(dest,'file')
    delete(dest);
  end
end
    
if ~exist(dest,'file'),
  try
    if rotateview ~= 0,
      rotatemov(src,dest,rotateview);
    else
      if isunix && ~ismac && dosoftlink
        fprintf('Soft-linking %s to %s...\n',src,dest);
        [~,srcn,srce] = fileparts(src);
        % use relative links so that maybe this works when mounting from windows/mac?
        srcn = [srcn,srce];
        [~,destn,deste] = fileparts(dest);
        destn = [destn,deste];
        syscmd = sprintf('ln -s ../%s %s',srcn,fullfile(jaabadir,destn));
        unix(syscmd);
      else
        fprintf('Copying %s to %s...\n',src,dest);
        copyfile(src,dest);
      end
    end
  catch ME,
    fprintf('Could not copy/link video file %s (%s)\n',src,ME.message);
  end
else
  fprintf('%s exists, not re-copying it.\n',dest);
end

function lclCreateTrx(jdir,trxname,movienames,doforce)
  
trxnamefull = fullfile(jdir,trxname);
if ~iscell(movienames),
  movienames = {movienames};
end
movienamesfull = cell(size(movienames));
for i = 1:numel(movienames),
  movienamesfull{i} = fullfile(jdir,movienames{i});
end

if doforce && exist(trxnamefull,'file')
  delete(trxnamefull);
end
  
if ~exist(trxnamefull,'file')
  try
    fprintf('Generating trx file %s...\n',trxnamefull);
    nframes = inf;
    for i = 1:numel(movienamesfull),
      [~,nframescurr,~,~] = get_readframe_fcn(movienamesfull{i});
      nframes = min(nframes,nframescurr);
    end
    if isinf(nframes),
      error('Could not read nframes for %s',rootdir);
    end
    T.trx = genTrack(nframes,[640 480]);
    save(trxnamefull,'-struct','T');
  catch ME,
    fprintf('Could not create trx file for %s (%s)\n',jdir,ME.message);
  end
else
  fprintf('Trx file exists for %s, not regenerating\n',jdir);
end

function lclCreatePFdir(jdir,trxname,doforce)
pfnamefull = fullfile(jdir,'perframe');
trxnamefull = fullfile(jdir,trxname);

if doforce && exist(pfnamefull,'dir')
  delete(fullfile(pfnamefull,'*'));
  [success,msg] = rmdir(pfnamefull);
  if ~success
    error('Failed to remove existing per-frame directory %s (%s)',pfnamefull,msg);
  end
end
  
if ~exist(pfnamefull,'dir')
  [success,msg,~] = mkdir(jdir,'perframe');
  if ~success,
    error('Error creating per-frame directory %s (%s)',pfnamefull,msg);
  end
  
  try
    fprintf('Generating perframe dir %s...\n',pfnamefull);
    t = load(trxnamefull);
    tmp = struct();
    tmp.data = {t.trx.a_mm};
    tmp.units = parseunits('mm');
    save(fullfile(pfnamefull,'a_mm.mat'),'-struct','tmp');
  catch ME
    fprintf('Could not create perframe dir (%s)\n',ME.message);
  end
else
  fprintf('Perframe dir %s exists, not regenerating\n',pfnamefull);
end


  

  