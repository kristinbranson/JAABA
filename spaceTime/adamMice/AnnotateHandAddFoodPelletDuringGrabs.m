function [annotations,expdirs,annotationi,tmpfile] = AnnotateHandAddFoodPelletDuringGrabs(varargin)

annotations = [];
annotationi = [];

persistent startdir;
if isempty(startdir),
  startdir = '.';
end

[expdirs,rootdatadirs,frontside,maxboutlengthignore,restartfile,tmpfile,windowradius] = myparse(varargin,...
  'expdirs',{},'rootdatadirs',{},...
  'frontside',true,...
  'maxboutlengthignore',2,...
  'restartfile','',...
  'tmpfile',sprintf('AutoSaveHandAndPelletAnnotations%s.mat',datestr(now,'yyyymmddTHHMMSS')),...
  'windowradius',5);

if ~isempty(restartfile),
  
  fprintf('Loading previous state from %s, not doing any checking!\n',restartfile);
  load(restartfile);
else
  

%% get experiment directories
if isempty(expdirs),
  if isempty(rootdatadirs),    
    rootdatadirs = uipickfiles('FilterSpec',startdir,'DirsOnly',true,'Output','cell');
  end
  if isnumeric(rootdatadirs) || isempty(rootdatadirs),
    fprintf('No root directories selected, returning\n');
    return;
  end
  startdir = fileparts(rootdatadirs{end});

  expdirs = {};
  for i = 1:numel(rootdatadirs),
    expdirscurr = findAllSubDirs(rootdatadirs{i});
    fprintf('%s: %d experiment directories found.\n',rootdatadirs{i},numel(expdirscurr));
    expdirs = [expdirs,expdirscurr]; %#ok<AGROW>
  end
  if isempty(expdirs),
    fprintf('No experiment directories found within any of the root directories, returning\n');
    return;
  end
end
nexps = numel(expdirs);

%% end part skipped if we are restarting

end

%% find grab scores files, filter out experiments without scores
foundscores = false(1,nexps);
scoresfiles = cell(1,nexps);
for i = 1:nexps,
  
  expdir = expdirs{i};
  tmp = dir(fullfile(expdir,'scores_Grab*.mat'));
  if isempty(tmp),
    continue;
  end
  maxtimestamp = 0;
  nfound = 0;
  for j = 1:numel(tmp),
    if ~isempty(regexp(tmp(j).name,'bak','once')),
      continue;
    end
    nfound = nfound + 1;
    if tmp(j).datenum > maxtimestamp,
      bestj = j;
      maxtimestamp = tmp(j).datenum;
    end
  end
  if nfound == 0,
    continue;
  end
  scoresfiles{i} = fullfile(expdir,tmp(bestj).name);
  if nfound > 1,
    fprintf('Found multiple grab scores files in %s, using %s\n',expdir,tmp(bestj).name);
  end
  foundscores(i) = true;
end

if any(~foundscores),
  fprintf('No grab scores files found for the following %d experiments, ignoring them.\n',nnz(~foundscores));
  fprintf('%s\n',expdirs{~foundscores});
end

expdirs = expdirs(foundscores);
scoresfiles = scoresfiles(foundscores);
nexps = nnz(foundscores);

%% find the first frame of grab

firstgrab = nan(1,nexps);
for i = 1:nexps,
  sd = load(scoresfiles{i});
  scorescurr = sd.allScores.postprocessed{1};
  [t0s,t1s] = get_interval_ends(scorescurr);
  idxlongenough = find(t1s-t0s > maxboutlengthignore);
  if isempty(idxlongenough),
    fprintf('No grabs found for %s\n',expdirs{i});
    continue;
  end
  firstgrab(i) = t0s(idxlongenough(1));
  fprintf('First grab for %s at frame %d in bout of length %d\n',expdirs{i},firstgrab(i),t1s(idxlongenough(1))-t0s(idxlongenough(1)));
end

%% read in these images

fprintf('Reading in frames for annotating...\n');
ims = cell(1,nexps);
nc = ceil(sqrt(nexps));
nr = ceil(nexps/nc);
hfig = figure;
hax = createsubplots(nr,nc,.01);
iplocs = [];
for i = 1:nexps,
  if isnan(firstgrab(i)),
    continue;
  end
  if frontside,
    moviefilestr = 'movie_comb.avi';
  else
    moviefilestr = 'movie.avi';
  end
  moviefile = fullfile(expdirs{i},moviefilestr);
  if ~exist(moviefile,'file'),
    fprintf('Could not find movie file %s, skipping\n',moviefile);
    continue;
  end
  [readframe,nframes] = get_readframe_fcn(moviefile);
  if nframes < firstgrab(i),
    warning('%s: first grab frame is at %d > movie length = %d, skipping\n',expdirs{i},firstgrab(i),nframes);
    continue;
  end
  fprintf('%d/%d: reading in first grab frame for %s\n',i,nexps,expdirs{i});
  for j = firstgrab(i)-windowradius:firstgrab(i)+windowradius,
    if j < 1 || j > nframes,
      continue;
    end
    ims{i}(:,:,j-firstgrab(i)+windowradius+1) = rgb2gray(readframe(j));
  end
  
  trxfile = fullfile(expdirs{i},'trx.mat');
  td = load(trxfile);
  iplocs = structarrayset(iplocs,i,td.trx(1).arena);
  
  imagesc(ims{i}(:,:,windowradius+1),'Parent',hax(i),[0,255]);
  colormap(hax(i),'gray');
  axis(hax(i),'image','off');
  text(1,1,expdirs{i}(max(1,numel(expdirs{i})-20):end),'Parent',hax(i),'Color','r','Interpreter','none');
  drawnow;
end

%% annotate all experiments

fprintf('Temporary results saved to file %s, restart from here if something crashes!\n',tmpfile);
save(tmpfile,'expdirs','firstgrab','foundscores','frontside','maxboutlengthignore','nexps','rootdatadirs','scoresfiles','startdir','windowradius');

if isempty(annotationi),
  annotationi = 0;
end

if ~exist('annotations','var'),
  annotations = [];
end
ipfns0 = setdiff(fieldnames(iplocs),{'face'});

figpos = [];

ipnames = {'PawSide','PelletSide','PawFront','PelletFront'};
for annotationi = annotationi+1:nexps,
  if isempty(ims{annotationi}),
    continue;
  end
  while true,
    fprintf('%d/%d: Annotate locations for %s\n',annotationi,nexps,expdirs{annotationi});
    
    annObj = HandleObj;
    waitfor(playfmf('moviename',{@(t) deal(ims{annotationi}(:,:,t),t/200),2*windowradius+1},...
      'startframe',windowradius+1,'ipnames',ipnames,...
      'doaskleftright',false,'annDataObj',annObj,'figpos',figpos));
    drawnow;
    if ~isempty(annObj.data) && isstruct(annObj.data) && isfield(annObj.data,'figpos'),
      figpos = annObj.data.figpos;
      annObj.data = rmfield(annObj.data,'figpos');
    end
    
    if isempty(annObj.data),
      res = questdlg('No interest points annotated. Redo, skip or quit?','No interest points annotated','Redo','Skip','Quit','Redo');
      if ~strcmpi(res,'Redo'),
        break;
      end
    elseif ~all(isfield(annObj.data,ipnames)),
      missingips = setdiff(ipnames,fieldnames(annObj.data));
      missingips = sprintf('%s ',missingips{:});
      res = questdlg(sprintf('The following interest points were not annotated: %s. Redo, skip or quit?',missingips),...
        'Interest points missing','Redo','Skip','Quit','Redo');
      if ~strcmpi(res,'Redo'),
        break;
      end
    else
      res = 'OK';
      break;
    end
  end
  
  if strcmpi(res,'Quit'),
    return;
  end

  iplocscurr = struct;
  for j = 1:numel(ipfns0),
    iplocscurr.(['ip_',ipfns0{j}]) = iplocs(annotationi).(ipfns0{j});
  end
  if strcmpi(res,'Skip'),
    annotationscurr = iplocscurr;
    for j = 1:numel(ipnames),
      annotationscurr.(ipnames{j}) = nan;
    end
    annotations = structarrayset(annotations,annotationi,annotationscurr);
  else
    annotations = structarrayset(annotations,annotationi,structmerge(annObj.data,iplocscurr));
  end

  save('-append',tmpfile,'annotationi','annotations');
end