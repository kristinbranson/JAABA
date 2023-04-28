function create_dummy_trx_file(expdir,x,y,names,varargin)

[movfilestr,trxfilestr,a,b,theta,fps,forcecompute,face] = myparse(varargin,...
  'movfilestr',[],'trxfilestr','trx.mat',...
  'a',5,'b',2,'theta',0,'fps',400,'forcecompute',false,'face','right');

defaultmovfilestrs = {
  'movie_comb.avi'
  {'movie_sde.avi','movie_frt.avi'}
  };


if isempty(movfilestr),
  movfilestrs = defaultmovfilestrs;
else
  movfilestrs = {movfilestr};
end
for i = 1:numel(movfilestrs),
  movfilestr = movfilestrs{i};
  nframes = get_nframes(expdir,movfilestr);
  if ~isinf(nframes),
    break;
  end
end

if isinf(nframes),
  error('Coult not find movie file in directory %s',expdir);
end

trxfile = fullfile(expdir,trxfilestr);
if exist(trxfile,'file')
  if forcecompute,
    fprintf('trxfile %s exists, overwriting...\n',trxfile);
  else
    fprintf('trxfile %s already exists, not computing.\n',trxfile);
    return;
  end
end

arena = struct;
for i = 1:numel(x),
  arena.(names{i}) = [x(i),y(i)];
end
arena.face = face;

npts = numel(x);
for i = 1:npts,

  trxcurr = struct;
  trxcurr.x = x(i) + zeros(1,nframes);
  trxcurr.y = y(i) + zeros(1,nframes);
  trxcurr.theta = theta + zeros(1,nframes);
  trxcurr.a = a + zeros(1,nframes);
  trxcurr.b = b + zeros(1,nframes);
  trxcurr.x_mm = x(i) + zeros(1,nframes);
  trxcurr.y_mm = y(i) + zeros(1,nframes);
  trxcurr.theta_mm = theta + zeros(1,nframes);
  trxcurr.a_mm = a + zeros(1,nframes);
  trxcurr.b_mm = b + zeros(1,nframes);
  trxcurr.id = i-1;
  trxcurr.moviename = '';
  trxcurr.annname = '';
  trxcurr.firstframe = 1;
  trxcurr.nframes = nframes;
  trxcurr.endframe = nframes;
  trxcurr.timestamps = (1:nframes)/fps;
  trxcurr.fps = fps;
  trxcurr.pxpermm = 1;
  trxcurr.README = 'arena is raw ips; .x, .y are massaged for JAABA-playback';
  trxcurr.arena = arena;
  trxcurr.concatmov = [];

  if i == 1,
    trx = repmat(trxcurr,[1,npts]);
  else
    trx(i) = trxcurr;
  end

end

save(trxfile,'trx');

function nframes = get_nframes(expdir,movfilestrs)

nframes = inf;
if ischar(movfilestrs),
  movfilestrs = {movfilestrs};
end

for i = 1:numel(movfilestrs),
  movfilestr = movfilestrs{i};
  if ~iscell(movfilestr),
    movfilestr = {movfilestr};
  end
  nframes = inf;
  for j = 1:numel(movfilestr),
    movfile = fullfile(expdir,movfilestr{j});
    if exist(movfile,'file'),
      [~,nframescurr,fid,~] = get_readframe_fcn(movfile);
      if fid > 0,
        fclose(fid);
      end
      nframes = min(nframes,nframescurr);
    end
  end
end
