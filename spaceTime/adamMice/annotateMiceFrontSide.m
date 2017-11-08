function success = annotateMiceFrontSide(expdir,varargin)

success = false; %#ok<NASGU>

[ips,doforce,imwidth] = myparse(varargin,...
  'ips',[],...
  'doforce',false,...
  'imwidth',640);

assert(~isempty(ips),'Interest pts must be specified.');

%%% Update trx.mat 

% IMPORTANT: trx.arena contains the raw interest points as created by user,
% applicable to raw front/side video. trx.x, trx.y are for JAABA-playback
% and assume i) concatenation of front/side videos and ii) fliplr of all
% videos for left-facing mice.

TRX_FILE = 'trx.mat';

T = load(fullfile(expdir,TRX_FILE));

if ~doforce && isfield(T,'trx') && isfield(T.trx(1),'arena') && ...
  all(isfield(T.trx(1).arena,{'food' 'mouth','perch','face','mouthfront','foodfront'})),
    fprintf('Skipping trx-file update for %s, annotations exist.\n',expdir);
else  
  fprintf('Updating trx structure for %s.\n',expdir);
  
  arena = ips;  
  T.trx(1).README = 'arena is raw ips; .x, .y are massaged for JAABA-playback';
  T.trx(1).arena = arena;
  T.trx(2) = T.trx(1);
  T.trx(3) = T.trx(1);
  T.trx(4) = T.trx(1);
  T.trx(5) = T.trx(1);
  
  % trx.x, trx.y contain coords for JAABA-playback. 

  xysde = [arena.food;arena.mouth;arena.perch];
  xyfrt = [arena.foodfront;arena.mouthfront];

  % first deal with fliplr
  fliplr = strcmpi(arena.face,'left'); % hardcoded, needs to match genFeatures.m
  if fliplr
    xysde(:,1) = imwidth - xysde(:,1);
    xyfrt(:,1) = imwidth - xyfrt(:,1);
  end
  
  % now deal with concatenation. playback movie is [sidemov frontmov].
  xyfrt(:,1) = xyfrt(:,1) + imwidth;
  
  % put it in trx  
  T.trx(1).x(:) = xysde(1,1);
  T.trx(1).y(:) = xysde(1,2);
  T.trx(2).x(:) = xysde(2,1);
  T.trx(2).y(:) = xysde(2,2);
  T.trx(3).x(:) = xysde(3,1);
  T.trx(3).y(:) = xysde(3,2);
  T.trx(4).x(:) = xyfrt(1,1);
  T.trx(4).y(:) = xyfrt(1,2);
  T.trx(5).x(:) = xyfrt(2,1);
  T.trx(5).y(:) = xyfrt(2,2);
  
  T.trx(1).concatmov.fliplr = fliplr;
  T.trx(1).concatmov.imwidth = imwidth;  

  save(fullfile(expdir,TRX_FILE),'-struct','T');
end

% moving this to gen features
% %%% Create concatenated movie for JAABA playback
% 
% COMBND_MOVIE_NAME = 'movie_comb.avi';
% 
% fprintf(1,'Concatenating front/side movies for playback...\n');
% catmovname = fullfile(expdir,COMBND_MOVIE_NAME);
% if exist(catmovname,'file')==2 && doforce
%   delete(catmovname);
% end
% if exist(catmovname,'file')==2
%   fprintf('Combined movie ''%s'' exists; not regenerating.\n',catmovname);
% else
%   catmov(fullfile(expdir,'movie_sde.avi'),fullfile(expdir,'movie_frt.avi'),catmovname);
% end

success = true;