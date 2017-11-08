function [success,ips] = annotateMice(expdir,varargin)

success = false; 

%%
[ips,doforce] = myparse(varargin,...
  'ips',[],...
  'doforce',false);

T = load(fullfile(expdir,'trx.mat'));

if ~doforce,  
  if isfield(T,'trx') && ...
      isfield(T.trx(1),'arena') && ...
      all(isfield(T.trx(1).arena,{'mouth','food','perch','face'})),
    success = true;
    fprintf('Skipping annotation for %s, annotations exist.\n',expdir);
    return;
  end
end

if isempty(ips)
  moviename = fullfile(expdir,'movie.avi');
  [tf ips] = genInterestPoints(moviename);
  if ~tf
    warning('annotateMice:annFailed','Failed to generate interest points for ''%s''.',moviename);
    ips = [];
    return;
  end
end

%%
T.trx(1).arena = ips;
T.trx(2) = T.trx(1);
T.trx(3) = T.trx(1);
T.trx(1).x(:) = T.trx(1).arena.food(1);
T.trx(1).y(:) = T.trx(1).arena.food(2);
T.trx(2).x(:) = T.trx(1).arena.mouth(1);
T.trx(2).y(:) = T.trx(1).arena.mouth(2);
T.trx(3).x(:) = T.trx(1).arena.perch(1);
T.trx(3).y(:) = T.trx(1).arena.perch(2);

save(fullfile(expdir,'trx.mat'),'-struct','T');

success = true;
