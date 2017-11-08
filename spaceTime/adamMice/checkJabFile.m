function [success,isproblem,msgs,expdirs] = checkJabFile(jabfile,varargin)

[frontside] = myparse(varargin,'frontside',false);

FEATURESFILESTR = 'features.mat';
TRXFILESTR = 'trx.mat';
if frontside,
  MOVIEFILESTR = 'movie_comb.avi';
else
  MOVIEFILESTR = 'movie.avi';
end

Q = load(jabfile,'-mat');
expdirs = Q.x.expDirNames;

isproblem = false(1,numel(expdirs));
msgs = repmat({{}},[1,numel(expdirs)]);
for i = 1:numel(expdirs),
  expdir = expdirs{i};
  if ~exist(expdir,'dir'),
    fprintf('Directory %s does not exist.\n',expdir);
    isproblem(i) = true;
    msgs{i}{end+1} = sprintf('Directory %s does not exist.',expdir);
  end

  file = fullfile(expdir,FEATURESFILESTR);
  if ~exist(file,'file'),
    fprintf('Feature file %s does note exist\n',file);
    isproblem(i) = true;
    msgs{i}{end+1} = sprintf('Feature file %s does note exist',file);
  end
  
  file = fullfile(expdir,TRXFILESTR);
  if ~exist(file,'file'),
    fprintf('Trx file %s does note exist\n',file);
    isproblem(i) = true;
    msgs{i}{end+1} = sprintf('Trx file %s does note exist',file);
  end
  
  file = fullfile(expdir,MOVIEFILESTR);
  if ~exist(file,'file'),
    fprintf('Movie file %s does note exist\n',file);
    isproblem(i) = true;
    msgs{i}{end+1} = sprintf('Movie file %s does note exist',file);
  end
  
end

success = ~any(isproblem);