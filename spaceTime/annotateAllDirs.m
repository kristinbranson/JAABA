function annotateAllDirs(rootdir,varargin)
% Annotates all the jaaba directories within the rootdir with interest
% points. If existing annotations exist, annotateAllDirs asks about reusing
% them. Otherwise, the user is asked to label the first movie and the rest
% of the movies are annotated with the same interest points. This function
% recursively annotates all the subfolders of the rootdir.

[ips,doforce,frontside] = myparse(varargin,...
  'ips',[],... %'firstmovie','',...
  'doforce',false,...
  'frontside',false);

if frontside
  isexpdirfcn = @isFrontSideExpDir;
  hasannfcn = @lclGetAnn;
  genipsfcn = @genInterestPointsFrontSide;
  annfcn = @annotateMiceFrontSide;
else
  isexpdirfcn = @(x)(exist(fullfile(x,'movie.avi'),'file')==2);
  hasannfcn = @lclGetAnn;
  genipsfcn = @(x)genInterestPoints(fullfile(x,'movie.avi'));
  annfcn = @annotateMice;
end

%%% recursively find all experiment dirs
m = depthFirstSearch(rootdir,isexpdirfcn);
alldirs = m.keys;
tf = cell2mat(m.values);
expdirs = alldirs(tf);
fprintf(1,'Found %d experiment directories:\n',numel(expdirs));
cellfun(@(x)fprintf(1,'  %s\n',x),expdirs);
expdirs = expdirs(:)';

assert(~isempty(expdirs),'No experiment directories found.');

%%% check expdirs for existing annotations/interest points
annvals = cellfun(hasannfcn,expdirs,'uni',0);
tfhasann = cellfun(@(x)~isempty(x),annvals);
dirswithanns = expdirs(tfhasann);
dirswoutanns = expdirs(~tfhasann);
anns = annvals(tfhasann);
tfannsexist = ~isempty(anns);
dirs4ips = [dirswoutanns(:);dirswithanns(:)]; % when considering existing annotations, prefer looking at unannotated/new videos

% short-circuit if all directories annotated
if ~doforce && isempty(dirswoutanns)
  fprintf(1,'All directories have existing annotations. Leaving annotation process.\n');
  return;
end

%%% do something hopefully smart
%
% Existing anns exist and user provided ips: check to see if ips matches 
% existing anns (probably it doesn't) and warn if ips is different
%
% Existing anns exist and ~doforce => either
% * reuse existing anns on unannotated movies, or
% * select new ips based on unannotated movies and apply to unannotated movies.
%
% Existing anns exist and doforce => either
% * reuse existing anns on all movies, or
% * select new ips and apply to all movies.
%
% No existing anns =>
% * select new ips and apply to all movies.

if tfannsexist
  fprintf(2,'%d/%d experiments have existing annotations:\n',nnz(tfhasann),numel(expdirs));
  cellfun(@(x,y)fprintf(1,'  %s (food=%s)\n',x,mat2str(round(y.food))),dirswithanns,anns);
  fprintf(1,'%d/%d experiments are unannotated:\n',nnz(~tfhasann),numel(expdirs));
  cellfun(@(x)fprintf(1,'  %s\n',x),dirswoutanns);
  
  if ~isempty(ips) % user provided ips
    if ~doforce 
      tf = cellfun(@(x)isequal(x,ips),anns);
      if ~all(tf)
        uiwait(warndlg('Annotations already exist in some experiments. Some/all of these annotations differ from the interest points you have specified.',...
          'Annotations Exist','modal'))
      end
    end
  else
    ips = lclAskToReuseAnn(dirs4ips,anns,frontside);
    if isempty(ips)
      ips = lclNewIps(dirs4ips,genipsfcn);
    end
  end
else
  if isempty(ips)
    ips = lclNewIps(dirs4ips,genipsfcn);
  end
end
assert(~isempty(ips),'Failed to identify interest points.');

%read first frame of video to get height and width info for annfcn
vidobj = VideoReader(fullfile(expdirs{1},"movie_sde.avi"));
frame1 = read(vidobj,1);
imheight = size(frame1,1);
imwidth  = size(frame1,2);

% Apply ips to expDirs
for edir = expdirs, edir = edir{1}; %#ok<FXSET>
  try
    
    tf = annfcn(edir,'ips',ips,'doforce',doforce,'imwidth',imwidth); % 'fno',30
    assert(tf);
  catch ME
    fprintf('Could not annotate %s (%s)\n',edir,getReport(ME));
  end
  fprintf(1,'Annotated exp: %s\n',edir);  
end

function ips = lclAskToReuseAnn(expdirs,anns,frontside)
tfannsareallsame = all(cellfun(@(x)isequal(x,anns{1}),anns));
if tfannsareallsame
  uiwait(msgbox('A single set of consistent annotations already exists. These will be displayed.',...
    'Annotations exist','modal'));
  firstexp = expdirs{1}; % for now just display on top of first movie
  foundips = anns{1};
  
  if frontside
    firstmovsde = fullfile(firstexp,'movie_sde.avi');
    firstmovfrt = fullfile(firstexp,'movie_frt.avi');
    firstmovcmb = fullfile(firstexp,'movie_comb.avi');
    
    if exist(firstmovsde,'file')==2 && exist(firstmovfrt,'file')==2
      foundipssde = rmfield(foundips,{'foodfront' 'mouthfront'});
      foundipsfrt = rmfield(foundips,{'food' 'mouth' 'perch'});
      
      htmp = [];
      htmp(1) = playfmf_interestpoints('moviename',firstmovsde,'ipsshow',foundipssde);
      htmp(2) = playfmf_interestpoints('moviename',firstmovfrt,'ipsshow',foundipsfrt);
      uiwaitvec(htmp);
    else
      T = load(fullfile(firstexp,'trx.mat'));
      imwidth = T.trx(1).concatmov.imwidth; % this info should exist as combined movie has been created in past
      foundipscmb = foundips;
      foundipscmb.foodfront(1) = foundipscmb.foodfront(1)+imwidth;
      foundipscmb.mouthfront(1) = foundipscmb.mouthfront(1)+imwidth;
      uiwait(playfmf('moviename',firstmovcmb,'ipsshow',foundipscmb));
    end
  else
    firstmov = fullfile(firstexp,'movie.avi');
    playfmf('moviename',firstmov,'ipsshow',foundips);
  end
  
  resp = questdlg('Do you want to reuse existing annotation?','Existing annotations','Yes','No','Cancel','No');
  if isempty(resp)
    resp = 'Cancel';
  end
  switch resp
    case 'Yes'
      ips = foundips;
    case 'No'
      ips = [];
    case 'Cancel'
      error('User canceled.');
  end
else 
  uiwait(warndlg('Multiple, different annotations already exist in this experiment tree. You will be generating yet another set of annotation points.',...
    'Annotations exist','modal'));
  ips = [];
end

function ips = lclNewIps(expdirs,genipsfcn)
% get interest points by trying movies in succession

for edir = expdirs(:)', edir = edir{1}; %#ok<FXSET>
  [tf ips] = genipsfcn(edir);
  if tf % success
    assert(~isempty(ips));
    break;
  end
end

function ips = lclGetAnn(edir)
ips = [];

trxfile = fullfile(edir,'trx.mat');
if exist(trxfile,'file')==2
  T = load(trxfile);
  if isfield(T.trx,'arena')
    ips = T.trx(1).arena;
  end
end

