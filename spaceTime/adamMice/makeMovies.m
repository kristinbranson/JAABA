function expdirs = makeMovies(rootdir,varargin)

[doforce,dosoftlink,~,rotatesideview] = ...
  myparse(varargin,'doforce',false,...
  'dosoftlink',false,'rotatefrontview',0,'rotatesideview',0);
% function makeMovies(rootdir)
% makes a jaaba folder for all the avi's nested in the rootdir.

expdirs = {};

dd = dir(rootdir);
for ndx =1:numel(dd)
  if strcmp(dd(ndx).name(1),'.'), continue; end
  curd = fullfile(rootdir,dd(ndx).name);
  if isdir(curd)
    
    if exist([curd '.avi'],'file') && exist(fullfile(curd,'movie.avi'),'file'),
      continue; 
    else
      expdirs1 = makeMovies(curd);
      expdirs = [expdirs,expdirs1]; %#ok<AGROW>
    end
    
  else
    
    [~,fname,ext] = fileparts(dd(ndx).name);
    if strcmpi(ext,'.avi')
      
      if ismember(dd(ndx).name,{'movie.avi','movie_sde.avi','movie_frt.avi'}),
        fprintf('Not creating directory for %s, %s is a reserved name.\n',fullfile(rootdir,dd(ndx).name),dd(ndx).name);
        continue;
      end
      
      m = regexp(dd(ndx).name,'_front_','once');
      if ~isempty(m),
        fprintf('%s is a front-facing movie, not creating a directory for it.\n',fullfile(rootdir,dd(ndx).name));
        continue;
      end

      
      curf  = fullfile(rootdir,fname);
      if ~exist(curf,'dir')
        mkdir(curf);
      end
%       scurd = regexprep(curd,' ','\\ ');
%       scurf = regexprep(curf,' ','\\ ');
      if exist(fullfile(curf,'movie.avi'),'file') && doforce,
        delete(fullfile(curf,'movie.avi'));
      end
      if ~exist(fullfile(curf,'movie.avi'),'file'),
        try
          if rotatesideview ~= 0,
            rotatemov(curd,fullfile(curf,'movie.avi'),rotatesideview);
          else
            if isunix && ~ismac && nargin > 1 && dosoftlink,
              fprintf('Soft-linking %s to %s...\n',curd,fullfile(curf,'movie.avi'));
              [~,n,e] = fileparts(curd);
              % use relative links so that maybe this works when mounting from windows/mac?
              n = [n,e];
              syscmd = sprintf('ln -s ''../%s'' ''%s/movie.avi''',n,curf);
              unix(syscmd);
            else
              fprintf('Copying %s to %s...\n',curd,fullfile(curf,'movie.avi'));
              copyfile(curd,fullfile(curf,'movie.avi'));
            end
          end
        catch ME,
          fprintf('Could not copy video file for %s (%s)\n',curf,ME.message);
          continue;
        end
      else
        fprintf('%s exists, not re-copying it.\n',fullfile(curf,'movie.avi'));
      end
      if ~exist(fullfile(curf,'trx.mat'),'file')
        try
          fprintf('Generating trx file %s...\n',fullfile(curf,'trx.mat'));
          [~,nframes,~,~] = get_readframe_fcn(fullfile(curf,'movie.avi'));
          T.trx = genTrack(nframes,[640 480]);
          save(fullfile(curf,'trx.mat'),'-struct','T');
        catch ME,
          fprintf('Could not create trx file for %s (%s)\n',curf,ME.message);
          continue;
        end
      else
        fprintf('Trx file exists for %s, not regenerating\n',curf);
      end
      expdirs{end+1} = curf; %#ok<AGROW>
    end
    
  end
  
end