function ips = genAllFeatures(rootdir,varargin)
% function ips = genAllFeatures(rootdir,varargin)
% Generates features for all the jaaba folder in rootdir.

ips = [];

if isdeployed,
  fprintf('Root directory: %s\n',rootdir);
  fprintf('Optional arguments:\n');
  fprintf('%s\n',varargin{:});
end

[doforce,dorecurse,framechunksize,frontside] = myparse(varargin,...
  'doforce',false,'dorecurse',true,'framechunksize',[],'frontside',false);

if ischar(doforce),
  doforce = str2double(doforce);
end

if ischar(dorecurse),
  dorecurse = str2double(dorecurse);
end

if ischar(frontside),
  frontside = str2double(frontside);
end

if ~isdir(rootdir),
  return;
end

if ~frontside && exist(fullfile(rootdir,'movie.avi'),'file') || ...
    frontside && isFrontSideExpDir(rootdir)
  
  try    
    curf  = rootdir;
    featuresfile = fullfile(curf,'features.mat');
    trxfile = fullfile(curf,'trx.mat');
    
    if frontside
      vidfile.front = fullfile(curf,'movie_frt.avi');
      vidfile.side = fullfile(curf,'movie_sde.avi');
      inputvidsexist = exist(vidfile.front,'file') && exist(vidfile.side,'file');

    else
      vidfile = fullfile(curf,'movie.avi');            
      inputvidsexist = exist(vidfile,'file');
    end
    
    
    ff = dir(featuresfile);
    tt = dir(trxfile);
    if exist(featuresfile,'file') && ...
        exist(trxfile,'file') && ...
        ff(1).datenum > tt(1).datenum && ...
        ~(doforce && inputvidsexist),
      fprintf('Features already computed for %s, skipping.\n',rootdir);
    else
      if exist(featuresfile,'file')
        delete(featuresfile);
      end
    
      trxfile = fullfile(curf,'trx.mat');
      T = load(trxfile);
      face = T.trx(1).arena.face;
    
      [F,H,info] = genFeatures(...
        'vidfile',vidfile,...
        'trxfile',trxfile,...
        'face',face,...
        'framechunksize',framechunksize,...
        'frontside',frontside);
            
      F = reshape(F,[],info.nframes)';
      H = reshape(H,[],info.nframes)';
      curFeatures = [F H]; %#ok<NASGU>
      
      save(fullfile(curf,'features'),'curFeatures');
      fprintf('Finished generating features for %s.\n',rootdir);
    end
    
    if frontside,
      
      %%% Create concatenated movie for JAABA playback
      
      COMBND_MOVIE_NAME = 'movie_comb.avi';
      
      fprintf(1,'Concatenating front/side movies for playback for %s...\n',rootdir);
      catmovname = fullfile(rootdir,COMBND_MOVIE_NAME);
      if exist(catmovname,'file')==2 && doforce
        delete(catmovname);
      end
      if exist(catmovname,'file')==2
        fprintf('Combined movie ''%s'' exists; not regenerating.\n',catmovname);
      else
        catmov(vidfile.side,vidfile.front,catmovname);
      end
      
      if exist(catmovname,'file'),
              
        % delete front/side movies (they must exist, we just used them).
        % Consider not copying them in first place.
        delete(vidfile.side);
        delete(vidfile.front);
        
      else
        fprintf('Combined movie %s does not exist!\n',catmovname);
      end
    end
    
    
  catch ME,
    fprintf('Could not generate features for %s (%s)\n',rootdir,getReport(ME));
  end
  
else
  
  if dorecurse == 0,
    return;
  end
  
  dd = dir(rootdir);
  for ndx =1:numel(dd)
    if strcmp(dd(ndx).name(1),'.'), continue; end
    curd = fullfile(rootdir,dd(ndx).name);
    if isdir(curd)
       genAllFeatures(curd,varargin{:});
    end
    
  end
  
end