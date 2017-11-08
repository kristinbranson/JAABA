function expdirs = genAllFeaturesExpDirList(rootdir,varargin)
% function expdirs = genAllFeaturesExpDirList(rootdir,varargin)
% Generates features for all the jaaba folder in rootdir.

[doforce,frontside,~] = myparse_nocheck(varargin,...
  'doforce',false,'frontside',false);

expdirs = {};

if ~isdir(rootdir),
  return;
end


if frontside,
  movieexists = exist(fullfile(rootdir,'movie_sde.avi'),'file') && ...
    exist(fullfile(rootdir,'movie_frt.avi'),'file');
else
  movieexists = exist(fullfile(rootdir,'movie.avi'),'file');
end



if movieexists,
  
  try
    
    curf  = rootdir;
    featuresfile = fullfile(curf,'features.mat');
    trxfile = fullfile(curf,'trx.mat');
    ff = dir(featuresfile);
    tt = dir(trxfile);
    if exist(featuresfile,'file') && ...
        exist(trxfile,'file') && ...
        ff(1).datenum > tt(1).datenum && ...
        ~doforce
      fprintf('Features have been computed for %s, skipping.\n',rootdir);
      return;
    end
    
    expdirs = {rootdir};
        
  catch ME,
    fprintf('Could not find all files for %s (%s)\n',rootdir,getReport(ME));
  end
  
else
  
  dd = dir(rootdir);
  for ndx =1:numel(dd)
    if strcmp(dd(ndx).name(1),'.'), continue; end
    curd = fullfile(rootdir,dd(ndx).name);
    if isdir(curd)
       expdirs1 = genAllFeaturesExpDirList(curd,varargin{:});
       expdirs = [expdirs,expdirs1]; %#ok<AGROW>
    end
    
  end
  
end