function classifyAllDirs(rootdir,jabfile)

if ~isdir(rootdir),
  return;
end

if exist(fullfile(rootdir,'movie.avi'),'file') || isFrontSideExpDir(rootdir)  
  try
    classifyMovie(rootdir,jabfile,'useGetClassifierFromJabFile',true,'doforce',true,'verbose',1);    
  catch ME,
    fprintf('Could not classify %s (%s)\n',rootdir,ME.message);
  end  
else
  
  dd = dir(rootdir);
  for ndx =1:numel(dd)
    if strcmp(dd(ndx).name(1),'.'), continue; end
    curd = fullfile(rootdir,dd(ndx).name);
    if isdir(curd)
      classifyAllDirs(curd,jabfile);
    end
    
  end
end