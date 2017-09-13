function clearDirectory(rootdir)

if exist(fullfile(rootdir,'movie.avi'),'file') 
  
  rmdir(rootdir,'s');
  
else
  
  dd = dir(rootdir);
  for ndx =1:numel(dd)
    if strcmp(dd(ndx).name(1),'.'), continue; end
    curd = fullfile(rootdir,dd(ndx).name);
    if isdir(curd)
       clearDirectory(curd);
    end
    
  end
  
end