function dirlist = findAllSubDirs(rootdir)
% function dirlist = findAllSubDirs(rootdir)
% Finds all the jaabs subdirectories under the rootdir

dirlist = {};
if exist(fullfile(rootdir,'movie.avi'),'file') || isFrontSideExpDir(rootdir)
  dirlist{end+1} = rootdir;
else
  dd = dir(rootdir);
  for ndx =1:numel(dd)
    if strcmp(dd(ndx).name(1),'.'), continue; end
    curd = fullfile(rootdir,dd(ndx).name);
    if isdir(curd)
      curlist = findAllSubDirs(curd);
      nnew = numel(curlist);
      dirlist(end+1:end+nnew) = curlist;
    end
    
  end
end
