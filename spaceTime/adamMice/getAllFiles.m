function [flist,dlist] = getAllFiles(rootdir,fname)

dd = dir(rootdir);
dlist = {};
flist = {};

for ndx = 1:numel(dd)
  if dd(ndx).name(1) == '.', continue; end
  if dd(ndx).isdir
    
    if regexp(dd(ndx).name,fname)
      dlist = [dlist fullfile(rootdir,dd(ndx).name)]; %#ok<*AGROW>
    end
    
    [tflist,tdlist] = getAllFiles(fullfile(rootdir,dd(ndx).name),fname);
    flist = [flist tflist];
    dlist = [dlist tdlist];
  else
    if regexp(dd(ndx).name,fname)
      flist = [flist fullfile(rootdir,dd(ndx).name)];
    end
    
    
  end
  
end

