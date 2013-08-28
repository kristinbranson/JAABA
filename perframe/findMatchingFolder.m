function fname = findMatchingFolder(rootdir,expname)

fname = '';
dd = dir(rootdir);
for ndx = 1:numel(dd);
  if strcmp(dd(ndx).name,expname) && dd(ndx).isdir;
    fname = fullfile(rootdir,expname);
    return;
  end
end

for ndx = 1:numel(dd)
  if dd(ndx).name(1) ==  '.',
    continue;
  end
  fname = findMatchingFolder(fullfile(rootdir,dd(ndx).name),expname);
  if ~isempty(fname)
    return;
  end
end