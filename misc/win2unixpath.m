function filepaths = win2unixpath(filepaths,rootdir)

if rootdir(end) ~= '/'
  rootdir(end+1) = '/';
end

filepaths = regexprep(filepaths,'^[A-Z]:\\',rootdir);
filepaths = strrep(filepaths,'\','/');