function [filename,didfind] = GetPCShortcutFileActualPath(filename)

didfind = exist(filename,'file');
if didfind,
  return;
end
if ~ispc,
  return;
end
if ~exist([filename,'.lnk'],'file'),
  return;
end
try
  x = java.io.File([filename,'.lnk']);
  y = sun.awt.shell.ShellFolder.getShellFolder(x);
  actualfilename = y.getLinkLocation();
  filename = char(actualfilename);
  didfind = exist(filename,'file');
catch %#ok<CTCH>
end