function testLoading(rootdir)

% rootdir = '/misc/public/Adamtokristin'

matfiles = dir(fullfile(rootdir,'*.mat'));
for i = 1:numel(matfiles),
  matfile = fullfile(rootdir,matfiles(i).name);
  tmp = load(matfile);
  fprintf('Loaded from %s:\n',matfile);
  disp(tmp);
end

subdirs = mydir(rootdir,'isdir',true);
for i = 1:numel(subdirs),
  testLoading(subdirs{i});
end