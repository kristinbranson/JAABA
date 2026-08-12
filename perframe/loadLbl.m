function lbl = loadLbl(lbl_file)

tname = tempname;
mkdir(tname);
try
  fprintf('Untarring project into %s\n',tname);
  untar(lbl_file,tname);
  fprintf('... done with untar.\n');
  rawLblFile = fullfile(tname,'label_file.lbl');
catch ME
  if strcmp(ME.identifier,'MATLAB:untar:invalidTarFile') || ...
      strcmp(ME.identifier,'MATLAB:io:archive:untar:invalidTarFile')
    warningNoTrace('Label file %s is not bundled. Using it in raw (mat) format.',lbl_file);
    rawLblFile = lbl_file;
  else
    ME.rethrow();
  end
end

warning('off','MATLAB:load:classNotFound');
warning('off','MATLAB:load:cannotInstantiateLoadedVariable');
lbl = load(rawLblFile,'-mat');
warning('on','MATLAB:load:classNotFound');
warning('on','MATLAB:load:cannotInstantiateLoadedVariable');

% Newer versions of APT save some Labeler properties under the name of their
% private backing field, which has a trailing underscore (e.g. movieFilesAll_
% rather than movieFilesAll).  Rename those to the public property names that
% the rest of JAABA expects, preferring the underscored version when both are
% present.  This is cheap despite rmfield() returning a new struct: only the
% small path-list fields are assigned, and copy-on-write means the big fields
% (labels, movieInfoAll, ...) are shared by reference rather than copied.
fns = fieldnames(lbl);
for i = 1:numel(fns)
  fn = fns{i};
  if numel(fn)>1 && fn(end)=='_'
    lbl.(fn(1:end-1)) = lbl.(fn);
    lbl = rmfield(lbl,fn);
  end
end

[success, message, ~] = rmdir(tname,'s');
if ~success
  error('Could not clear the temp directory %s\n',message);
else
  fprintf('Cleared out temp directory %s\n',tname);
end
