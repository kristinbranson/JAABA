function labels_plot(expdirs,jabs,nframesplot,varargin) 
% Labelgram

if ~exist('expdirs','var')
  expdirs = {};
end
assert(iscellstr(expdirs));
if ~exist('jabs','var')
  jabs = [];
end
if isempty(jabs)
  [tfsuccess,jabs] = ExpPP.uiGetJabFiles;
  if ~tfsuccess
    return;
  end
end
if ~exist('nframesplot','var')
  nframesplot = [];
elseif ~isempty(nframesplot)
  assert(isnumeric(nframesplot));
  assert(isscalar(nframesplot) || numel(nframesplot)==2,'Invalid nframesplot');
end

assert(isscalar(jabs),'Currently expect single jab/macguffin.');
if iscellstr(jabs)
  jabfiles = jabs;
  jabs = loadAnonymous(jabs{1});
elseif isa(jabs,'Macguffin');
  jabfiles = {[]};
else
  assert(false,'''jabs'' argument must be either a jab filename or object.');
end

if isempty(expdirs)
  expdirs = jabs.expDirNames;
end
Nexp = numel(expdirs);

jab = jabs;
%jab.modernize();

% get the labels for expdirs
assert(numel(jab.expDirNames)==numel(jab.labels));
[tf,loc] = ismember(expdirs,jab.expDirNames);
lbls = Labels.labels(Nexp);
lbls(tf) = jab.labels(loc(tf));

% Trx/T0/T1/nframesplot
trxfname = jab.file.trxfilename;
trxstarts = nan(Nexp,1);
trxends = nan(Nexp,1);
for i = 1:Nexp
  trxfile = fullfile(expdirs{i},trxfname);
  tmp = load(trxfile,'trx');
  trxstarts(i) = min([tmp.trx.firstframe]);
  trxends(i) = max([tmp.trx.endframe]);
end
T0 = min(trxstarts);
T1 = max(trxends)-1;
if ~isempty(nframesplot)
  if numel(nframesplot)==1,
    T1 = min(T1,T0+nframesplot-1);
  else
    assert(numel(nframesplot)==2);
    T0 = nframesplot(1);
    T1 = min(T1,nframesplot(2));
  end
end

% go
behnames = Labels.verifyBehaviorNames(jab.behaviors.names);
m = Labels.labelMatrix(lbls,T0,T1,behnames);
ethogram_plot_core(figure,m,T0,T1,expdirs,behnames,repmat(jabfiles,numel(behnames),1));
