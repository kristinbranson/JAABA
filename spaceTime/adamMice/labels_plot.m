function labels_plot(expdirs,jabs,nframesplot,varargin) 
% Labelgram

% AL 20150128: This is looking more like ethogram_plot 
% AL 20150202: even worse, xlsExport

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

[xlsexport,xlsfile,ethcoreargexport] = myparse(varargin,...
  'xlsexport',false,...
  'xlsfile',[],...
  'ethcoreargexport',false);

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
jab.modernize();

% Special case for AH exps. Call ExpPP.loadexps on all experiments now.
% Any experiments that do not load are excluded.
behs = Labels.verifyBehaviorNames(jab.behaviors.names);
doautomarks = ExpPP.doNamesSpanAllBasicBehaviors(behs);
if doautomarks
  % guess in this case require jabfiles to be nonempty, ie jabs cannot have
  % been obj
  dPP = ExpPP.loadexps(expdirs,jabfiles,'includelabels',true);
  if isempty(dPP)
    error('labels_plot:noExps','ExpPP failed to load any experiments.');
  end
  expsLoaded = {dPP.expfull}';
  tfloaded = ismember(expdirs,expsLoaded);
  if ~all(tfloaded)
    tmpstr = expdirs(~tfloaded);
    tmpstr = sprintf('%s\n',tmpstr{:});
    warningNoTrace('labels_plot:expsExcluded','Excluding %d incomplete experiments from labelgram: %s.',...
      nnz(~tfloaded),tmpstr);
  end
  expdirs = expsLoaded;
  Nexp = numel(expdirs);
end

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


if doautomarks
  % in this case:
  % - use exptags from jab if possible
  
  dPP = ExpPP.ensureCustomGroupInit(dPP); 
  assert(isequal({dPP.expfull}',expdirs));
  [tfExpInJab,loc] = ismember(expdirs,jab.expDirNames);
  if all(tfExpInJab)
    if isstruct(jab) && isfield(jab,'expDirTags') || ...
       isobject(jab) && isprop(jab,'expDirTags')
      exptags = jab.expDirTags(loc);
      exptags = ExperimentTags.cleanLegacyNotTags(exptags);
    else 
      exptags = [];
    end
    if isempty(exptags)
      exptags = ExperimentTags.expTags(expdirs);
    end
    
    allUniqueTags = ExperimentTags.allUniqueTags(exptags);
    allUniqueTagsValid = matlab.lang.makeValidName(allUniqueTags);
    if numel(unique(allUniqueTagsValid))~=numel(allUniqueTagsValid)
      error('ethogram_plot:tag','Repeated tags.'); % very rare
    end
    Ntag = numel(allUniqueTags);
    for iTag = 1:Ntag
      tagRaw = allUniqueTags{iTag};
      cgrpname = allUniqueTagsValid{iTag};
      tfTag = ExperimentTags.findTag(exptags,tagRaw);
      dPP = ExpPP.createCustomGroup(dPP,tfTag,cgrpname);
    end
  end

  ethplotcoreargs = {'doautomarks',true,'automarkdata',dPP,'expppstatprefix','labl','expppcontrolargs',{'labelgram' true}};
else
  ethplotcoreargs = cell(1,0);
end
  
% go
m = Labels.labelMatrix(lbls,T0,T1,behs);
ethplotcoreargs = [{figure,m,T0,T1,expdirs,behs,repmat(jabfiles,numel(behs),1)} ethplotcoreargs];
if ethcoreargexport
  fname = sprintf('ethcoreargs.%s.mat',datestr(now,'yyyymmddTHHMMSS'));
  save(fname,'ethplotcoreargs');  
end
[boutmat,line_names] = ethogram_plot_core(ethplotcoreargs{:});

if xlsexport
  if isempty(xlsfile)
    [xlsfile,xlspath] = uiputfile('bouts.txt','Export to file');
    if isequal(xlsfile,0)
      xlsfile = [];
    else
      xlsfile = fullfile(xlspath,xlsfile);
    end
  elseif exist(xlsfile,'file')
    warning('ethogram_plot:fileExists','Overwriting file ''%s''.',xlsfile);    
  end      
  
  if ~isempty(xlsfile)
    Ethogram.boutExport(xlsfile,boutmat,line_names,behs);
    fprintf(1,'Exported to tab-delimited file ''%s''.\n',xlsfile);
  end
end

