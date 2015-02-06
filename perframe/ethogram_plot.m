function ethogram_plot(expdirs,jabs,nframesplot,varargin)
% ethogram_plot(rootexpdirs,jabfiles,nframesplot,varargin)
% Plots experiment ethograms.
%
% expdirs: cellstr, experiment directories.
% jabs: cellstr, jab filenames, or array of Macguffins.
% nframesplot: Show this many frames in ethogram.
% varargin: optional P-V pairs:
%     
% Double clicking: If you double click on any patch, the corresponding
% experiment and the behavior detector will be opened up in JAABA.
% Shift clicking: Shift-clicking zooms in on an experiment.

% [xlsexport,xlsfile] = myparse(varargin,...
%   'xlsexport',false,...
%   'xlsfile',[]);

if nargin < 1
  expdirs = {};
end
if nargin < 2
  jabs = {};
else
  assert(iscellstr(jabs) || isa(jabs,'Macguffin'));
end
if nargin < 3
  nframesplot = [];
elseif ~isempty(nframesplot)
  assert(isnumeric(nframesplot));
  assert(isscalar(nframesplot) || numel(nframesplot)==2,'Invalid nframesplot');
end

if isempty(expdirs)
  expdirs = uipickfiles('Prompt','Select experiment directories','DirsOnly',true);
  if ~iscellstr(expdirs) || isempty(expdirs)
    return;
  end
end

if isempty(jabs)
  jabs = uipickfiles('Prompt','Select jab files');
  if ~iscellstr(jabs) || isempty(jabs)
    return;
  end
end

% % compile expdirs
% expdirs = {};
% for ndx = 1:numel(expdirs)
%   dd = findAllSubDirs(expdirs{ndx});
%   nn = numel(dd);
%   expdirs(end+1:end+nn) = dd;
% end
% nexps = numel(expdirs);
% if nexps==0
%   warning('ethogram_plot:noExps','No experiments.');
%   return;
% end

[scorefns,behaviornames,jabsforscorefiles,trxfilename] = compileJabInfo(jabs);
assert(isequal(numel(scorefns),numel(behaviornames),numel(jabsforscorefiles)));

[labels,T0,T1] = ...
  formLabelMatrix(expdirs,trxfilename,scorefns,behaviornames,nframesplot);

hfig = figure;
ethplotcoreargs = {hfig,labels,T0,T1,expdirs,behaviornames,jabsforscorefiles};
ethogram_plot_core(ethplotcoreargs{:});

% if xlsexport
%   if isempty(xlsfile)
%     [xlsfile,xlspath] = uiputfile('bouts.txt','Export to file');
%     if isequal(xlsfile,0)
%       xlsfile = [];
%     else
%       xlsfile = fullfile(xlspath,xlsfile);
%     end
%   elseif exist(xlsfile,'file')
%     warning('ethogram_plot:fileExists','Overwriting file ''%s''.',xlsfile);    
%   end      
%   
%   if ~isempty(xlsfile)
%     xlsExport(xlsfile,boutmat,line_names,behaviornames);
%     fprintf(1,'Exported to tab-delimited file ''%s''.\n',xlsfile);
%   end
% end

function [scorefns,behaviornames,jabsout,trxfilename] = compileJabInfo(jabs)
% Compile scorefiles/behaviors from jabfiles.
%
% jabs: cellstr of jabfilenames or cell vec of Macguffins
%
% scorefns: cellstr of (short) score filenames
% behaviornames: cellstr of behaviornames, one per scorefn
% jabsout: cellstr of jabfiles, one per scorefn. eg, jabs{i} is the jabfile
% containing scorefns{i} and behaviornames{i}. jabs may contain repeated
% entries when multi-classifier jabs are present.
% trxfilename: char, unique trxfilename in jabs (all jabs are expected to
% be consistent)

tfLoad = iscellstr(jabs);

Njab = numel(jabs);
scorefns = cell(0,1);
behaviornames = cell(0,1);
jabsout = cell(0,1);
trxfilename = cell(Njab,1);

for iJab = 1:Njab
  if tfLoad
    Q = loadAnonymous(jabs{iJab});
  else
    Q = jabs(iJab);
  end
  
  sfns = Q.file.scorefilename;
  if ischar(sfns)
    sfns = {sfns};
  end
  behs = Labels.verifyBehaviorNames(Q.behaviors.names);
  assert(isequal(numel(sfns),numel(behs)),...
    'Inconsistent scores/behaviors in JAB file.');
  Nbeh = numel(behs);
  scorefns = cat(1,scorefns,sfns(:));
  behaviornames = cat(1,behaviornames,behs(:));
  jabsout = cat(1,jabsout,repmat(jabs(iJab),Nbeh,1));
  
  trxfilename{iJab} = Q.file.trxfilename;
end

trxfilename = unique(trxfilename);
assert(isscalar(trxfilename),'Jabs contain inconsistent trxfilenames.');
trxfilename = trxfilename{1};


function [labels,T0,T1] = formLabelMatrix(expdirs,trxfilename,scorefns,behaviornames,nframesplot)
% Form label matrix using first fly/target in each scorefile.allScores.
%
% labels: nexps x nsamps x nbehaviors array
% T0,T1: first/last frame indices for 2nd dimension of labels, ie
% nsamps=T1-T0+1.
%
% behaviornames just used to check against behaviornames in scorefiles

nexps = numel(expdirs);
nbehaviors = numel(behaviornames);
assert(numel(scorefns)==nbehaviors);
labels = false(nexps,0,nbehaviors);
tStarts = nan(1,nexps);
tEnds = nan(1,nexps);

for i = 1:nexps
  trxfile = fullfile(expdirs{i},trxfilename);
  load(trxfile,'trx');
  tStarts(i) = min([trx.firstframe]);
  tEnds(i) = max([trx.endframe]);

  % load score data  
  allScores = cell(nbehaviors,1);
  for j = 1:nbehaviors,
    scorefile = fullfile(expdirs{i},scorefns{j});
    allScores{j} = load(scorefile);
    if isfield(allScores{j},'behaviorName') && ...
       ~strcmp(allScores{j}.behaviorName,behaviornames{j})
        warningNoTrace('ethogram_plot:behaviorMismatch',...
          'Behavior label mismatch in experiment %s. Scorefile: %s. Expected: %s.',...
          expdirs{i},allScores{j}.behaviorName,behaviornames{j});
    end
    allScores{j} = allScores{j}.allScores;
    if iscell(allScores{j})
      assert(isscalar(allScores{j}));
      allScores{j} = allScores{j}{1};
    end
  end
  
  for j = 1:nbehaviors
    t0 = tStarts(i);
    t1 = min(tEnds(i),numel(allScores{j}.postprocessed{1}));
    if t1 < tEnds(i),
      fprintf('Number of frames for behavior idx%d = %d < number of frames in trajectory file %s = %d\n',...
        j,t1,trxfile,tEnds(i));
      tEnds(i) = t1;
    end    
    
    if t1 > size(labels,2),
      labels = cat(2,labels,nan(nexps,t1-size(labels,2),nbehaviors));
    end
    
    labels(i,t0:t1,j) = allScores{j}.postprocessed{1}(t0:t1);
  end
end

T0 = min(tStarts);
T1 = max(tEnds)-1;

if ~isempty(nframesplot)
  if numel(nframesplot) == 1,
    T1 = min(T1,T0+nframesplot-1);
  else
    assert(numel(nframesplot)==2);
    T0 = nframesplot(1);
    T1 = min(T1,nframesplot(2));
  end
end

labels = labels(:,T0:T1,:);

% function xlsExport(xlsfile,boutmat,expnames,behnames)
% Nexp = numel(expnames);
% Nbeh = numel(behnames);
% assert(isequal(size(boutmat),[Nexp Nbeh]));
% 
% fh = fopen(xlsfile,'w');
% assert(fh~=-1,'Cannot open file ''%s'' for writing.',xlsfile);
% 
% for i = 1:Nexp
%   fprintf(fh,'%s\n',expnames{i});
%   for j = 1:Nbeh
%     fprintf(fh,'\t%s\n',behnames{j});
%     bouts = boutmat{i,j};
%     Nbout = size(bouts,1);
%     assert(size(bouts,2)==2);
%     for k = 1:Nbout
%       fprintf(fh,'\t\t%d\t%d\n',bouts(k,1),bouts(k,2));
%     end
%   end
% end
% 
% fclose(fh);
