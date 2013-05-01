
addpath ../misc;
addpath ../filehandling;
addpath params;

%% parameters

%projectfile = '/groups/branson/home/robiea/Projects_data/JAABA/ProjectFiles/Backup.mat';
projectfile = '/groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/experiments/wingextension/WingExtensionJAABAProject.mat';
classifierfile = '/groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/experiments/wingextension/WingExtensionClassifier_KB20130424.mat';

%classifierfile = '/groups/branson/home/robiea/Code_versioned/Jdetect_github/Jdetect/perframe/Backup_AR_v11.mat';
guffdir = '/groups/branson/home/bransonk/behavioranalysis/code/Jdetect/JAABA_guff/perframe';
masterdir = '/groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/perframe';
outguffdir = fullfile(guffdir,'testguff');

%% convert

cd(guffdir);

if ~exist(outguffdir,'dir'),
  mkdir(outguffdir);
end

[~,projectname] = fileparts(projectfile);
outgufffile = fullfile(outguffdir,[projectname,'.jab']);

everythingFileFromOldStyleProjectAndClassifierFiles(...
  outgufffile, ...
  projectfile, ...
  classifierfile, ...
  {});

%% train a classifier using master version

cd(masterdir);

jdmaster = JLabelData(projectfile,'openmovie',false);
classifierdata = load(classifierfile);
jdmaster.classifier_params = classifierdata.classifier_params;
jdmaster.SetExpDirs(classifierdata.expdirs);
jdmaster.Train();

allscores_master = {};
for i = 1:jdmaster.nexps,
  allScores_master{i} = jdmaster.PredictWholeMovie(i);
end

jdmasterstruct = CopyClassProperties(jdmaster,struct);

clear jdmaster;

save(fullfile(outguffdir,sprintf('scores_%s_master.mat',projectname)),'allScores_master','jdmasterstruct');

%% create guff JLabelData

cd(guffdir);

everythingParams=loadAnonymous(outgufffile);
jdguff = JLabelData('macguffin',everythingParams,...
  'setstatusfn',@statusfn_fprintf);

%% compare window data, and reorder if necessary

idxmaster = [jdmasterstruct.windowdata.exp,jdmasterstruct.windowdata.flies,jdmasterstruct.windowdata.t];
idxguff = [jdguff.windowdata.exp,jdguff.windowdata.flies,jdguff.windowdata.t];

nmaster = size(idxmaster,1);
nguff = size(idxguff,1);
if nmaster ~= nguff,
  fprintf('Sizes of windowdata do not match. guff=%d, master=%d\n',nguff,nmaster);
end
n = min(nmaster,nguff);
nmismatch = nnz(any(idxmaster(1:n,:)~=idxguff(1:n,:),2));
if nmismatch == 0,
  fprintf('windowdata matches\n');
end
if nmismatch > 0,
  fprintf('%d index mismatches\n',nmismatch);
  
  [ismguff,idxguff1] = ismember(idxguff,idxmaster,'rows');
  nmissingmaster = nnz(~ismguff);
  if nmissingmaster > 0,
    fprintf('%d indices in guff not in master\n',nmissingmaster);
  end
  
  [ismmaster,idxmaster1] = ismember(idxmaster,idxguff,'rows');
  nmissingguff = nnz(~ismmaster);
  if nmissingguff > 0,
    fprintf('%d indices in master not in guff\n',nmissingguff);
  end
  
  % check that window data is the same for ones that do match
  Xmaster2guff = nan(size(jdguff.windowdata.X));
  Xmaster2guff(ismguff,:) = jdmasterstruct.windowdata.X(idxguff1(ismguff),:);
  nxmismatch1 = nnz(~all(Xmaster2guff(ismguff,:) == jdguff.windowdata.X(ismguff,:),2));
  
  if nxmismatch1 == 0,
    fprintf('X matches after reordering\n');
  else
    fprintf('%d mismatches in X after reordering\n',nxmismatch1);
  end
  
  % check that window data is the same for ones that do match
  Xguff2master = nan(size(jdmasterstruct.windowdata.X));
  Xguff2master(ismmaster,:) = jdguff.windowdata.X(idxmaster1(ismmaster),:);
  nxmismatch2 = nnz(~all(Xguff2master(ismmaster,:) == jdmasterstruct.windowdata.X(ismmaster,:),2));
  
  if nxmismatch2 == 0,
    fprintf('X matches after reordering\n');
  else
    fprintf('%d mismatches in X after reordering\n',nxmismatch2);
  end

  fnscheck = {'labelidx_new','labelidx_imp'};
  nymismatch = 0;
  for k = 1:numel(fnscheck),
    fn = fnscheck{k};
    ymaster2guff = nan(size(jdguff.windowdata.(fn)));
    ymaster2guff(ismguff) = jdmasterstruct.windowdata.(fn)(idxguff1(ismguff));
    nymismatch1 = nnz(ymaster2guff(ismguff,:)~=jdguff.windowdata.(fn)(ismguff));
    if nymismatch1 > 0,
      fprintf('%d mismatches for windowdata.%s\n',nymismatch1,fn);
    else
      fprintf('windowdata.%s matches\n',fn);
    end
    nymismatch = max(nymismatch,nymismatch1);
  end
  
  if nmissingguff == 0 && nmissingmaster == 0 && nxmismatch1 == 0 && nxmismatch2 == 0 && nymismatch == 0,
    
    fprintf('Only ordering is different, rearranging guff windowdata...\n');
    jdguff.windowdata.X = jdguff.windowdata.X(idxmaster1,:);
    jdguff.windowdata.exp = jdguff.windowdata.exp(idxmaster1);
    jdguff.windowdata.flies = jdguff.windowdata.flies(idxmaster1,:);
    jdguff.windowdata.t = jdguff.windowdata.t(idxmaster1);
    jdguff.windowdata.labelidx_cur = jdguff.windowdata.labelidx_cur(idxmaster1);
    jdguff.windowdata.labelidx_new = jdguff.windowdata.labelidx_new(idxmaster1);
    jdguff.windowdata.labelidx_old = jdguff.windowdata.labelidx_old(idxmaster1);
    jdguff.windowdata.labelidx_imp = jdguff.windowdata.labelidx_imp(idxmaster1);
    jdguff.windowdata.predicted = jdguff.windowdata.predicted(idxmaster1);
    jdguff.windowdata.scores = jdguff.windowdata.scores(idxmaster1);
    jdguff.windowdata.scores_old = jdguff.windowdata.scores_old(idxmaster1);
    jdguff.windowdata.scores_validated = jdguff.windowdata.scores_validated(idxmaster1);
    if isfield(jdguff.windowdata,'scores_postprocessed'),
      jdguff.windowdata.scores_postprocessed = jdguff.windowdata.scores_postprocessed(idxmaster1);
    end
    
  end
  

end

%% train guff classifier

jdguff.windowdata.X = single(jdguff.windowdata.X);

jdguff.Train();

allScores_guff = {};
for i = 1:jdguff.nexps,
  allScores_guff{i} = jdguff.PredictWholeMovie(i);
end

jdguffstruct = CopyClassProperties(jdguff,struct);

clear jdguff;

save(fullfile(outguffdir,sprintf('scores_%s_guff.mat',projectname)),'allScores_guff','jdguffstruct');


%% compare scores  

for i = 1:jdmasterstruct.nexps,
  
  nguff = numel(allScores_guff{i}.scores);
  nmaster = numel(allScores_master{i}.scores);
  n = min(nguff,nmaster);
  if nguff ~= nmaster,
    fprintf('Different numbers of flies for experiment %d: guff=%d, master=%d\n',i,nguff,nmaster);
  end
  for j = 1:n,
    nguff1 = numel(allScores_guff{i}.scores{j});
    nmaster1 = numel(allScores_master{i}.scores{j});
    n1 = min(nguff1,nmaster1);
    isnan_guff = isnan(allScores_guff{i}.scores{j}(1:n1));
    isnan_master = isnan(allScores_master{i}.scores{j}(1:n1));
    if nguff1 ~= nmaster1,
      fprintf('Different numbers of frames for experiment %d, fly %d: guff=%d, master=%d\n',i,j,nguff1,nmaster1);
    end
    x = nnz(isnan_guff ~= isnan_master);
    if x > 0,
      fprintf('nan mismatch for experiment %d, fly %d nmismatches = %d\n',i,j,x);
    end
    isgood = ~isnan_guff & ~isnan_master;
    x = max(abs(allScores_guff{i}.scores{j}(isgood)-allScores_master{i}.scores{j}(isgood)));
    if x > 0,
      fprintf('Max mismatch in scores for experiment %d, fly %d = %f\n',i,j,x);
    end
  end
end

%% compare windowdata, to try to figure out why scores didn't match

idxmaster = [jdmasterstruct.windowdata.exp,jdmasterstruct.windowdata.flies,jdmasterstruct.windowdata.t];
idxguff = [jdguffstruct.windowdata.exp,jdguffstruct.windowdata.flies,jdguffstruct.windowdata.t];

nmaster = size(idxmaster,1);
nguff = size(idxguff,1);
if nmaster ~= nguff,
  fprintf('Sizes of windowdata do not match. guff=%d, master=%d\n',nguff,nmaster);
end
n = min(nmaster,nguff);
nmismatch = nnz(any(idxmaster(1:n,:)~=idxguff(1:n,:),2));
if nmismatch > 0,
  fprintf('%d index mismatches\n',nmismatch);
else
  fprintf('windowdata index matches\n');
end

[ismguff,idxguff1] = ismember(idxguff,idxmaster,'rows');
nmissingmaster = nnz(~ismguff);
if nmissingmaster > 0,
  fprintf('%d indices in guff not in master\n',nmissingmaster);
end
  
[ismmaster,idxmaster1] = ismember(idxmaster,idxguff,'rows');
nmissingguff = nnz(~ismmaster);
if nmissingguff > 0,
  fprintf('%d indices in master not in guff\n',nmissingguff);
end

% check that window data is the same for ones that do match
Xmaster2guff = nan(size(jdguffstruct.windowdata.X));
Xmaster2guff(ismguff,:) = jdmasterstruct.windowdata.X(idxguff1(ismguff),:);
nxmismatch1 = nnz(~all(Xmaster2guff(ismguff,:) == jdguffstruct.windowdata.X(ismguff,:),2));

if nxmismatch1 == 0,
  fprintf('X matches after reordering\n');
else
  fprintf('%d mismatches in X after reordering\n',nxmismatch1);
end

% check that window data is the same for ones that do match
Xguff2master = nan(size(jdmasterstruct.windowdata.X));
Xguff2master(ismmaster,:) = jdguffstruct.windowdata.X(idxmaster1(ismmaster),:);
nxmismatch2 = nnz(~all(Xguff2master(ismmaster,:) == jdmasterstruct.windowdata.X(ismmaster,:),2));

if nxmismatch2 == 0,
  fprintf('X matches after reordering\n');
else
  fprintf('%d mismatches in X after reordering\n',nxmismatch2);
end

fnscheck = {'labelidx_new','labelidx_imp'};
nymismatch = 0;
for k = 1:numel(fnscheck),
  fn = fnscheck{k};
  ymaster2guff = nan(size(jdguffstruct.windowdata.(fn)));
  ymaster2guff(ismguff) = jdmasterstruct.windowdata.(fn)(idxguff1(ismguff));
  nymismatch1 = nnz(ymaster2guff(ismguff,:)~=jdguffstruct.windowdata.(fn)(ismguff));
  if nymismatch1 > 0,
    fprintf('%d mismatches for windowdata.%s\n',nymismatch1,fn);
  else
    fprintf('windowdata.%s matches\n',fn);
  end
  nymismatch = max(nymismatch,nymismatch1);
end

%% compare classifiers, to try to figure out why scores didn't match

nguff = numel(jdguffstruct.classifier);
nmaster = numel(jdmasterstruct.classifier);
if nguff ~= nmaster,
  fprintf('Classifier sizes differ. nguff=%d, nmaster=%d\n',nguff,nmaster);
end
n = min(nguff,nmaster);
fns = fieldnames(jdmasterstruct.classifier);
for j = 1:numel(fns),
  ismismatch = false(1,n);
  fn = fns{j};
  for i = 1:n,
    if abs(jdmasterstruct.classifier(i).(fn)-jdguffstruct.classifier(i).(fn)) > 1e-5,
      ismismatch(i) = true;
    end
  end
  if nnz(ismismatch) > 0,
    fprintf('%d mismatches in weak learners for %s, first mismatch at %d\n',nnz(ismismatch),fn,find(ismismatch,1));
  end

end

%% dissect, to try to figure out why scores didn't match

islabeled_guff = (jdguffstruct.windowdata.labelidx_new ~= 0) & (jdguffstruct.windowdata.labelidx_imp);
islabeled_master = (jdmasterstruct.windowdata.labelidx_new ~= 0) & (jdmasterstruct.windowdata.labelidx_imp);

if ~all(islabeled_guff==islabeled_master),
  fprintf('mismatch in islabeled\n');
end

stream = RandStream.getGlobalStream;
reset(stream);

X = single(jdguffstruct.windowdata.X);
[binVals_guff] = findThresholds(X(islabeled_guff,:),jdguffstruct.classifier_params);

stream = RandStream.getGlobalStream;
reset(stream);

[binVals_master] = findThresholds(jdmasterstruct.windowdata.X(islabeled_master,:),jdmasterstruct.classifier_params);

if ~all(binVals_master(:) == binVals_guff(:)),
  fprintf('binVals mismatch\n');
end

bins_master = findThresholdBins(jdmasterstruct.windowdata.X(islabeled_master,:),binVals_master);
bins_guff = findThresholdBins(jdguffstruct.windowdata.X(islabeled_guff,:),binVals_guff);

if ~all(bins_master(:) == bins_guff(:)),
  fprintf('binVals mismatch\n');
end

stream = RandStream.getGlobalStream;
reset(stream);

classifier_guff = boostingWrapper( jdguffstruct.windowdata.X(islabeled_guff,:), ...
  jdguffstruct.windowdata.labelidx_new(islabeled_guff),[],...
  binVals_guff,...
  bins_guff,jdguffstruct.classifier_params);

stream = RandStream.getGlobalStream;
reset(stream);

classifier_master = boostingWrapper( jdmasterstruct.windowdata.X(islabeled_master,:), ...
  jdmasterstruct.windowdata.labelidx_new(islabeled_master),[],...
  binVals_master,...
  bins_master,jdmasterstruct.classifier_params);

n = numel(classifier_master);

for j = 1:numel(fns),
  ismismatch = false(1,n);
  fn = fns{j};
  for i = 1:n,
    if abs(classifier_master(i).(fn)-classifier_guff(i).(fn)) > 1e-5,
      ismismatch(i) = true;
    end
  end
  if nnz(ismismatch) > 0,
    fprintf('%d mismatches in weak learners for %s, first mismatch at %d\n',nnz(ismismatch),fn,find(ismismatch,1));
  end
end

%%
