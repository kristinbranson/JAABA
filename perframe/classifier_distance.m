
%% Set the classifier and find a bagging model.

SetUpJAABAPath;
clfile = '/groups/branson/home/robiea/Code_versioned/Jdetect_github/Jdetect/perframe/ChaseAR_newv9.mat';
cl = load(clfile);
configfile = cl.configfilename;
data = JLabelData(configfile,'openmovie',false);
data.SetClassifierFileName(clfile);

%% Keep aside an experiment
exp2rem = 4;
oldlabels = data.windowdata.labelidx_new;
idx2rem = data.windowdata.exp==exp2rem;
data.windowdata.labelidx_new(idx2rem) = 0;
%% Do normal bagging to normalize scores

bouts = data.getLabeledBouts();
bouts.ndx = bouts.ndx>0.5;

islabeled = data.windowdata.labelidx_new ~=0;
binVals = findThresholds(data.windowdata.X(islabeled,:),data.classifier_params);

[oldBagModels, ~] =...
  doBaggingBouts( data.windowdata.X, ...
  data.windowdata.labelidx_new,data,...
  binVals,...
  data.classifier_params,bouts);

bscores = zeros(size(data.windowdata.X,1),numel(oldBagModels));
for ndx = 1:numel(oldBagModels)
  bscores(:,ndx) = myBoostClassify(data.windowdata.X,oldBagModels{ndx});
end

tot = 0;
for ndx = 1:numel(oldBagModels)
  tot = tot+sum(abs([oldBagModels{ndx}(:).alpha]));
  
end

stdScores = std(bscores,[],2);

%% For fast bagging stuff..
islabeled = data.windowdata.labelidx_new ~= 0;
bins = findThresholdBins(data.windowdata.X(islabeled,:),binVals);
bmodel = fastBag(data.windowdata.X(islabeled,:),...
  data.windowdata.labelidx_new(islabeled),...
  binVals,bins,data.classifier_params);

distmatfast = BagDistMat(data.windowdata.X,{bmodel});

%% Find errors in the experiment that was kept aside
numTrain = 25;
oldCls = {};
oldScores = [];
for ndx = 1:numTrain;
  data.Train();
  oldCls{ndx} = data.classifier;
  oldScores(:,ndx) = myBoostClassify(data.windowdata.X,oldCls{ndx});
end

%%

modLabels = sign(1.5-oldlabels);
nt = size(oldScores,2);
mistakes = sum( (oldScores.*repmat(modLabels,1,nt))<0,2);
mostmistakes = find( mistakes>=nt & idx2rem);

%% Select all

mostmistakes = find(idx2rem);

%% 
 
numrep = 1;
numTrain = 25;
scores_slow = zeros(numTrain,numrep);
scores_fast = zeros(numTrain,numrep);
allselex = [];
dist_slow = {};
dist_fast = {};
for repndx = 1:numrep
  
data.windowdata.labelidx_new = oldlabels;
data.windowdata.labelidx_new(idx2rem) = 0;


%% Compute the distance of the frames with most mistakes to the training.

selex = randsample(mostmistakes,1);
allselex(end+1) = selex;
mean_dist = zeros(size(data.windowdata.X,1),1);
for bndx = 1:numel(oldBagModels)
  distMat = BagDistMat(data.windowdata.X,oldBagModels(bndx));
  cur_dist = sum(abs(distMat - repmat(distMat(selex,:),size(distMat,1),1) ),2);
  mean_dist = mean_dist+ cur_dist;
end
max_train = max(mean_dist(~idx2rem));
mean_dist(idx2rem) = max_train;


% Find the training bout that have the closest ex..

[sort_dist training_ord] = sort(mean_dist);

nbouts = numel(bouts.label);
bout_dist = inf(1,nbouts);
selndx = false(size(mean_dist));
ndx = 0; count = 0;
while count < 20,
  ndx = ndx + 1;
  if oldlabels(training_ord(ndx))==oldlabels(selex), 
    continue;
  end
  curbout = bouts.ndx(:,training_ord(ndx));
  if nnz(curbout) == 0, continue; end % Bout is unimportant.
  selndx = selndx | bouts.ndx(curbout,:)';
  count = nnz(selndx);
end

sel2change = find(selndx);

% 
% for ndx = 1:numel(training_ord),
%   if oldlabels(training_ord(ndx)) == oldlabels(selex),
%     continue;
%   else
%     sel2change = training_ord(ndx);
%     break;
%   end
% end


mistakes(sel2change)'
dist_slow{repndx} = mean_dist(sel2change);


%
idx = false(size(data.windowdata.labelidx_new));

idx(sel2change) = true;

newCls = {}; newScores = [];
data.windowdata.labelidx_new(idx) = 3 -data.windowdata.labelidx_new(idx);
% hf = figure;
for ndx = 1:numTrain
  data.Train();
  newCls{ndx} = data.classifier;
  cursc = myBoostClassify(data.windowdata.X,newCls{ndx});
  newScores(:,ndx) = cursc;
  scores_slow(repndx,ndx) = cursc(selex);
  normdiff = (mean(newScores,2) - mean(oldScores,2))./stdScores;
  newmistakes = sum( (newScores.*repmat(modLabels,1,ndx))<0,2);
%   sttnew = sprintf('%.3f ',newmistakes(sel2change)/ndx);
%   sttold = sprintf('%.3f ',mistakes(sel2change)/numTrain);
%   figure(hf); 
%   subplot(2,2,1:2),hist(normdiff);
%   title( {sprintf('Iter:%d dist:%.3f',ndx,normdiff(selex)),...
%         sprintf('new mistakes:%s ',sttnew),...
%         sprintf('oldmistakes:%s',sttold)});
%   subplot(2,2,3); hist(normdiff(idx2rem))
%   subplot(2,2,4); hist([oldScores(selex,1:ndx);newScores(selex,1:ndx)]');
end
data.windowdata.labelidx_new(idx) = 3 -data.windowdata.labelidx_new(idx);

%% Reset labels
data.windowdata.labelidx_new = oldlabels;
data.windowdata.labelidx_new(idx2rem) = 0;

%% For fast bagging.
mean_distfast = sum(abs(distmatfast-repmat(distmatfast(selex,:),size(distmatfast,1),1)),2);
max_train = max(mean_distfast);
mean_distfast(idx2rem) = max_train;

% Find the training bout that have the closest ex..

[sort_dist training_ord] = sort(mean_distfast);

nbouts = numel(bouts.label);
bout_dist = inf(1,nbouts);
selndx = false(size(mean_dist));
ndx = 0; count = 0;
while count < 40,
  ndx = ndx + 1;
  if oldlabels(training_ord(ndx))==oldlabels(selex), 
    continue;
  end
  curbout = bouts.ndx(:,training_ord(ndx));
  if nnz(curbout) == 0, continue; end % Bout is unimportant.
  selndx = selndx | bouts.ndx(curbout,:)';
  count = nnz(selndx);
end

sel2change = find(selndx);
% 
% for ndx = 1:numel(training_ord),
%   if oldlabels(training_ord(ndx)) == oldlabels(selex),
%     continue;
%   else
%     sel2change = training_ord(ndx);
%     break;
%   end
% end


mistakes(sel2change)'

dist_fast{repndx} = mean_distfast(sel2change);

%
idx = false(size(data.windowdata.labelidx_new));

idx(sel2change) = true;

newCls = {}; newScores = [];
data.windowdata.labelidx_new(idx) = 3 -data.windowdata.labelidx_new(idx);
% hf = figure;
for ndx = 1:numTrain
  data.Train();
  newCls{ndx} = data.classifier;
  cursc = myBoostClassify(data.windowdata.X,newCls{ndx});
  newScores(:,ndx) = cursc;
  scores_fast(repndx,ndx) = cursc(selex);
  normdiff = (mean(newScores,2) - mean(oldScores,2))./stdScores;
  newmistakes = sum( (newScores.*repmat(modLabels,1,ndx))<0,2);
%   sttnew = sprintf('%.3f ',newmistakes(sel2change)/ndx);
%   sttold = sprintf('%.3f ',mistakes(sel2change)/numTrain);
%   figure(hf); 
%   subplot(2,2,1:2),hist(normdiff);
%   title( {sprintf('Iter:%d dist:%.3f',ndx,normdiff(selex)),...
%         sprintf('new mistakes:%s ',sttnew),...
%         sprintf('oldmistakes:%s',sttold)});
%   subplot(2,2,3); hist(normdiff(idx2rem))
%   subplot(2,2,4); hist([oldScores(selex,1:ndx);newScores(selex,1:ndx)]');
end
data.windowdata.labelidx_new(idx) = 3 -data.windowdata.labelidx_new(idx);
data.windowdata.labelidx_new = oldlabels;
data.windowdata.labelidx_new(idx2rem) = 0;

end

%%
mdistfast = [];
for ndx = 1:numel(dist_fast),
  zz = sort(dist_fast{ndx},'ascend');
  mdistfast(ndx) = zz(1);
end

mdistslow = [];
for ndx = 1:numel(dist_slow),
  zz = sort(dist_slow{ndx},'ascend');
  mdistslow(ndx) = zz(1);
end
mfast = mean(scores_fast(:,1:numTrain),2) - mean(oldScores(allselex,:),2);
mslow = mean(scores_slow(:,1:numTrain),2) - mean(oldScores(allselex,:),2);

figure; 
ff1 = mfast.*sign(1.5 - oldlabels(allselex));
ff2 = mslow.*sign(1.5 - oldlabels(allselex));
mm1 = min([ff1(:);ff2(:)]);
mm2 = max([ff1(:);ff2(:)]);
subplot(2,2,1:2); scatter(ff1,ff2,'.');
hold on; plot(mm1:mm2,mm1:mm2,'r'); axis equal;
subplot(2,2,3); scatter(mdistfast(oldlabels(allselex)==2),mfast(oldlabels(allselex)==2),'.');
subplot(2,2,4); scatter(mdistslow(oldlabels(allselex)==2),mslow(oldlabels(allselex)==2),'.');


%% Compare 2 distances..

for ndx = 1:numel(bmodel);
  bmodel(ndx).alpha = bmodel(ndx).newalpha;
end
distmatfast_new = BagDistMat(data.windowdata.X,{bmodel});

for ndx = 1:numel(bmodel);
  bmodel(ndx).alpha = bmodel(ndx).oldalpha;
end
distmatfast_old = BagDistMat(data.windowdata.X,{bmodel});

%%
z = figure;
corr_newalpha = [];
corr_oldalpha = [];

for rep = 1:5
  selext = randsample(mostmistakes,1);
  mean_dist = zeros(size(data.windowdata.X,1),1);
  for bndx = 1:numel(oldBagModels)
    distMat = BagDistMat(data.windowdata.X,oldBagModels(bndx));
    cur_dist = sum(abs(distMat - repmat(distMat(selext,:),size(distMat,1),1) ),2);
    mean_dist = mean_dist+ cur_dist;
  end

  mean_distfast = sum(abs(distmatfast_new-repmat(distmatfast_new(selext,:),size(distmatfast_new,1),1)),2);
  
  mean_dist(idx2rem) = [];
  mean_distfast(idx2rem) = [];
  figure(z); subplot(2,5,rep); scatter(mean_dist,mean_distfast,'.');
  title(sprintf('Corrcoef: %5f',corr(mean_dist,mean_distfast)));
  corr_newalpha(rep) = corr(mean_dist,mean_distfast);

  mean_distfast = sum(abs(distmatfast_old-repmat(distmatfast_old(selext,:),size(distmatfast_old,1),1)),2);
  mean_distfast(idx2rem) = [];
  figure(z); subplot(2,5,5+rep); scatter(mean_dist,mean_distfast,'.');
  title(sprintf('Corrcoef: %5f',corr(mean_dist,mean_distfast)));
  corr_oldalpha(rep) = corr(mean_dist,mean_distfast);
end

figure; scatter(corr_newalpha,corr_oldalpha,'.'); hold on;
plot(0:0.5:1,0:0.5:1,'r');

%%

z = figure;
corr_newalpha = [];

for rep = 1:4
  selext = randsample(mostmistakes,1);
  mean_dist = zeros(size(data.windowdata.X,1),1);
  for bndx = 1:numel(oldBagModels)
    distMat = BagDistMat(data.windowdata.X,oldBagModels(bndx));
    cur_dist = sum(abs(distMat - repmat(distMat(selext,:),size(distMat,1),1) ),2);
    mean_dist = mean_dist+ cur_dist;
  end

  mean_distfast = sum(abs(distmatfast_new-repmat(distmatfast_new(selext,:),size(distmatfast_new,1),1)),2);
  
  mean_dist(idx2rem) = [];
  mean_distfast(idx2rem) = [];
  figure(z); subplot(2,2,rep); scatter(mean_dist,mean_distfast,'.');
  title(sprintf('Corrcoef: %5f',corr(mean_dist,mean_distfast)));
  corr_newalpha(rep) = corr(mean_dist,mean_distfast);
end


%%
for ndx = 1:numel(bmodel);
bmodel(ndx).alpha = bmodel(ndx).oldalpha;
end
%%
for ndx = 1:numel(bmodel);
bmodel(ndx).alpha = bmodel(ndx).newalpha;
end

%%
for ndx = 1:numel(bmodel);
bmodel(ndx).newalpha = bmodel(ndx).alpha;
end



%% New features.

modLabels = sign(1.5-oldlabels);
nt = size(oldScores,2);
mistakes = sum( (oldScores.*repmat(modLabels,1,nt))<0,2);
mostmistakes = find( mistakes>=nt & idx2rem);

ex_order = zeros(numel(oldlabels),numel(mostmistakes));
ex_dist =  zeros(numel(oldlabels),numel(mostmistakes));
for ndx = 1:numel(mostmistakes)
selex = mostmistakes(ndx);
mean_dist = zeros(size(data.windowdata.X,1),1);
for bndx = 1:numel(oldBagModels)
  distMat = BagDistMat(data.windowdata.X,oldBagModels(bndx));
  cur_dist = sum(abs(distMat - repmat(distMat(selex,:),size(distMat,1),1) ),2);
  mean_dist = mean_dist+ cur_dist;
end
max_train = max(mean_dist(~idx2rem));
mean_dist(idx2rem) = max_train;
[sort_dist training_ord] = sort(mean_dist);
ex_order(training_ord,ndx) = 1:numel(oldlabels);
ex_dist(:,ndx) = mean_dist;
if mod(ndx,10)==0, fprintf('.'); end
end
fprintf('\n');


%%
ii = ex_order<100;
zz = sum(ii,2);
ex_train = find(zz>30);
tt = sum(ii(ex_train,:));
ex_mistakes = mostmistakes(find(tt>9));

yet_another_var = sum(ii(ex_train,tt>9),2);
ex_train(yet_another_var< numel(ex_train)/2) = [];

dd_mistakes = BagDistMat(data.windowdata.X(ex_mistakes,:),oldBagModels);
dd_train = BagDistMat(data.windowdata.X(ex_train,:),oldBagModels);

label_mistakes = 2;
train_same = find( (oldlabels==label_mistakes) & (oldlabels>0) & ~idx2rem);
rand_same = randsample(train_same,20);
dd_rand_same = BagDistMat(data.windowdata.X(rand_same,:),oldBagModels);
ss_rand_same = mean(sign(dd_rand_same),1);

label_mistakes = 2;
train_opposite = find( (oldlabels~=label_mistakes) & (oldlabels>0) & ~idx2rem);
rand_train = randsample(train_opposite,20);
dd_rand_opp = BagDistMat(data.windowdata.X(rand_train,:),oldBagModels);
ss_rand_opp = mean(sign(dd_rand_opp),1);

ss_mistakes = mean(sign(dd_mistakes),1);
ss_train = mean(sign(dd_train),1);
wt = abs(dd_mistakes(1,:));

%%

f_err = abs(ss_mistakes-ss_train)>1.5;
f_same = abs(ss_rand_same - ss_mistakes)<0.10;
f_opp = abs(ss_rand_opp - ss_mistakes)>1.50;

f_sel = f_err & f_same & f_opp;
sum(wt(f_sel))/sum(wt)


%%

numBins = size(binVals,1)+1;
% for histogramming
edges = .5:numBins+.5;
% indices with positive labels
% always normalize histograms by Z


numS = data.classifier_params.numSample;

selX = data.windowdata.labelidx_new ~=  0;
labels = sign(1.5- data.windowdata.labelidx_new(selX));
curX = data.windowdata.X(selX,:);
curbins = bins;

numDim = size(data.windowdata.X,2);
numEx = size(curX,1); 

origScores = myBoostClassify(curX,oldCls{1});
tt = origScores.*labels;
initWt = getWeights(labels);
dist = initWt./(1+exp(tt));
dist = dist/sum(dist);

curSel = randsample(numEx,numS,true,dist);

allBinErr = zeros(numBins,numDim);
allDir = zeros(numBins,numDim);

allSelErr = zeros(numBins,numDim);
subX_cur = data.windowdata.X(ex_mistakes,:);
subX_train = data.windowdata.X(ex_train,:);
num_mistakes = numel(ex_mistakes);
num_train = numel(ex_train);

parfor dim = 1:numDim
  
  curBins = curbins(dim,curSel);
  curLabels = labels(curSel);
  idxpos = curLabels > 0;
  % pre-split up curBins for parfor
  curBinsPos = curBins(idxpos);
  curBinsNeg = curBins(~idxpos);
  
  Zpos = nnz(idxpos); Zneg = nnz(~idxpos);
  posCount = histc(curBinsPos,edges);
  posCount = posCount(1:end-1) / Zpos;
  negCount = histc(curBinsNeg,edges);
  negCount = negCount(1:end-1) / Zneg;
  
  dir = ones(1,numBins);
  posLeft = cumsum(posCount);
  posRight = sum(posCount)-posLeft;
  negLeft = cumsum(negCount);
  negRight = sum(negCount)-negLeft;
  binErr = (posLeft+negRight)/2;
  negErr = binErr>0.5;
  binErr(negErr) = 1-binErr(negErr);
  dir(negErr) = -1;

  for ndx = 1:numBins-1
    curVal = binVals(ndx,dim);
    pred_train = sum(sign(subX_train(:,dim)-curVal))/num_train;
    pred_mistakes = sum(sign(subX_cur(:,dim)-curVal))/num_mistakes;
    
    allSelErr(ndx,dim) = (pred_train-pred_mistakes)*dir(ndx)*sign(label_mistakes-1.5);
  end
  
  allBinErr(:,dim) = binErr;
  allDir(:,dim) = dir;
  
  
  
  err = binErr;
  [curBestErr(dim) binNo(dim)]= min(err(1:end-1));
  bestDir(dim) = dir(binNo(dim));
end



%% Compare to random switching of labels..

idx = false(size(data.windowdata.labelidx_new));

randbouts = randsample(find(bouts.label~=oldlabels(selex)),10);
for ndx = 1:10
  idx = idx | bouts.ndx(randbouts(ndx),:)';
end

newCls = {}; newScores = [];
data.windowdata.labelidx_new(idx) = 3 -data.windowdata.labelidx_new(idx);
hf = figure;
for ndx = 1:numTrain
  data.Train();
  newCls{ndx} = data.classifier;
  newScores(:,ndx) = myBoostClassify(data.windowdata.X,newCls{ndx});
  normdiff = (mean(newScores,2) - mean(oldScores,2))./stdScores;
  newmistakes = sum( (newScores.*repmat(modLabels,1,ndx))<0,2);
  sttnew = sprintf('%.3f ',newmistakes(sel2change)/ndx);
  sttold = sprintf('%.3f ',mistakes(sel2change)/numTrain);
  figure(hf); 
  subplot(2,2,1:2),hist(normdiff);
  title( {sprintf('Random!! Iter:%d dist:%.3f',ndx,normdiff(selex)),...
        sprintf('new mistakes:%s ',sttnew),...
        sprintf('oldmistakes:%s',sttold)});
  subplot(2,2,3); hist(normdiff(idx2rem))
  subplot(2,2,4); hist([oldScores(selex,1:ndx);newScores(selex,1:ndx)]');
end

data.windowdata.labelidx_new(idx) = 3 -data.windowdata.labelidx_new(idx);

%% Compute the distance to all frames for a fly.

dstruct = struct('dist',[],'exp',[],'flies',[]);
blockSize = 5000;
 
for expi = 1
  for flies = 3
    first_frame = data.firstframes_per_exp{expi}(flies);
    end_frame = data.GetTrxEndFrame(expi,flies);
    X = JLabelData.ComputeWindowDataChunkStatic(data.curperframefns,data.allperframefns,...
      data.GetPerframeFiles(expi),flies,data.windowfeaturescellparams,1,end_frame-first_frame+1);
  end
end
%% Get scores
oldScores = myBoostClassify(X,oldCl);
newScores = myBoostClassify(X,newCl);

%%
oldScores = []; newScores = [];
for ndx = 1:numTrain
oldScores(:,ndx) = myBoostClassify(X,oldCls{ndx});
newScores(:,ndx) = myBoostClassify(X,newCls{ndx});
end
%% fast bag dist mat.

distMat = zeros(size(X,1),numel(bmodel));
for ndx = 1:numel(bmodel)
  curWk = bmodel(ndx);
  dd = X(:,curWk.dim)*curWk.dir;
  tt = curWk.tr*curWk.dir;
  distMat(:,ndx) = sign( (dd>tt)-0.5)*curWk.alpha;
end

%% Distance to the bouts. fast bag..

bidx = bouts.ndx(randBout,:)>0;
boutMat = zeros(nnz(bidx),numel(bmodel));
for ndx = 1:numel(bmodel)
  curWk = bmodel(ndx);
  dd = data.windowdata.X(bidx,curWk.dim)*curWk.dir;
  tt = curWk.tr*curWk.dir;
  boutMat(:,ndx) = sign( (dd>tt)-0.5)*curWk.alpha;
end


%% Mean distance for all examples to the choosen bouts.

nframes = nnz(sel2change);
mean_dist = zeros(size(X,1),1);
for bndx = 1:numel(oldBagModels)
  distMat = BagDistMat(X,oldBagModels(bndx));
  boutMat = BagDistMat(data.windowdata.X(sel2change,:),oldBagModels(bndx));
  for ndx = 1:nframes
    cur_dist = sum(abs(distMat - repmat(boutMat(ndx,:),size(distMat,1),1) ),2);
    mean_dist = mean_dist+ cur_dist;
  end
end

%%

[sort_dist ord ]  = sort(mean_dist,'ascend');
[newScores(ord(1:40))' ;oldScores(ord(1:40))';sort_dist(1:40)'/10000 ]

%%
[sort_dist ord ]  = sort(mean_dist,'ascend');

bscores = zeros(size(X,1),numel(oldBagModels));
for ndx = 1:numel(oldBagModels)
  bscores(:,ndx) = myBoostClassify(X,oldBagModels{ndx});
end

tot = 0;
for ndx = 1:numel(oldBagModels)
  tot = tot+sum(abs([oldBagModels{ndx}(:).alpha]));
  
end

stdScores = std(bscores,[],2);
normdiff = (sum(newScores,2) - sum(oldScores,2))./stdScores/numTrain;
[normdiff(ord(1:40))'; sort_dist(1:40)'/tot/nframes/2 ]

%%

figure; plot(cumsum(normdiff(ord)<-1))

%% Intra training distances.

sel1 = randsample(find(oldlabels>0),70);
sel2 = randsample(find(oldlabels>0),70);

dist = zeros(numel(sel2change),70);
for bndx = 1:numel(oldBagModels)
  distMat1 = BagDistMat(data.windowdata.X(sel1,:),oldBagModels(bndx));
  distMat2 = BagDistMat(data.windowdata.X(sel2change,:),oldBagModels(bndx));
  for ndx = 1:numel(sel1)
    dist(:,ndx) = dist(:,ndx) + sum(abs(distMat2 - repmat(distMat1(ndx,:),size(distMat2,1),1) ),2);
  end
end

%% Mean distance for all the training set to the top 40.
nums = 50;
traindist = zeros(nums,size(data.windowdata.X,1));
for bndx = 1:numel(oldBagModels)
  trainMat = BagDistMat(data.windowdata.X,oldBagModels(bndx));
  testMat = BagDistMat(X(ord(1:nums),:),oldBagModels(bndx));
  for ndx = 1:nums
    traindist(ndx,:) = traindist(ndx,:)+sum(abs(trainMat - repmat(testMat(ndx,:),size(trainMat,1),1) ),2)';
  end
end


%% Features used for bagging and the classifiers

bfeatures = [];bwts = [];
for ndx = 1:numel(oldBagModels) 
  bfeatures = [bfeatures [oldBagModels{ndx}(:).dim]];
  bwts = [bwts [oldBagModels{ndx}(:).alpha]];

end

oldfeatures = [];oldwts = [];
for ndx = 1:numel(oldCls)
  oldfeatures = [oldfeatures [oldCls{ndx}(:).dim]];
  oldwts = [oldwts [oldCls{ndx}(:).alpha]];
end

newfeatures = []; newwts = [];
for ndx = 1:numel(newCls)
  newfeatures = [newfeatures [newCls{ndx}(:).dim]];
  newwts = [newwts [newCls{ndx}(:).alpha]];
end

[ubfeat idx1 idx2] = unique(bfeatures);
ubwt = [];
for ndx = 1:numel(ubfeat)
  ubwt(ndx) = sum(bwts(idx2 == ndx));
end
[sortbwt,sort_ord] = sort(ubwt);
ubfeat = ubfeat(sort_ord);

[uoldfeat idx1 idx2] = unique(oldfeatures);
uoldwt = [];
for ndx = 1:numel(uoldfeat)
  uoldwt(ndx) = sum(oldwts(idx2 == ndx));
end
[sortoldwt,sort_ord] = sort(uoldwt);
uoldfeat = uoldfeat(sort_ord);

[unewfeat idx1 idx2] = unique(newfeatures);
unewwt = [];
for ndx = 1:numel(unewfeat)
  unewwt(ndx) = sum(newwts(idx2 == ndx));
end

[sortnewwt,sort_ord] = sort(unewwt);
unewfeat = unewfeat(sort_ord);

%%
[com comnew comold] = intersect(unewfeat,uoldfeat);
sum(sortnewwt(comnew))
sum(sortoldwt(comold))
tt = false(size(unewfeat));
tt(comnew) = true;
sum(sortnewwt(~tt))

[com comnew combag] = intersect(unewfeat,ubfeat);
tt = false(size(unewfeat));
tt(comnew) = true;
sum(sortnewwt(~tt))

%% 

for ndx = 1:50
  [vv oo] = sort(traindist(ndx,:));
  jj = find( ismember(oo,find(bidx)),3);
  fprintf('%.3f %d %d %d\n',normdiff(ord(ndx)),jj(1),jj(2),jj(3));
end

%% Choose a random bout from test and find the closest in training.

len = 2;

meanoldscores = mean(oldScores,2);

selB = randsample(find(abs(meanoldscores)<0.5),1);
mean_dist_sel = zeros(size(data.windowdata.X,1),1);
for bndx = 1:numel(oldBagModels)
  trainMat = BagDistMat(data.windowdata.X,oldBagModels(bndx));
  selMat = BagDistMat(X(selB+[-len:len],:),oldBagModels(bndx));
  for ndx = 1:2*len+1
    mean_dist_sel = mean_dist_sel+sum(abs(trainMat - repmat(selMat(ndx,:),size(trainMat,1),1) ),2);
  end
end

[dist_sel seldistord] = sort(mean_dist_sel);


%% Change the label for the closest and train.

sel2change = [seldistord(1)+[-3:3] seldistord(3)+[-3:3]] ;
numTrain = 20; newSelCls = {}; oldCls = {};
numChange = 6;

idx = false(size(data.windowdata.labelidx_new));

idx(sel2change) = true;

data.windowdata.labelidx_new(idx) = 3 -data.windowdata.labelidx_new(idx);
for ndx = 1:numTrain
  data.Train();
  newSelCls{ndx} = data.classifier;
end
data.windowdata.labelidx_new(idx) = 3 -data.windowdata.labelidx_new(idx);

%%

newselscores = [];
for ndx = 1:numel(newSelCls),
  newselscores(:,ndx) = myBoostClassify(X,newSelCls{ndx});
end

scoresdiff = (sum(newselscores,2)-sum(oldScores(:,1:numTrain),2))./stdScores/numel(newSelCls);
scoresdiff(selB+[-len:len])
figure; hist(scoresdiff);

