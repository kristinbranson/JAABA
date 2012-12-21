
%% Set the classifier and find a bagging model.


clfile = '/groups/branson/home/robiea/Code_versioned/Jdetect_github/Jdetect/perframe/ChaseAR_newv9.mat';
cl = load(clfile);
configfile = cl.configfilename;
data = JLabelData(configfile,'openmovie',false);
data.SetClassifierFileName(clfile);
data.Train();

oldCl = data.classifier;


%%
islabeled = data.windowdata.labelidx_new ~= 0;
bins = findThresholdBins(data.windowdata.X(islabeled,:),data.windowdata.binVals);
bagmodel = fastBag(data.windowdata.X(islabeled,:),...
  data.windowdata.labelidx_new(islabeled),...
  data.windowdata.binVals,bins,data.classifier_params);

%%
clear bmodel

count = 1;
for ndx = 1:numel(bagmodel)
  numwk = numel(bagmodel{ndx}.dim);
  for wndx = 1:numwk
    bmodel(count).dim = bagmodel{ndx}.dim(wndx);
    bmodel(count).dir = bagmodel{ndx}.dir(wndx);
    bmodel(count).error = bagmodel{ndx}.error(wndx);
    bmodel(count).tr = bagmodel{ndx}.tr(wndx);
    bmodel(count).alpha = bagmodel{ndx}.alpha(wndx);
    count = count+1;
  end
end

%% Switch labels for a random bout

bouts = data.getLabeledBouts();
randBout = randsample(numel(bouts.label),1);
idx = bouts.ndx(randBout,:)>0;
data.windowdata.labelidx_new(idx) = 3 -data.windowdata.labelidx_new(idx);
data.Train();
newCl = data.classifier;
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
%%

oldScores = myBoostClassify(X,oldCl);
newScores = myBoostClassify(X,newCl);
distMat = zeros(size(X,1),numel(bmodel));
for ndx = 1:numel(bmodel)
  curWk = bmodel(ndx);
  dd = X(:,curWk.dim)*curWk.dir;
  tt = curWk.tr*curWk.dir;
  distMat(:,ndx) = sign( (dd>tt)-0.5)*curWk.alpha;
end

%% Distance to the bouts.

bidx = bouts.ndx(randBout,:)>0;
boutMat = zeros(nnz(bidx),numel(bmodel));
for ndx = 1:numel(bmodel)
  curWk = bmodel(ndx);
  dd = data.windowdata.X(bidx,curWk.dim)*curWk.dir;
  tt = curWk.tr*curWk.dir;
  boutMat(:,ndx) = sign( (dd>tt)-0.5)*curWk.alpha;
end


%% Mean distance.

nframes = nnz(bidx);
mean_dist = zeros(size(X,1),1);
for ndx = 1:nframes
  cur_dist = sum(abs(distMat - repmat(boutMat(ndx,:),size(distMat,1),1) ),2);
  mean_dist = mean_dist+ cur_dist;
end

%%

[sort_dist ord ]  = sort(mean_dist,'ascend');
[newScores(ord(1:40))' ;oldScores(ord(1:40))';sort_dist(1:40)'/10000 ]

