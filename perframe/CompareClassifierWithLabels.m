function CompareClassifierWithLabels(clfile,ncomparisons,explist)
% Compares the classifier's predictions with a classifier trained with the current labels.

if nargin <2
  ncomparisons = 3;
end

if nargin<3,
  explist = {};
end

Q = load(clfile);
configfile = Q.configfilename;

data = JLabelData(configfile,'openmovie',false);

if isempty(explist),
  data.SetClassifierFileName(clfile);
  allOrig = getCurrentScores(data);

else
  data.SetClassifierFileNameWoExp(clfile);
  for ndx = 1:numel(explist)
    [success,msg] = data.AddExpDir(explist{ndx});
    if ~success,
      warning(msg);
      continue;
    end
  end

  scores = [];
  for endx = 1:data.nexps
    sfn = data.GetFile('scores',endx);
    Q = load(sfn);
    allScores = Q.allScores;
    
    for ndx = 1:numel(allScores.scores);
      ts = allScores.tStart(ndx);
      te = allScores.tEnd(ndx);
      scores = [scores allScores.scores{ndx}(ts:te)];
    end
  end
  allOrig = scores;
  
end




allNew = {};
for ndx = 1:ncomparisons
  data.Train();
  allNew{ndx} = getCurrentScores(data);
end

origdiff = [];
newdiff = [];
origdiffpos = [];
newdiffpos = [];
origdiffneg = [];
newdiffneg = [];
sn = data.windowdata.scoreNorm;
bins = linspace(-sn,sn,20);
origmat = zeros(2,2,ncomparisons);
newmat = zeros(2,2,ncomparisons);
for ndx = 1:ncomparisons
  t = histc(allOrig - allNew{ndx},[-inf bins inf]);
  t(1) = []; t(end) = [];
  origdiff(end+1,:) = t;
  
  origmat(1,1,ndx) = nnz(allOrig>0 & allNew{ndx}>0);
  origmat(1,2,ndx) = nnz(allOrig>0 & allNew{ndx}<0);
  origmat(2,1,ndx) = nnz(allOrig<0 & allNew{ndx}>0);
  origmat(2,2,ndx) = nnz(allOrig<0 & allNew{ndx}<0);
  
  popul = setdiff(1:ncomparisons,ndx);
  if isempty(popul), 
    ndx1 = 1;
  else
    ndx1 = randsample(setdiff(1:ncomparisons,ndx),1);
  end
  t = histc(allNew{ndx1} - allNew{ndx},[-inf bins inf]);
  t(1) = []; t(end) = [];
  newdiff(end+1,:) = t;

  newmat(1,1,ndx) = nnz(allNew{ndx1}>0 & allNew{ndx}>0);
  newmat(1,2,ndx) = nnz(allNew{ndx1}>0 & allNew{ndx}<0);
  newmat(2,1,ndx) = nnz(allNew{ndx1}<0 & allNew{ndx}>0);
  newmat(2,2,ndx) = nnz(allNew{ndx1}<0 & allNew{ndx}<0);


  pos = allOrig>0;
  t = histc(allOrig(pos) - allNew{ndx}(pos),[-inf bins inf]);
  t(1) = []; t(end) = [];
  origdiffpos(end+1,:) = t;
  
  t = histc(allNew{ndx1}(pos) - allNew{ndx}(pos),[-inf bins inf]);
  t(1) = []; t(end) = [];
  newdiffpos(end+1,:) = t;

  pos = allOrig<0;
  t = histc(allOrig(pos) - allNew{ndx}(pos),[-inf bins inf]);
  t(1) = []; t(end) = [];
  origdiffneg(end+1,:) = t;
  
  t = histc(allNew{ndx1}(pos) - allNew{ndx}(pos),[-inf bins inf]);
  t(1) = []; t(end) = [];
  newdiffneg(end+1,:) = t;


end

origmat = round(mean(origmat,3));
newmat = round(mean(newmat,3));


figure;
subplot(2,2,[1 2]); hold on;
for ndx = 1:ncomparisons
  plot(bins/sn,origdiff(ndx,:),'r');
  plot(bins/sn,newdiff(ndx,:),'b');
end

numPosOrig = nnz(allOrig>0);
numNegOrig = nnz(allOrig<0);
numPosNew = 0; numNegNew = 0;
for ndx = 1:numel(allNew)
  numPosNew = numPosNew + nnz(allNew{ndx}>0);
  numNegNew = numNegNew + nnz(allNew{ndx}<0);
  
end
numPosNew = numPosNew/numel(allNew);
numNegNew = numNegNew/numel(allNew);


title({'Red: Orig - new, Blue: new - new'});
msg = {...
 sprintf('Frames predicted by the orig, pos :%d neg:%d',numPosOrig,numNegOrig),...
 sprintf('Frames predicted by the new (avg), pos :%d neg:%d',round(numPosNew),round(numNegNew)),...
 sprintf('Orig Vs New'),...
 sprintf('           New Pos          New Neg     '),...
 sprintf('Orig Pos     %d          %d     ',origmat(1,1),origmat(1,2)),...
 sprintf('Orig Neg     %d          %d     ',origmat(2,1),origmat(2,2)),...
 sprintf('New Vs New'),...
 sprintf('           New Pos          New Neg     '),...
 sprintf('New Pos     %d          %d     ',newmat(1,1),newmat(1,2)),...
 sprintf('New Neg     %d          %d     ',newmat(2,1),newmat(2,2)),...
};
warndlg(msg);

subplot(2,2,3); hold on;
for ndx = 1:ncomparisons
  plot(bins/sn,origdiffpos(ndx,:),'r');
  plot(bins/sn,newdiffpos(ndx,:),'b');
end
title('On frames predicted positive by orig');

subplot(2,2,4); hold on;
for ndx = 1:ncomparisons
  plot(bins/sn,origdiffneg(ndx,:),'r');
  plot(bins/sn,newdiffneg(ndx,:),'b');
end
title('On frames predicted negative by orig');



function scores = getCurrentScores(data)
scores = [];
for expi = 1:data.nexps
  allScores = data.PredictWholeMovie(expi);
  for ndx = 1:numel(allScores.scores);
    ts = allScores.tStart(ndx);
    te = allScores.tEnd(ndx);
    scores = [scores allScores.scores{ndx}(ts:te)];
  end
end


