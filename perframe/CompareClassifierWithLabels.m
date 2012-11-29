function CompareClassifierWithLabels(clfile,ncomparisons)
% Compares the classifier's predictions with a classifier trained with the current labels.

if nargin <2
  ncomparisons = 3;
end

Q = load(clfile);
configfile = Q.configfilename;

data = JLabelData(configfile,'openmovie',false);
data.SetClassifierFileName(clfile);

allOrig = getCurrentScores(data);


allNew = {};
for ndx = 1:ncomparisons
  data.Train();
  allNew{ndx} = getCurrentScores(data);
end

origdiff = [];
newdiff = [];
sn = data.windowdata.scoreNorm;
bins = linspace(-sn,sn,20);
for ndx = 1:ncomparisons
  t = histc(allOrig - allNew{ndx},[-inf bins inf]);
  t(1) = []; t(end) = [];
  origdiff(end+1,:) = t;
  
  ndx1 = randsample(setdiff(1:ncomparisons,ndx),1);
  t = histc(allNew{ndx1} - allNew{ndx},[-inf bins inf]);
  t(1) = []; t(end) = [];
  newdiff(end+1,:) = t;
end

figure; hold on;
for ndx = 1:ncomparisons
  plot(bins/sn,origdiff(ndx,:),'r');
  plot(bins/sn,newdiff(ndx,:),'b');
end
title('Red: Orig - new, Blue: new - new');

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
