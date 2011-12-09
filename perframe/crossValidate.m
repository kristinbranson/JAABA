function crossError = crossValidate(data,labels,obj,binVals,bins)

k = 3;

boostIterations = 100;
% Learn classifier with all the data.

numEx = size(data,1);
posEx = labels == 1;
negEx = ~posEx;
numPos = sum(posEx);
numNeg = sum(negEx);

if numPos<1 || numNeg<1,
  obj.setStatus('Training set has only one kind of labels');
  crossError = 1;
  return;
end

modLabels = sign( (labels==1)-0.5);

rr = 1:numEx;
splitPt = randi(numEx);
rr = [rr(splitPt:end) rr(1:splitPt-1)];
bStarts = round(linspace(1,numEx+1,k+1));
crossError = 0;

for bno = 1:k
  curBlock = rr(bStarts(bno):bStarts(bno+1)-1);
  curTest = ismember(1:numEx,curBlock);
  curTrain = ~curTest;
  
  curTrainLabels = modLabels(curTrain);
  curTestLabels = modLabels(curTest);
  
  wt = getWeights(curTrainLabels);  
  tt = tic;
  [~,curModel] = loglossboostLearnMod(data(curTrain,:),curTrainLabels,...
    boostIterations,wt,binVals,bins(:,curTrain));
  tScores = myBoostClassify(data(curTest,:),curModel);
  
  crossError = crossError + sum( (sign(tScores).*sign(curTestLabels))<0);
  
  etime = toc(tt);
  obj.SetStatus('%d%% cross validation done. Time Remaining:%ds ',...
    round( bno/k*100), round((k-bno)*etime));
  drawnow();
end

crossError = crossError/numEx;