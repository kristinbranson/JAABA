function scores= crossValidate(data,labels,obj,binVals,bins,params)

k = params.CVfolds;

boostIterations = 100;
% Learn classifier with all the data.

numEx = size(data,1);
posEx = labels == 1;
posCum = cumsum(posEx)/nnz(posEx);
negEx = ~posEx;
negCum = cumsum(negEx)/nnz(negEx);
numPos = sum(posEx);
numNeg = sum(negEx);

scores = zeros(1,numEx);
if numPos<1 || numNeg<1,
  obj.setStatus('Training set has only one kind of labels');
  return;
end

modLabels = sign( (labels==1)-0.5);

% rr = 1:numEx;
% splitPt = randi(numEx);
% rr = [rr(splitPt:end) rr(1:splitPt-1)];
% bStarts = round(linspace(1,numEx+1,k+1));

for bno = 1:k
%   curBlock = rr(bStarts(bno):bStarts(bno+1)-1);
%   curTest = ismember(1:numEx,curBlock);
%   curTrain = ~curTest;
  
  curPos = (posCum <= (bno/k)) & (posCum > ( (bno-1)/k)) & posEx;
  curNeg = (negCum <= (bno/k)) & (negCum > ( (bno-1)/k)) & negEx;
  
  curTest = curPos | curNeg ;
  curTrain = ~curTest;
  
  curTrainLabels = modLabels(curTrain);
  
  wt = getWeights(curTrainLabels);  
  tt = tic;
  [~,curModel] = loglossboostLearnMod(data(curTrain,:),curTrainLabels,...
    boostIterations,wt,binVals,bins(:,curTrain),params);
  tScores = myBoostClassify(data(curTest,:),curModel);
  scores(curTest) = tScores;
  
  etime = toc(tt);
  obj.SetStatus('%d%% cross validation done. Time Remaining:%ds ',...
    round( bno/k*100), round((k-bno)*etime));
  drawnow();
end

