function [allDataModel bagModels trainDistMat] = boostingWrapper(data,labels)

boostIterations = 100;
% Learn classifier with all the data.

posEx = labels == 1;
negEx = ~posEx;
numPos = sum(posEx);
numNeg = sum(negEx);

wt = ones(length(posEx),1);
posWt = numNeg/(numPos+numNeg);
negWt = numPos/(numPos+numNeg);
wt(posEx) = posWt;
wt(negEx) = negWt;
wt = wt./sum(wt);

modLabels = sign( (labels==1)-0.5);

[~, allDataModel] = loglossboostLearnMod(data,modLabels,boostIterations,wt);


% Do bagging.

numRepeat = 2;
numEx = size(data,1);

count = 1;
outOfBag = ones(numEx,numRepeat*6);
bagModels = {};

for numIter = 1:numRepeat
  rr = randperm(numEx);
  for bno = 1:4
    block{bno} = rr(bno:4:end);
  end
  
  for bno = 1:3
    curTrain{1} = ismember(1:numEx,[block{1} block{1+bno}]);
    curTrain{2} = ~curTrain{1};
    
    for ndx =1:2
      curTrainLabels = modLabels(curTrain{ndx});
      numPos = sum(curTrainLabels>0); 
      numNeg = sum(curTrainLabels<0);
      wt = zeros(length(curTrainLabels),1);
      negWt = numPos/(numNeg+numPos);
      posWt = numNeg/(numPos+numNeg);
      wt(curTrainLabels==1)=posWt;
      wt(curTrainLabels~=1)=negWt;
      wt = wt./sum(wt);

      [scores curModel] = loglossboostLearnMod(data(curTrain{ndx},:),curTrainLabels,boostIterations,wt);
      outOfBag(curTrain{ndx},count)=0;
      bagModels{count} = curModel;
      count = count+1;
    end    
  end
  
end

% Distance matrix based on bagging.

trainDistMat = zeros(numEx,length(bagModels)*length(bagModels{1}));
count = 1;
for mno = 1:length(bagModels)
  curModel = bagModels{mno};
  validEx = outOfBag(:,mno);
  for j = 1:length(curModel)
    curWk = curModel(j);
    dd = data(:,curWk.dim)*curWk.dir;
    tt = curWk.tr*curWk.dir;
    trainDistMat(:,count) = (dd>tt)*curWk.alpha;
    trainDistMat(~validEx,count)=nan;
    count = count+1;
  end
end
