function [allDataModel bagModels trainDistMat outScores] = boostingWrapper(data,labels,obj,binVals,bins)

boostIterations = 10;
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

[outScores, allDataModel] = loglossboostLearnMod(data,modLabels,boostIterations,wt,obj,binVals,bins);

% bagModels =[]; trainDistMat = []; 
% return;
% Do bagging.

numRepeat = 2;
numEx = size(data,1);

count = 1;
outOfBag = ones(numEx,numRepeat*6);
bagModels = {};

% outScores = zeros(numEx,1);
for numIter = 1:numRepeat
%   rr = randperm(numEx);
  rr = 1:numEx;
  splitPt = randi(numEx);
  rr = [rr(splitPt:end) rr(1:splitPt-1)];
  bStarts = round(linspace(1,numEx+1,5));
  for bno = 1:4
    block{bno} = rr(bStarts(bno):bStarts(bno+1)-1);
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

      [scores curModel] = loglossboostLearnMod(data(curTrain{ndx},:),curTrainLabels,boostIterations,wt,obj,binVals,bins);
      tScores = myBoostClassify(data(curTrain{3-ndx},:),curModel);
%       outScores(curTrain{3-ndx}) = outScores(curTrain{3-ndx}) + tScores;
      outOfBag(curTrain{ndx},count)=0;
      bagModels{count} = curModel;
      count = count+1;
    end    
  end
  
end

% outScores = outScores/6;
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
