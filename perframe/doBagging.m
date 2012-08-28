function [bagModels trainDistMat] = doBagging(data,labels,obj,binVals,bins,params)

boostIterations = 25;
% Learn classifier with all the data.
numRepeat = 20;

numEx = size(data,1);
posEx = labels == 1;
negEx = ~posEx;
numPos = sum(posEx);
numNeg = sum(negEx);

if numPos<1 || numNeg<1,
  allDataModel = struct('dim',1,'error',0.5,'dir',1,'tr',0,'alpha',0);
  for ndx = 1:(6*numRepeat)
    bagModels{ndx} = allDataModel;
  end
  trainDistMat = zeros(numEx,6*numRepeat*boostIterations);
  return;
end

modLabels = sign( (labels==1)-0.5);

count = 1;
outOfBag = ones(numEx,numRepeat*6);
bagModels = {};

for numIter = 1:numRepeat
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
      wt = getWeights(curTrainLabels);
      tt = tic;
      [scores curModel] = loglossboostLearnRandomFeatures(data(curTrain{ndx},:),curTrainLabels,...
        boostIterations,wt,binVals,bins(:,curTrain{ndx}),params);
      tScores = myBoostClassify(data(curTrain{3-ndx},:),curModel);
      outOfBag(curTrain{ndx},count)=0;
      bagModels{count} = curModel;
      count = count+1;
      etime = toc(tt);
      if ~isempty(obj)
        obj.SetStatus('%d%% training done. Time Remaining:%ds ',...
          round( count/(numRepeat*6+1)*100), round((numRepeat*6-count+1)*etime)); 
        drawnow();
      else
        fprintf('%d%% training done. Time Remaining:%ds \n',...
          round( count/(numRepeat*6+1)*100), round((numRepeat*6-count+1)*etime)); 
      end
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
    dd = data(:,curWk.dim);
    tr = curWk.tr;
    if curWk.dir>0,
      tt = ((dd>tr)-0.5)*2;
    else
      tt = ((dd<=tr)-0.5)*2;
    end
    trainDistMat(:,count) = tt*curWk.alpha;
    trainDistMat(~validEx,count)=nan;
    count = count+1;
  end
end
