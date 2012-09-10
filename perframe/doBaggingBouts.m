function [bagModels trainDistMat] = doBaggingBouts(data,labels,obj,binVals,params,bouts)

boostIterations = 25;
% Learn classifier with all the data.
numRepeat = 200;

posBouts = bouts.label == 1;
negBouts = ~posBouts;

numPosBouts = nnz(posBouts);
numNegBouts = nnz(negBouts);

modLabels = sign( (labels==1)-0.5);
trainNdx = labels~=0;

if numPosBouts<4 || numNegBouts<4,

  allDataModel = struct('dim',1,'error',0.5,'dir',1,'tr',0,'alpha',0);
  for ndx = 1:(6*numRepeat)
    bagModels{ndx} = allDataModel;
  end
  trainDistMat = zeros(nnz(trainNdx),6*numRepeat*boostIterations);
  
  if numPosBouts<4
    msg = 'Too few bouts of behavior to do bagging';
  end
  if numNegBouts <4
    msg = 'Too few bouts of not behavior to do bagging';
  end
  
  return;
end

bins = findThresholdBins(data(trainNdx,:),binVals);

count = 1;
outOfBag = ones(nnz(trainNdx),numRepeat*6);
bagModels = {};

for numIter = 1:numRepeat
  
  sel = rand(1,numel(bouts.label))>0.5;
  
  curbouts = find( (posBouts|negBouts) &sel);
  curTrain = false(1,numel(bouts.ndx(1,:)));
  for ndx = curbouts(:)'
    curTrain = curTrain | bouts.ndx(ndx,:);
  end

  curTrainLabels = modLabels(curTrain);
  wt = getWeights(curTrainLabels);
  tt = tic;
  [~, curModel] = loglossboostLearnRandomFeatures(data(curTrain,:),curTrainLabels,...
    boostIterations,wt,binVals,bins(:,curTrain(trainNdx)),params);
  
  outOfBag(curTrain(trainNdx),count)=0;
  bagModels{count} = curModel;
  count = count+1;
  etime = toc(tt);
  if ~isempty(obj)
    obj.SetStatus('%d%% training done. Time Remaining:%ds ',...
      round( count/(numRepeat)*100), round((numRepeat-count+1)*etime));
    drawnow();
  else
    fprintf('%d%% training done. Time Remaining:%ds \n',...
      round( count/(numRepeat)*100), round((numRepeat-count+1)*etime));
  end
  
end


trainDistMat = zeros(nnz(trainNdx),length(bagModels)*length(bagModels{1}));
count = 1;
for mno = 1:length(bagModels)
  curModel = bagModels{mno};
  validEx = outOfBag(:,mno);
  for j = 1:length(curModel)
    curWk = curModel(j);
    dd = data(trainNdx,curWk.dim)*curWk.dir;
    tt = curWk.tr*curWk.dir;
    trainDistMat(:,count) = sign( (dd>tt)-0.5)*curWk.alpha;
    trainDistMat(~validEx,count)=nan;
    count = count+1;
  end
end
