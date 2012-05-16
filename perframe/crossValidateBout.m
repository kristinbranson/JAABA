function [success,msg,scores,tlabels] = crossValidateBout(data,labels,bouts,obj,binVals,params,timed)

success = true; msg = '';
k = params.CVfolds;

% Learn classifier with all the data.
if nargin< 8,
  timed = false;
end

posBouts = bouts.label == 1;
negBouts = ~posBouts;

numPosBouts = nnz(posBouts);
numNegBouts = nnz(negBouts);

if numPosBouts<k || numNegBouts<k,
  scores = zeros(1,size(data,1));
  scores(labels==0) = [];
  tlabels = {};
  success = false;
  
  if numPosBouts<k
    msg = 'Too few bouts of behavior to do cross validation';
  end
  if numNegBouts <k
    msg = 'Too few bouts of not behavior to do cross validation';
  end
  
  return;
end

posBlocks = linspace(0,numPosBouts+1,k+1);
negBlocks = linspace(0,numNegBouts+1,k+1);
posCum = cumsum(posBouts);
negCum = cumsum(negBouts);
% Randomly permute the bouts.
posCum = posCum(randperm(numel(posCum)));
negCum = negCum(randperm(numel(negCum)));

modLabels = sign( (labels==1)-0.5);
tlabels = {};
if timed,
  tpoints(1) = max(bouts.timestamp);
  tlabels{1} = datestr(tpoints(1));
  tpoints(end+1) = addtodate(tpoints(1),-5,'minute');
  tlabels{end+1} = '-5m';
  tpoints(end+1) = addtodate(tpoints(1),-10,'minute');
  tlabels{end+1} = '-10m';
  tpoints(end+1) = addtodate(tpoints(1),-15,'minute');
  tlabels{end+1} = '-15m';
  tpoints(end+1) = addtodate(tpoints(1),-20,'minute');
  tlabels{end+1} = '-20m';
  tpoints(end+1) = addtodate(tpoints(1),-25,'minute');
  tlabels{end+1} = '-25m';
  tpoints(end+1) = addtodate(tpoints(1),-30,'minute');
  tlabels{end+1} = '-30m';
  tpoints(end+1) = addtodate(tpoints(1),-45,'minute');
  tlabels{end+1} = '-45m';
  tpoints(end+1) = addtodate(tpoints(1),-1,'hour');
  tlabels{end+1} = '-1h';
  tpoints(end+1) = addtodate(tpoints(1),-2,'hour');
  tlabels{end+1} = '-2h';
  tpoints(end+1) = addtodate(tpoints(1),-3,'hour');
  tlabels{end+1} = '-3h';
  tpoints(end+1) = addtodate(tpoints(1),-4,'hour');
  tlabels{end+1} = '-4h';
  tpoints(end+1) = addtodate(tpoints(1),-1,'day');
  tlabels{end+1} = '-1d';
  tpoints(end+1) = addtodate(tpoints(1),-2,'day');
  tlabels{end+1} = '-2d';
  tpoints(end+1) = addtodate(tpoints(1),-3,'day');
  tlabels{end+1} = '-3d';
  tpoints(end+1) = addtodate(tpoints(1),-7,'day');
  tlabels{end+1} = '-1w';
  tpoints(end+1) = addtodate(tpoints(1),-7*2,'day');
  tlabels{end+1} = '-2w';
  tpoints(end+1) = addtodate(tpoints(1),-1,'month');
  tlabels{end+1} = '-1mo';
  tpoints(end+1) = addtodate(tpoints(1),-2,'month');
  tlabels{end+1} = '-2mo';
  tpoints(end+1) = addtodate(tpoints(1),-6,'month');
  tlabels{end+1} = '-6mo';
  tpoints(end+1) = addtodate(tpoints(1),-1,'year');
  tlabels{end+1} = '-1y';
  tpoints(end+1) = addtodate(tpoints(1),-2,'year');
  tlabels{end+1} = '-2y';
  tpoints(end+1) = addtodate(tpoints(1),-10,'year');
  tlabels{end+1} = '-10y';
  
  tooOld = tpoints< min(bouts.timestamp) ;
  tpoints(tooOld) = [];
  tlabels(tooOld) = [];
else
  tpoints = now;
end


scores = zeros(numel(tpoints),size(data,1));
bins = findThresholdBins(data(labels~=0,:),binVals);

for bno = 1:k
  curPosTest = posCum >= posBlocks(bno) & ...
    posCum < posBlocks(bno+1) & ...
    posBouts;
  
  curNegTest = negCum >= negBlocks(bno) & ...
    negCum < negBlocks(bno+1) & ...
    negBouts;
  
  curTestNdx = false(1,size(data,1));
  for posNdx = find(curPosTest)
    curTestNdx = curTestNdx | bouts.ndx(posNdx,:);
  end
  
  for negNdx = find(curNegTest)
    curTestNdx = curTestNdx | bouts.ndx(negNdx,:);
  end
  
  
  for tndx = 1:numel(tpoints)
  
    curPosTrain = find(~curPosTest & posBouts & bouts.timestamp<=tpoints(tndx));
    curNegTrain = find(~curNegTest & negBouts & bouts.timestamp<=tpoints(tndx));

    curTrainNdx = false(1,size(data,1));
    for posNdx = curPosTrain
      curTrainNdx = curTrainNdx | bouts.ndx(posNdx,:);
    end
    for negNdx = curNegTrain
      curTrainNdx = curTrainNdx | bouts.ndx(negNdx,:);
    end

    curTrainLabels = modLabels(curTrainNdx);

    wt = getWeights(curTrainLabels);  
    tt = tic;
    [~,curModel] = loglossboostLearnRandomFeatures(data(curTrainNdx,:),curTrainLabels,...
      params.iter,wt,binVals,bins,params);
    tScores = myBoostClassify(data(curTestNdx,:),curModel);
    scores(tndx,curTestNdx) = tScores;
    etime = toc(tt);
    done = ((bno-1)*numel(tpoints) + tndx);
    obj.SetStatus('%d%% cross validation done. Time Remaining:%ds ',...
      round( done/(numel(tpoints)*k)*100), ...
      round( ((numel(tpoints)*k)-done)*etime));
    drawnow();

  end
  
 
end

scores(:,labels==0) = [];
