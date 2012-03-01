function scores= crossValidateBout(data,labels,bouts,obj,binVals,bins,params)

k = params.CVfolds;

% Learn classifier with all the data.

posBouts = bouts.label == 1;
negBouts = ~posBouts;

numPosBouts = nnz(posBouts);
numNegBouts = nnz(negBouts);

if numPosBouts<k || numNegBouts<k,
  warndlg('Too few bouts to do cross validation');
  scores = zeros(1,size(data,1));
  scores(labels==0) = [];
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
scores = zeros(1,size(data,1));

for bno = 1:k
  curPosTest = posCum >= posBlocks(bno) & ...
              posCum < posBlocks(bno+1) & ...
              posBouts;
  
  curNegTest = negCum >= negBlocks(bno) & ...
              negCum < negBlocks(bno+1) & ...
              negBouts;
            
  curPosTrain = find(~curPosTest & posBouts);
  curNegTrain = find(~curNegTest & negBouts);
  
  curTrainNdx = false(1,size(data,1));
  for posNdx = curPosTrain
    curTrainNdx = curTrainNdx | bouts.ndx(posNdx,:);
  end
  for negNdx = curNegTrain
    curTrainNdx = curTrainNdx | bouts.ndx(negNdx,:);
  end
  
  curTestNdx = ~curTrainNdx & labels' ~=0;
  
  curTrainLabels = modLabels(curTrainNdx);
  
  wt = getWeights(curTrainLabels);  
  tt = tic;
  [~,curModel] = loglossboostLearnMod(data(curTrainNdx,:),curTrainLabels,...
    params.iter,wt,binVals,bins(:,curTrainNdx),params);
  tScores = myBoostClassify(data(curTestNdx,:),curModel);
  scores(curTestNdx) = tScores;

  etime = toc(tt);
  obj.SetStatus('%d%% cross validation done. Time Remaining:%ds ',...
    round( bno/k*100), round((k-bno)*etime));
  drawnow();
 
end

scores(labels==0) = [];
