function [newModel,scores] = boostingUpdate(data,labels,oldModel,binVals,bins)

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

updatedModel = updateOldWeights(data,modLabels,oldModel,wt);
[newModel,scores] = addNewRules(data,modLabels,updatedModel,binVals,bins,wt);


function updatedModel = updateOldWeights(data,labels,oldModel,exWt)

wt = exWt;
updatedModel = oldModel;
scores = zeros(size(data,1),1);
for ndx = 1:numel(oldModel)
  curWkRule = oldModel(ndx);
  tr = curWkRule.tr;
  dir = curWkRule.dir;
  dim = curWkRule.dim;
  tt = (((data(:,dim)*dir)>(dir*tr))-0.5)*2;
  curError = sum( (tt.*labels).*wt);
  if(curError<0); 
    curError = - curError;
    updatedModel(ndx).dir = -updatedModel(ndx).dir;
  end
  updatedModel(ndx).error = 0.5-curError/2;
  updatedModel(ndx).alpha = 1-2*updatedModel(ndx).error;
  scores = scores + myBoostClassify(data,updatedModel(ndx));
  tt = scores.*labels;
  wt = exWt./(1+exp(tt));
  wt = wt./sum(wt);
end


function [newModel,scores] = addNewRules(data,labels,updatedModel,binVals,bins,exWt)

extraIters = 10;

newModel = updatedModel;
scores = myBoostClassify(data,updatedModel);
tt = scores.*labels;
wt = exWt./(1+exp(tt));
wt = wt./sum(wt);
for itt = 1:extraIters
  wkRule = findWeakRuleSamples(data,labels,wt,binVals,bins);

  tr = wkRule.tr;
  dir = wkRule.dir;
  dim = wkRule.dim;
  tt = (((data(:,dim)*dir)>(dir*tr))-0.5)*2;
  curError = sum( (tt.*labels).*wt);
  if(curError<0); 
    curError = - curError;
    wkRule.dir = -wkRule.dir;
  end
  
  wkRule.error = 0.5-curError/2;
  wkRule.alpha = 1-2*wkRule.error;

  newModel(end+1) = wkRule;
  scores = scores + myBoostClassify(data,newModel(end));
  tt = scores.*labels;
  wt = exWt./(1+exp(tt));
  wt = wt./sum(wt);
end
