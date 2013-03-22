function [newModel,scores,stats] = boostingUpdate(data,labels,oldModel,binVals,bins,params)

% Learn classifier with all the data.

posEx = labels == 1;
negEx = ~posEx;
numPos = sum(posEx);
numNeg = sum(negEx);

if numPos<1 || numNeg<1,
  newModel = struct('dim',1,'error',0.5,'dir',1,'tr',0,'alpha',0);
  scores = zeros(1,numel(labels));
  return;
end

wt = ones(length(posEx),1);
posWt = numNeg/(numPos+numNeg);
negWt = numPos/(numPos+numNeg);
wt(posEx) = posWt;
wt(negEx) = negWt;
wt = wt./sum(wt);

modLabels = sign( (labels==1)-0.5);

updatedModel = updateOldWeights(data,modLabels,oldModel,wt);
[newModel,scores] = addNewRules(data,modLabels,updatedModel,binVals,bins,wt,params);

if nargout >= 3,
  % compute some statistics of how well training worked
  stats = ComputeBoostingStats(scores,labels);
end



function updatedModel = updateOldWeights(data,labels,oldModel,exWt)

wt = exWt;
updatedModel = oldModel;
scores = zeros(size(data,1),1);
for ndx = 1:numel(oldModel)
  curWkRule = oldModel(ndx);
  tr = curWkRule.tr;
  dir = curWkRule.dir;
  dim = curWkRule.dim;
  if dir>0,
    tt = ((data(:,dim)> tr)-0.5)*2;
  else
    tt = ((data(:,dim)<= tr)-0.5)*2;
  end
  curError = sum( (tt.*labels).*wt);
  if(curError<0); 
    curError = 0; 
  end
  updatedModel(ndx).error = 0.5-curError/2;
  updatedModel(ndx).alpha = 1-2*updatedModel(ndx).error;
  scores = scores + myBoostClassify(data,updatedModel(ndx));
  tt = scores.*labels;
  wt = exWt./(1+exp(tt));
  wt = wt./sum(wt);
end


function [newModel,scores] = addNewRules(data,labels,updatedModel,binVals,bins,exWt,params)

extraIters = params.iter_updates;
newModel = updatedModel;
scores = myBoostClassify(data,updatedModel);
tt = scores.*labels;
wt = exWt./(1+exp(tt));
wt = wt./sum(wt);
for itt = 1:extraIters
  wkRule = findWeakRuleSamples(data,labels,wt,binVals,bins,params);

  tr = wkRule.tr;
  dir = wkRule.dir;
  dim = wkRule.dim;
  if dir>0,
    tt = ((data(:,dim)> tr)-0.5)*2;
  else
    tt = ((data(:,dim)<= tr)-0.5)*2;
  end
  curError = sum( (tt.*labels).*wt);
  if(curError<0); 
    curError = 0;
  end
  
  wkRule.error = 0.5-curError/2;
  wkRule.alpha = 1-2*wkRule.error;

  newModel(end+1) = wkRule;
  scores = scores + myBoostClassify(data,newModel(end));
  tt = scores.*labels;
  wt = exWt./(1+exp(tt));
  wt = wt./sum(wt);
end
