function [allDataModel outScores] = boostingWrapper(data,labels,obj,binVals,bins)

boostIterations = 100;
% Learn classifier with all the data.
numRepeat = 2;

numEx = size(data,1);
posEx = labels == 1;
negEx = ~posEx;
numPos = sum(posEx);
numNeg = sum(negEx);

if numPos<1 || numNeg<1,
  allDataModel = struct('dim',1,'error',0.5,'dir',1,'tr',0,'alpha',0);
  outScores = zeros(1,numEx);
  return;
end

wt = ones(length(posEx),1);
posWt = numNeg/(numPos+numNeg);
negWt = numPos/(numPos+numNeg);
wt(posEx) = posWt;
wt(negEx) = negWt;
wt = wt./sum(wt);

modLabels = sign( (labels==1)-0.5);
tt = tic;
[outScores, allDataModel] = loglossboostLearnMod(data,modLabels,boostIterations,wt,binVals,bins);
etime = toc(tt);
obj.SetStatus('%d%% training done. Time Remaining:%ds ',...
  round(1/(numRepeat*6+1)*100),round(numRepeat*6*etime)); 
drawnow();

return;


