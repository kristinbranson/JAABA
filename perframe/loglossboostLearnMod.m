function [scores model binVals bins] = loglossboostLearnMod(data,labels,numIters,initWt,obj,binVals,bins)

numEx = size(data,1);
if(nargin<4)
  wt = ones(numEx,1)/numEx;
else
  wt = initWt;
end
model = struct('dim',{},'error',{},'dir',{},'tr',{},'alpha',{});
scores = zeros(numEx,1);

for itt = 1:numIters
  wkRule = findWeakRuleSamples(data,labels,wt,binVals,bins);
  
  tr = wkRule.tr;
  dir = wkRule.dir;
  dim = wkRule.dim;
  tt = (((data(:,dim)*dir)> (dir*tr))-0.5)*2;
  curError = sum( (tt.*labels).*wt);
  if(curError<0); 
    curError = - curError;
    wkRule.dir = -wkRule.dir;
  end
  
  wkRule.error = 0.5-curError/2;
  wkRule.alpha = 1-2*wkRule.error;
  model(itt) = wkRule;
  scores = myBoostClassify(data,model(1:itt));
  tt = scores.*labels;
  wt = initWt./(1+exp(tt));
  wt = wt./sum(wt);
  if(mod(itt,10)==0)  
    obj.SetStatus('%d%% training done.. ',round(itt/numIters*100)); 
    drawnow();
  end
end
end

