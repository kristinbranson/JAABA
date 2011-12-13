function [scores model binVals bins] = loglossboostLearnMod(data,labels,numIters,initWt,binVals,bins,obj)

numEx = size(data,1);
wt = initWt;
model = struct('dim',{},'error',{},'dir',{},'tr',{},'alpha',{});
scores = zeros(numEx,1);

clk = tic;
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
  
  if( nargin>6 && mod(itt,10)==0)
    etime = toc(clk);
    obj.SetStatus('%d%% training done. Time Remaining:%ds ',...
      round(itt/numIters*100),round((numIters-itt)/10*etime));
    drawnow();
    clk = tic;
  end
  
end
dump = toc(clk);
end

