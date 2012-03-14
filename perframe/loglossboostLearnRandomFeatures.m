function [scores model binVals bins] = loglossboostLearnRandomFeatures(data,labels,numIters,initWt,binVals,bins,params,obj)

numEx = size(data,1);
wt = initWt;
model = struct('dim',{},'error',{},'dir',{},'tr',{},'alpha',{});
scores = zeros(numEx,1);

clk = tic;
etimehist = [];
etimehistlength = 5;
nittoutput = 3;
origbins = bins; origdata = data; origbinVals = binVals;
for itt = 1:numIters
  if mod(itt,10)==1
    selFeatures = rand(1,size(origbins,1))>0.3;
    feature_map = find(~selFeatures);
    bins = origbins;
    data = origdata;
    binVals = origbinVals;
    bins(selFeatures,:) = [];
    data(:,selFeatures) = [];
    binVals(:,selFeatures) = [];
  end
  wkRule = findWeakRuleSamples(data,labels,wt,binVals,bins,params);
  wkRule.dim = feature_map(wkRule.dim);
  tr = wkRule.tr;
  dir = wkRule.dir;
  dim = wkRule.dim;
  if dir>0,
    tt = ((origdata(:,dim)> tr)-0.5)*2;
  else
    tt = ((origdata(:,dim)<= tr)-0.5)*2;
  end
  curError = sum( (tt.*labels).*wt);
  if(curError<0); 
    curError = - curError;
    wkRule.dir = -wkRule.dir;
  end
  
  wkRule.error = 0.5-curError/2;
  wkRule.alpha = 1-2*wkRule.error;
  model(itt) = wkRule;
  scores = myBoostClassify(origdata,model(1:itt));
  tt = scores.*labels;
  wt = initWt./(1+exp(tt));
  wt = wt./sum(wt);
  
  if( nargin>7 && mod(itt,nittoutput)==0)
    etime = toc(clk);
    etimehist(end+1) = etime;
    if numel(etimehist) > etimehistlength,
      etimehist = etimehist(end-etimehistlength+1:end);
    end
    obj.SetStatus('%d%% training done. Time Remaining:%ds ',...
      round(itt/numIters*100),round((numIters-itt)/nittoutput*mean(etimehist)));
    drawnow();
    clk = tic;
  end
  
end
dump = toc(clk);
end

