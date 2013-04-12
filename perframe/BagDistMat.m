function distMat = BagDistMat(data,bagModels)

distMat = zeros(size(data,1),length(bagModels)*length(bagModels{1}),'single');
count = 1;
for mno = 1:length(bagModels)
  curModel = bagModels{mno};
  for j = 1:length(curModel)
    curWk = curModel(j);
    dd = data(:,curWk.dim)*curWk.dir;
    tt = curWk.tr*curWk.dir;
    distMat(:,count) = sign( (dd>tt)-0.5)*curWk.alpha;
    count = count+1;
  end
end
