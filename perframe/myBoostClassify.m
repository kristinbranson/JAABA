function scores = myBoostClassify(data,model)

numEx = size(data,1);
scores = zeros(numEx,1);
for j = 1:length(model)
  curWk = model(j);
  dd = data(:,curWk.dim)*curWk.dir;
  tt = curWk.tr*curWk.dir;
  addScores = (2*(dd > tt)-1)*curWk.alpha;
  scores = scores + addScores;
end
