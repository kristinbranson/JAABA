function stats = ComputeBoostingStats(outScores,labels)
  
  stats = struct;
  
  % number of training examples
  stats.npos = numel(labels==1);
  stats.nneg = numel(labels==2);
  
  % number of errors
  stats.nfalsepos_train = nnz(outScores(labels==1)<=0);
  stats.nfalseneg_train = nnz(outScores(labels==2)>=0);

  stats.fracfalsepos_train = stats.nfalsepos_train / stats.npos;
  stats.fracfalseneg_train = stats.nfalseneg_train / stats.nneg;
  
  % worst score
  stats.worstscore_pos = min(outScores(labels==1));
  stats.worstscore_neg = -max(outScores(labels==2));
  stats.worstscore = min(stats.worstscore_pos,stats.worstscore_neg);

  % normalization for scores
  err = outScores;
  err(labels==2) = -err(labels==2);
  stats.norm = prctile(err,80);
  
  % normalized worst score
  stats.worstscorenorm_pos = stats.worstscore_pos / stats.norm;
  stats.worstscorenorm_neg = stats.worstscore_neg / stats.norm;
  stats.worstscorenorm = min(stats.worstscorenorm_pos,stats.worstscorenorm_neg);
  
  % other stats of scores
  stats.meanscore_pos = mean(outScores(labels==1));
  stats.meanscore_neg = -mean(outScores(labels==1));
  stats.meanscorenorm_pos = stats.meanscore_pos / stats.norm;
  stats.meanscorenorm_neg = stats.meanscore_neg / stats.norm;

  stats.stdscore_pos = std(outScores(labels==1),1);
  stats.stdscore_neg = std(outScores(labels==1),1);
  stats.stdscorenorm_pos = stats.stdscore_pos / stats.norm;
  stats.stdscorenorm_neg = stats.stdscore_neg / stats.norm;
