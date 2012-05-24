function wt = getWeights(labels)

numPos = sum(labels>0);
numNeg = sum(labels<0);
wt = zeros(length(labels),1);
negWt = numPos/(numNeg+numPos);
posWt = numNeg/(numPos+numNeg);
wt(labels>0)=posWt;
wt(labels<0)=negWt;
wt = wt./sum(wt);
