function [sorted_weights,feature_order,bins,scores] = ...
  sortWindowFeaturesByWeight(classifier,data)

alpha = [classifier.alpha];
dim = [classifier.dim];
nfeatures = size(data,2);
nbins = 11;

if numel(alpha) ~= numel(dim),
  error('Sanity check: alpha and dim should be the same size');
end

% total absolute weight of all rules using each feature
totalweight = zeros(1,nfeatures);

scores = zeros(nbins,nfeatures);
bins = prctile(data,linspace(0,100,nbins),1);

% loop through features
for i = 1:nfeatures,
  
  idx = dim == i;
  if ~any(idx),
    continue;
  end
  
  % add up absolute weights
  totalweight(i) = sum(abs(alpha(idx)));
  
  % weak rules using just this feature
  classifier_curr = classifier(idx);
  [classifier_curr.dim] = deal(1);
    
  % classify using just this feature
  scores(:,i) = myBoostClassify(bins(:,i),classifier_curr);
  
end
[sorted_weights,feature_order] = sort(totalweight,2,'descend');
bins = bins(:,feature_order);
scores = scores(:,feature_order);
