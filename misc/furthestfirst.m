% [mu,idx] = furthestfirst(x,k,mu0)
% inputs: 
% x: N x D matrix of data to be clustered. N is the number of data
% points, D is the number of dimensions
% k: number of centers
% mu0: [optional] first center
% if mu0 is not input, a point is chosen at random
function [mu,idx,mu_idx] = furthestfirst(x,k,varargin)

[n,d] = size(x);

[start,distance,weights] = myparse(varargin,'Start','mean','Distance','sqEuclidean','weights',ones(n,1));

%isweight = ~all(weights==weights(1));

% output centers
mu = zeros(k,d);
mu_idx = zeros(k,1);

% if the first center was not input, then choose a random point
if isnumeric(start)
  mu(1,:) = start;
elseif strcmpi(start,'mean'),
  mu(1,:) = sum( x.*repmat(weights,[1,d]), 1 ) / sum(weights);
  %mu(1,:) = mean(x,1);
  mu_idx(1) = nan;
else
  i = randsample(n,1,true,weights);
  %i = ceil(rand(1)*n);
  mu(1,:) = x(i,:);
  mu_idx(1) = i;
end;

Dall = inf(k,n);
for i = 2:k,

  % compute the minimum distance from all points to the centers
  Dall(i-1,:) = distfun(x,mu(i-1,:),distance);
  D = min(Dall(1:i-1,:),[],1);

  % choose the point furthest from the centers as the next center
  j = argmax(D);
  mu(i,:) = x(j,:);
  mu_idx(i) = j;

end;

% set the labels for each point
Dall(k,:) = distfun(x,mu(k,:),distance);
[D,idx] = min(Dall,[],1);
idx = idx(:);

function D = distfun(X, C, dist)
%DISTFUN Calculate point to cluster centroid distances.
[n,p] = size(X);
D = zeros(n,size(C,1));
nclusts = size(C,1);

switch lower(dist)
case 'sqeuclidean'
    for i = 1:nclusts
        D(:,i) = sum((X - C(repmat(i,n,1),:)).^2, 2);
    end
case 'cityblock'
    for i = 1:nclusts
        D(:,i) = sum(abs(X - C(repmat(i,n,1),:)), 2);
    end
case {'cosine','correlation'}
    % The points are normalized, centroids are not, so normalize
    % them
    normC = sqrt(sum(C.^2, 2));
    if any(normC < eps(class(normC))) % small relative to
                                      % unit-length data poin\
ts
        error('furthestfirst:ZeroCentroid', ...
              'Zero cluster centroid created.');
    end
    % This can be done without a loop, but the loop saves memory
    % allocations
    for i = 1:nclusts
        D(:,i) = 1 - (X * C(i,:)') ./ normC(i);
    end
case 'hamming'
    for i = 1:nclusts
        D(:,i) = sum(abs(X - C(repmat(i,n,1),:)), 2) / p;
    end
end

