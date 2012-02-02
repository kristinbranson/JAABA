% mu = weighted_mean(x,w)
% x is n x d, w is n x 1
function mu = weighted_mean(x,w)

%[n,d] = size(x);

z = sum(w);
mu = sum(bsxfun(@times,x,w(:)),1)/z;