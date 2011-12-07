function res = MeanWindowCore(x,w)

% average
fil = ones(1,w);
% full: res(t+r) corresponds to frame t
res = imfilter(x,fil,'full',0);
% normalize
res(w:end-w+1) = res(w:end-w+1) / w;
% boundary conditions
res(1:w-1) = bsxfun(@rdivide,res(1:w-1),1:w-1);
res(end-w+2:end) = bsxfun(@rdivide,res(end-w+2:end),w-1:-1:1);
