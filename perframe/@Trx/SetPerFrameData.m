function SetPerFrameData(obj,fn,x,varargin)

nidx = numel(varargin);
if nidx < 1,
  error('Index not specified');
end
if nidx > 2,
  error('At most two indexing parameters allowed');
end
if ischar(varargin{1}),
  n = find(strcmp(obj.expdirs,varargin{1}),1);
  if isempty(n),
    error('Unknown experiment %s',varargin{1});
  end
elseif isnumeric(varargin{1}),
  n = varargin{1};
  if numel(n) > 1 || n < 1 || round(n) ~= n || n > obj.nexpdirs,
    error('Illegal experiment index %s',mat2str(n));
  end
else
  error('Illegal experiment index');
end
if nidx == 1,
  flies = 1:obj.nfliespermovie(n);
else
  flies = varargin{2};
  if ~isnumeric(flies) || any(flies < 1) || any(flies > obj.nfliespermovie(n)) || ...
      any(flies ~= round(flies)),
    error('Illegal fly indices');
  end  
end

if iscell(x),
  ndataadd = sum(cellfun(@numel,x(1:numel(flies))));
else
  ndataadd = numel(x);
end
obj.FreeDataCache(ndataadd);

obj.nfnscached(n) = obj.nfnscached(n) + 1;
j = obj.nfnscached(n);
%fprintf('Incremented nfnscached(%d) to %d for %s\n',n,obj.nfnscached(n),fn);
obj.fnscached{n}{j} = fn;
obj.datacached{n}{j} = cell(1,numel(flies));

for flyidx = 1:numel(flies),
  fly = flies(flyidx);
  if iscell(x),
    xcurr = x{flyidx};
  else
    xcurr = x;
  end
%   % delete data from cache if necessary
%   ndataadd = numel(xcurr);
%   obj.FreeDataCache(ndataadd);

  % add to cache
  %fprintf('Adding %s to data cache for fly %d at idx %d\n',fn,fly,j);
  obj.datacached{n}{j}{fly} = xcurr;
  
end
  
% update cache size
obj.ndatacached = obj.ndatacached + ndataadd;
obj.ndatacachedperexp(n) = obj.ndatacachedperexp(n) + ndataadd;

% update access time for this property
j = find(strcmp(obj.perframehistory(:,1),fn),1);
if isempty(j),
  j = size(obj.perframehistory,1)+1;
end
obj.perframehistory{j,1} = fn;
obj.perframehistory{j,2} = now;
