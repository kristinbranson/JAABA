% ... = obj.GetPerFrameData(fn,flyidx)
% ... = obj.GetPerFrameData(fn,expidx,flyidx)
function [varargout] = GetPerFrameData(obj,fn,varargin)

nidx = numel(varargin);
if nidx < 1,
  error('Index not specified');
end
if nidx > 2,
  error('At most two indexing parameters allowed');
end
if nidx == 1,
  flyidxes = varargin{1};
  [ns,flies] = obj.getExpFly(flyidxes);
else
  if ischar(varargin{1}),
    ns = find(strcmp(obj.expdirs,varargin{1}),1);
    if isempty(ns),
      error('Unknown experiment %s',varargin{1});
    end
  elseif iscell(varargin{1}),
    [didfind,ns] = ismember(varargin{1},obj.expdirs);
    if any(~didfind),
      error(['Unknown experiments ',sprintf('%s ',varargin{1}{~didfind})]);
    end
  elseif isnumeric(varargin{1}),
    ns = varargin{1};
    if any(ns < 1) || any(round(ns) ~= ns) || any(ns > obj.nexpdirs),
      error('Illegal experiment indices %s',mat2str(ns));
    end
  end
  
  flies = varargin{2};
  if ~isnumeric(flies) || any(flies < 1) || any(flies > obj.nfliespermovie(ns)) || ...
      any(flies ~= round(flies)),
    error('Illegal fly indices');
  end
  
  [ns,flies] = meshgrid(ns(:),flies(:));
  ns = ns(:)';
  flies = flies(:)';

end

res = cell(1,numel(ns));
for i = 1:numel(ns),
  
  n = ns(i);
  fly = flies(i);

  if ~isfield(obj.datacached{n},fn)
    x = obj.LoadPerFrameData(fn,n);
    if iscell(x),
      res{i} = x{fly};
    else
      res{i} = x;
    end
  else
    res{i} = obj.datacached{n}(fly).(fn);
  end

end

if nargout < numel(res),
  if numel(res) == 1,
    varargout{1:nargout} = res{1};
  else
    varargout{1:nargout} = res;
  end
else
  varargout = res;
end

% update access time for this property
obj.UpdatePerFrameAccessTime(fn);
