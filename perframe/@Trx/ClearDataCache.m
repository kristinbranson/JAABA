function ClearDataCache(obj,fn,ns)

if nargin < 3,
  ns = 1:obj.nexpdirs;
end

for i = ns,

  ndatacurr = 0;
  j = find(strcmp(fn,obj.fnscached{i}),1);
  if ~isempty(j),
    for k = 1:numel(obj.datacached{i}{j}),
      ndatacurr = ndatacurr + numel(obj.datacached{i}{j}{k});
    end
    obj.datacached{i}(j) = [];
    obj.fnscached{i}(j) = [];
    obj.ndatacachedperexp(i) = obj.ndatacachedperexp(i) - ndatacurr;
    obj.ndatacached = obj.ndatacached - ndatacurr;
  end
  
end