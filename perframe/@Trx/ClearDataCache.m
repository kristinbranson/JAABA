function ClearDataCache(obj,fn,ns)

if nargin < 3,
  ns = 1:obj.nexpdirs;
end

for i = ns,

  ndatacurr = 0;
  if isfield(obj.datacached{i},fn),
    for k = 1:numel(obj.datacached{i}),
      ndatacurr = ndatacurr + numel(obj.datacached{i}(k).(fn));
    end
    obj.datacached{i} = rmfield(obj.datacached{i},fn);
    obj.ndatacachedperexp(i) = obj.ndatacachedperexp(i) - ndatacurr;
    obj.ndatacached = obj.ndatacached - ndatacurr;
  end
  
end