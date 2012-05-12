function ClearDataCache(obj,fn,ns)

if nargin < 3,
  ns = 1:obj.nexpdirs;
end

for i = ns,

  ndatacurr = 0;
  j = find(strcmp(fn,obj.fnscached{i}),1);
  if ~isempty(j),
    %fprintf('Clearing %s from cache at idx %d for %d\n',fn,j,i);
    for k = 1:numel(obj.datacached{i}{j}),
      ndatacurr = ndatacurr + numel(obj.datacached{i}{j}{k});
    end
    obj.datacached{i}(j) = [];
    obj.fnscached{i}(j) = [];
    obj.nfnscached(i) = obj.nfnscached(i)-1;
    obj.ndatacachedperexp(i) = obj.ndatacachedperexp(i) - ndatacurr;
    obj.ndatacached = obj.ndatacached - ndatacurr;
    %fprintf('Done clearing %s from cache at idx %d for %d\n',fn,j,i);
  end
  
end