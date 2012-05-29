function FreeDataCache(obj,ndataadd)

while true,
  
  if obj.ndatacached + ndataadd <= obj.maxdatacached || obj.ndatacached == 0,
    break;
  end
  if isempty(obj.perframehistory),
    break;
  end
  [~,j] = min(cell2mat(obj.perframehistory(:,2)));
  fn = obj.perframehistory{j,1};
  ClearDataCache(obj,fn);
%   for i = 1:obj.nexpdirs,
%     if isfield(obj.datacached{i},fn),
%       ndatacurr = 0;
%       for k = 1:numel(obj.datacached{i}),
%         ndatacurr = ndatacurr + numel(obj.datacached{i}(k).(fn));
%       end
%       obj.datacached{i} = rmfield(obj.datacached{i},fn);
%       obj.ndatacachedperexp(i) = obj.ndatacachedperexp(i) - ndatacurr;
%       obj.ndatacached = obj.ndatacached - ndatacurr;
%     end
%   end
  
  obj.perframehistory(j,:) = [];

end