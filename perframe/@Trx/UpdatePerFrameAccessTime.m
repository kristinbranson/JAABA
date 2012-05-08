function UpdatePerFrameAccessTime(obj,fn)

% update access time for this property
j = find(strcmp(obj.perframehistory(:,1),fn),1);
if isempty(j),
  j = size(obj.perframehistory,1)+1;
end
obj.perframehistory{j,1} = fn;
obj.perframehistory{j,2} = now;
