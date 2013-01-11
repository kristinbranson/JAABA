function stringOut=upperFirstLowerRest(stringIn)

stringOut=lower(stringIn);
if ~isempty(stringOut)
  stringOut(1)=upper(stringIn(1));
end
  
end
