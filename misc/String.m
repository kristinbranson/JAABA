classdef String
  % String utils  

  % comma-separated list of items
  methods (Static)
    
    function c = commaSepList2CellStr(s)
      c = regexp(s,',','split');
      c = strtrim(c);
    end
    
    function s = cellstr2CommaSepList(c)
      if isempty(c)
        s = '';
      else
        s = [sprintf('%s,',c{1:end-1}) c{end}];
      end
    end
    
    % see civilizedStringFromCellArrayOfStrings
  end
  
end