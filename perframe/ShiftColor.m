
classdef ShiftColor < handle
  
  properties (Access = public)
    
  end
  
  methods (Access = public, Static = true)
    
    function newColor = shiftColorFwd(oldColor)
      oldSize = size(oldColor);
      oldColor = reshape(oldColor,[1 1 3]);
      hh = rgb2hsv(oldColor);
      hh(1) = mod(hh(1)+0.085,1);
      newColor = hsv2rgb(hh);
      newColor = reshape(newColor,oldSize);
    end
    
    function newColor = shiftColorBkwd(oldColor)
      oldSize = size(oldColor);
      oldColor = reshape(oldColor,[1 1 3]);
      hh = rgb2hsv(oldColor);
      hh(1) = mod(hh(1)-0.065,1);
      newColor = hsv2rgb(hh);
      newColor = reshape(newColor,oldSize);
      
    end

    function newColor = increaseIntensity(oldColor)
      oldSize = size(oldColor);
      oldColor = reshape(oldColor,[1 1 3]);
      hh = rgb2hsv(oldColor);
      hh(3) = min(hh(3)+0.2,1);
      newColor = hsv2rgb(hh);
      newColor = reshape(newColor,oldSize);
      
    end

  end
end