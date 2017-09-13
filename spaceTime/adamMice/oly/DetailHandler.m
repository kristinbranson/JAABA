classdef DetailHandler < OlyDat.ExperimentDetailHandler
  % Concrete ExperimentDetailHandler for Hantman MouseJAABA
  
  methods
    
    function open(obj,data) %#ok<INUSL>
      if ~isempty(data)
        jabs = data(1).jabs;
        jab = jabs{1};
        expfull = data(1).expfull;
        
        JLabelInterf.openJabExpFrame(jab,expfull,1);       
      end
    end    
   
  end
  
  
end