classdef ExperimentDetailHandler < handle
%ExperimentDetailHandler Provide experiment detail plots
%   If experiment-level views are desired when using the OlyDat Browser,
%   write a concrete ExperimentDetailHandler to implement those views. The
%   Browser will call the methods below during operation; override them as
%   appropriate.
%
%   See also OlyDat.ExperimentDetailHandler/open,
%   OlyDat.ExperimentDetailHandler/refresh,
%   OlyDat.ExperimentDetailHandler/close, OlyDat.Browser, OlyDat.Plot
   
    methods
                
        function open(obj,data) %#ok<MANU,INUSD>
            % open() is called when the user pushes the 'Experiment Detail'
            % button. The data argument is the struct (row) for the currently
            % selected experiment, or [] if there is currently no experiment
            % selected.
        end
        
        function refresh(obj,data) %#ok<INUSD,MANU>
            % refresh() is called when the current experiment selection
            % changes. The data argument is the struct (row) for the currently
            % selected experiment, or [] if there is currently no experiment
            % selected.
        end
        
        function close(obj) %#ok<MANU>
            % close() is called when the Browser is shutting down.
        end
        
    end   
    
end