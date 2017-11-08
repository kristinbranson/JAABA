classdef DataPuller < handle
%DataPuller Pull olympiad data from SAGE.
%   A DataPuller knows how, for a particular assay, to pull data from SAGE
%   based on a user-specified query, and merge/clean the results into a
%   single struct array. This is an assay-dependent operation; each assay
%   has differerent SAGE views and different requirements for merging and
%   cleaning.
%
%   DataPuller is meant to be used in conjunction with the
%   OlyDat.DataSelector. 
%
%   See also OlyDat.DataSelector, OlyDat.Browser.

    properties (Abstract)
        % if true, pull from FS by default.
        tfPullFSDefault; 
        
        % cellstr of options displayed in DataSelector pull-down and passed
        % to pullData function. can be empty.
        options; 
    end
    
    methods (Abstract)
        %pullData/pullDataFS Pull Data from SAGE
        %   pullData is a SAGE data-pulling function for use with
        %   OlyDat.DataSelector. queryObjs is a cell array of query objects
        %   created by the user using the DataSelector. datapath is the
        %   user-specified folder containing experimental directories.
        %   option is a string, selected by the user in the DataSelector;
        %   it is one of the options specified in the DataPuller .options
        %   property.
        %
        %   data should be a struct array that results from pulling data
        %   from the appropriate dataSets and cleaning and merging that
        %   data. This struct array should have one element per experiment
        %   for use with the OlyDat.Browser.
        data = pullData(obj,queryObjs,option)
        data = pullDataFS(obj,queryObjs,datapath,option)
    end
    
    methods % hook methods
               
        function s = getCustomFieldPulldownInfo(obj) %#ok<MANU>
            % s = getCustomFieldPulldownInfo(obj)
            % Return info to resolve ambiguities when two or more datasets
            % share fields with the same name. Information returned by this
            % method is used by the DataSelector when populating pulldowns for
            % specifying the field in a query.
            %
            % This is a hook method intended to be overridden.
            %
            % s: struct array with fields
            %   * .Name: name of field
            %   * .DataSetName: name of DataSet whose DataSetField object
            %   should be used
            %   * .PrettyName: (optional). Currently not used/implemented.
           
            % Default implementation: none
            s = [];
        end
        
    end
    
    methods (Static) % utilities
        
        function fsp = reseatFileSystemPath(fsp,datapath)
            % fsp = reseatFileSystemPath(fsp,datapath)            
            [~,fsp] = fileparts(fsp);
            fsp = fullfile(datapath,fsp);
        end
        
    end        
end