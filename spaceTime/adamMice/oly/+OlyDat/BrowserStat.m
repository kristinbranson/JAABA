classdef BrowserStat < handle
%BrowserStat OlyDat Browser Statistic
%   A BrowserStat object represents a fly olympiad assay statistic.
%   BrowserStats are used by the OlyDat Browser and provide metadata
%   relevant to plotting.
%
%   See also OlyDat.BrowserStat/BrowserStat, OlyDat.Browser, OlyDat.Plot

% TODO relationship to SAGE.DataSetField
    
    properties
        Name % Name of statistic. Typically, this is the fieldname in SAGE, or the fieldname in the assay's data structure.
        PrettyName % Name of the statistic as listed in Browser pull-down menus, etc.
        ValidValues % List of valid values for enumerated stats. Either a col cellstr, or a col numeric array.
        tfBuriedStat % (logical) If true, the statistic does not appear at the top-level of the assay's data structure.
    end
    
    methods
        
        function name = getDataFieldName(obj) 
            % name = getDataFieldName(obj)
            % Returns the fieldname the OlyDat Browser will use to access 
            % data structure. 
            %
            % This method is intended to be overridden,
            % typically/especially when tfBuriedStat is true. The default
            % implementation returns obj's .Name property.

            name = obj.Name;
        end

    end
    
    methods
        
        function obj = BrowserStat(name,pname,validvalues,tfBuried)
        %BrowserStat BrowserStat constructor.
        %   obj = BrowserStat(name) constructs a BrowserStat object with
        %   the given Name. The object's PrettyName is set to the same
        %   name. ValidValues is set to empty and tfBuriedStat is set to
        %   false.
        % 
        %   obj = BrowserStat(name,pname) constructs a BrowserStat object
        %   with the given Name/PrettyName. ValidValues and tfBuried
        %   default as above.
        %
        %   obj = BrowserStat(name,pname,validvalues)
        %   obj = BrowserStat(name,pname,validvalues,tfBuried)
        %
        %   Note: any subclasses of BrowserStat should accept a constructor
        %   signature obj = ctor(name,prettyname). This is used by the
        %   OlyDat Browser.
        
            assert(ischar(name));
            obj.Name = name;
            if nargin < 2
                obj.PrettyName = name;
            else
                assert(ischar(pname));
                obj.PrettyName = pname;
            end
            if nargin < 3
                obj.ValidValues = cell(0,1);
            else                    
                assert(iscellstr(validvalues) || isnumeric(validvalues));
                obj.ValidValues = validvalues(:);
            end
            if nargin < 4
                obj.tfBuriedStat = false;
            else
                validateattributes(tfBuried,{'logical'},{'scalar'});
                obj.tfBuriedStat = tfBuried;
            end
        end
        
    end
end
        