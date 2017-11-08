classdef CurationField < handle
%CurationField OlyDat Browser Curation Field
%   A CurationField represents a field used in manual curation.
%
%   See also OlyDat.BrowserStat/BrowserStat, OlyDat.Browser

    
    properties
        Name % Name of the statistic in SAGE
        TFEditable; % Logical. If true, stat can be marked during curation.
        
        TableName % Name of the statistic to be presented in Curation Manager Table. To leave field out of table, use an empty name ''.
        TableColWidth % Column width in curation manager table. Defaults to 0 indicating auto-size.
        
        TSVName % Name used in Curation TSV file header. To leave field out of the TSV file, use an empty name ''.
        ValidValues; % (For enumerated stats) Cellstr of valid options
    end
    
    methods
    
        function obj = CurationField(name,tfEditable,tableName,tableColWidth,tsvName,validvals)
        %CurationField CurationField constructor.
        %   obj = CurationField(name,tfEditable)
        %   obj = CurationField(name,tfEditable,tableName,tableColWidth,tsvName,validvals)        
        
            error(nargchk(2,6,nargin,'struct'));
            if nargin < 3 || isnumeric(tableName) && isequal(tableName,[])
                tableName = name;
            end
            if nargin < 4 || isempty(tableColWidth)
                tableColWidth = 0;
            end
            if nargin < 5 || isnumeric(tsvName) && isequal(tsvName,[])
                tsvName = name;
            end            
            if nargin < 6 || isempty(validvals)
                validvals = cell(0,1);
            end
            assert(ischar(name));
            validateattributes(tfEditable,{'logical'},{'scalar'});
            assert(ischar(tableName));
            validateattributes(tableColWidth,{'numeric'},{'integer' 'nonnegative' 'scalar'});
            assert(ischar(tsvName));
            assert(iscellstr(validvals));
            
            obj.Name = name;
            obj.TFEditable = tfEditable;
            obj.TableName = tableName;
            obj.TableColWidth = tableColWidth;
            obj.TSVName = tsvName;
            obj.ValidValues = validvals;
        end
                
        function tfValid = validateValue(obj,val)
            assert(isscalar(obj));
            if isempty(obj.ValidValues)
                tfValid = true;
            else
                tfValid = any(strcmp(val,obj.ValidValues));
            end
        end
    end
end
        