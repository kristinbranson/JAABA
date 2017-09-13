function y = getNumericScalarFieldSafe(data,fld,mode)
% y = getNumericScalarFieldSafe(data,fld,mode)
% mode: optional, either 'safe' or 'safest' (default).
% y: col vector

if nargin < 3
    mode = 'safest';
end

if isempty(data)
    y = zeros(size(data));
else
    switch mode
        case 'safest'
            y = {data.(fld)}';
            cellfun(@(x)assert((isnumeric(x)||islogical(x))&&isscalar(x)),y);
            y = cell2mat(y);
            assert(isequal(size(y),[numel(data) 1]));
        case 'safe'
            y = cat(1,data.(fld));
            assert(isequal(size(y),[numel(data) 1]));            
        otherwise
            assert(false);
    end
end