function y = getFieldSafe(data,fld)
% y = getFieldSafe(data,fld)
% mode: optional, either 'safe' or 'safest' (default).
% y: col cell

if isempty(data)
    y = cell(0,1);
else
    y = {data.(fld)}';
    assert(isequal(size(y),[numel(data) 1]));
end