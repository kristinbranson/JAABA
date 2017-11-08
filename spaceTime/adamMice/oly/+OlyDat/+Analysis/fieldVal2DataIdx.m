function m = fieldVal2DataIdx(data,fld)
% m = fieldVal2DataIdx(data,fld)
% data: data struct
% fld: data fieldname
% 
% m: containers.map. Keys: distinct fieldnames. Vals: index arrays into
% data
%
% Example:
%   m = OlyDat.Analysis.fieldVal2DataIdx(data,'line_name');
%

assert(ischar(fld));
assert(isstruct(data) && isfield(data,fld));

if isempty(data)
    m = containers.Map();
    % assume char keytype
    return;
end

f1 = data(1).(fld);
if ischar(f1)
    fldType = 'char';
elseif isnumeric(f1) && isscalar(f1)
    fldType = 'double';
else
    assert(false,'Field has unsupported value type.');
end

m = containers.Map('KeyType',fldType,'ValueType','any');

for dIdx = 1:numel(data)
    v = data(dIdx).(fld);
    switch fldType
        case 'char'
            assert(ischar(v),'Field values have mixed types.');
        case 'double'
            assert(isnumeric(v)&&isscalar(v),'Field values have mixed types.');
    end
    
    if ~m.isKey(v)
        m(v) = zeros(0,1);
    end
    
    idx = m(v);
    idx(end+1,1) = dIdx; %#ok<AGROW>
    m(v) = idx;
end
               