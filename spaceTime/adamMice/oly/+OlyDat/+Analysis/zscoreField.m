function d = zscoreField(d,fld,type,postfix)
% d = zscoreField(d,fld,type,postfix)
% compute zscored field.
%
% d: data struct
% fld: field
% type: either 'zscore' or 'iqr'. Defaults to 'zscore'.
% postfix: if type is 'zscore', this defaults to 'zs'. if type is 'iqr',
% this defaults to 'iqr'.

assert(isstruct(d)&&isfield(d,fld));
if nargin < 3
    type = 'zscore';
end
switch type
    case 'zscore'
        postfix = 'zs';
    case 'iqr'
        postfix = 'iqr';
    otherwise
        assert(false);
end

y = OlyDat.Analysis.getNumericScalarFieldSafe(d,fld);
switch type
    case 'zscore'
        zsy = (y-nanmean(y))/nanstd(y);
    case 'iqr'
        zsy = (y-nanmean(y))/iqr(y);
end

newfld = [fld '_' postfix];
assert(~isfield(d,newfld));
zsy = num2cell(zsy);
[d.(newfld)] = deal(zsy{:});
