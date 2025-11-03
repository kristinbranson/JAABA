function [pfdata,pfunits,update] = readPFData(pfname,flies,params)

% disabling mfile version, not faster, seems to be slowing things down even
% to check

if nargin < 3
    params = struct;
    do_check = false;
else
    do_check = true;
end

saved_params = struct;
update = false;
m = load(pfname);
pfdata = m.data(flies);
pfunits = m.units;
if isfield(m,'params')
  saved_params = m.params;
else
  do_check = false;
end

if ~isequal(saved_params,params) && do_check ...
    && ~isempty(saved_params) && ~isempty(params)
    update = true;
end
