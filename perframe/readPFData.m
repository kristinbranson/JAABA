function [pfdata,pfunits,update] = readPFData(pfname,flies,params)

ss = fopen(pfname,'r');
qq = fread(ss,10,'*char')';
fclose(ss);
if nargin < 3
    params = struct;
    do_check = false;
else
    do_check = true;
end

% Doing matfile is slower for more than 4 flies
saved_params = struct;
update = false;
if strcmp(qq(8:10),'7.3') && numel(flies)<4
  m = matfile(pfname);
  pfdata = cell(1,numel(flies));
  for ndx = 1:numel(flies);
    pfdata(ndx) = m.data(1,flies(ndx));
  end
  pfunits = m.units;
  if isprop(m,'params')
      saved_params = m.params;
  else
      do_check = false;
  end
else
  m = load(pfname);
  pfdata = m.data(flies);
  pfunits = m.units;
  if isfield(m,'params')
      saved_params = m.params;
  else
      do_check = false;
  end
end

if ~isequal(saved_params,params) && do_check ...
    && ~isempty(saved_params) && ~isempty(params)
    update = true;
end
