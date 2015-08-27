function [pfdata,pfunits] = readPFData(pfname,flies)

ss = fopen(pfname,'r');
qq = fread(ss,10,'*char')';
fclose(ss);

% Doing matfile is slower for more than 4 flies
if strcmp(qq(8:10),'7.3') && numel(flies)<4
  m = matfile(pfname);
  pfdata = cell(1,numel(flies));
  for ndx = 1:numel(flies);
    pfdata(ndx) = m.data(1,flies(ndx));
  end
  pfunits = m.units;
else
  m = load(pfname);
  pfdata = m.data(flies);
  pfunits = m.units;
end