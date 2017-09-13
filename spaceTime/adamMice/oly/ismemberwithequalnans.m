function tf = ismemberwithequalnans(a,b)

if isnumeric(a) && isvector(a) && isnumeric(b) && isvector(b)
  tf = ismember(a,b) | isnan(a) & any(isnan(b));
else
  tf = ismember(a,b);
end