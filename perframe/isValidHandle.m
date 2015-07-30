function valid = isValidHandle(harray)
if verLessThan('matlab','8.4.0'),
  valid = ~isnan(harray);
else
  valid = ishandle(harray);
end