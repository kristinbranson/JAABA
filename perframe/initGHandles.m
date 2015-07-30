function hh = initGHandles(sz)

if verLessThan('matlab','8.4.0'),
 hh = nan(sz);
else
  hh = gobjects(sz);
end