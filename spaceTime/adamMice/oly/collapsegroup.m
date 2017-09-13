function [gg h] = collapsegroup(g)
% [gg h] = collapsegroup(g)
% g: group values. col vec of ints, or col cellstr. 
% gg: col vec of ints taking values in the interval [1,Ngroups].
% h: if g is a col int vec, h is an int col vec. if g is a col cellstr, h
% is a col cellstr. h is an Ngroups-by-1 vector.
%
% gg(i) is the new grouping value for g(i). h(gg(i)) is equal to g(i).
%
% This is similar to stats/grp2idx, but does something I prefer with nans.

assert(isnumeric(g)||iscellstr(g));
g = g(:);

h = unique(g);

% nans are not uniqued; if there is a nan, include one in h
if isnumeric(h) && any(isnan(h))
    h(isnan(h),:) = [];
    h(end+1,1) = nan;
end

[~,gg] = ismember(g,h);

% if there are nans in g/h, gg will have 0 values since a nan in g is not
% considered to match up with a nan in h.
if isnumeric(h) && any(isnan(h))
    gg(gg==0) = numel(h);
end
