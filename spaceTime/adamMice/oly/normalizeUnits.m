function normalizeUnits(gdata,fontUnitsExcludeTags)
% normalizeUnits(gdata)
% Change 'Units' properties of all handles contained in gdata to
% 'Normalized' (the 'FontUnits' property is unaffected). gdata is a
% structure with fields whose values are handles to GUI components (eg as
% returned by guidata).

if nargin < 2
    fontUnitsExcludeTags = cell(0,1);
end

comps = fieldnames(gdata);
for c = 1:numel(comps)
    tag = comps{c};
    h = gdata.(tag);
    if isprop(h,'Units')
        set(h,'Units','Normalized');
    end
    if isprop(h,'FontUnits') && ~any(strcmp(tag,fontUnitsExcludeTags))
        set(h,'FontUnits','Normalized');
    end
end

end