function tweakGUILayout(handles)

switch computer
    case 'MACI64'
        zlclTweakPUMsForMacI64(handles);
end

end

function zlclTweakPUMsForMacI64(handles)
tags = fieldnames(handles);
pumTags = tags(strncmp(tags,'pum',3));
for c = 1:numel(pumTags)
    h = handles.(pumTags{c});
    origUnits = get(h,'Units');
    set(h,'Units','pixels');
    pos = get(h,'Position');
    pos(1) = pos(1)-2; % shift PUM to the left by 2 pixels
    set(h,'Position',pos);
    set(h,'Units',origUnits);
end

end