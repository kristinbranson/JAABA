function makePlatformSpecificBrowserFig(figFile,newFigFile,targetPlatform)
% figFile: name of original .fig file (with extension)
% newFigFile: name of new .fig file (with extension)
% targetPlatform: platform of target.one of {'pc' 'mac' 'unix'}

% Notes
%
% The reason this file exists is because we must simultaneously enable:
% * cross-platform use of the Browser (PC, Unix, Mac)
% * proportional resizing of Browser after startup
% * proportional resizing-to-screen of Browser at startup (for small
% screens)
%
% The latter two points require that we set Units and FontUnits in the
% BrowserGUI fig file to be 'normalized'. However, using a single .fig file
% in this way across platforms leads to cosmetic annoyances. On the other
% hand, in my experience, during GUI development it is easiest to work with
% a single .fig file using 'pixels' for all units.
%
% This suggests a workflow where
% i) GUI Development is done in a single .fig file (with units==pixels)
% ii) After development, platform-specific, resizeable, fig files are
% generated (by this script.)
%
% This file takes a base, "development" .fig file and generates a
% platform-specific .fig file that is resizeable and actually used when
% running the Browser. The base .fig file should have Units='pixels' and
% FontUnits='pixels' for all uicontrols, uipanels, etc. It should be "as
% good as it gets" for cross-platform usage from a single .fig file. This
% script will convert all Units and FontUnits to 'normalized' and resave
% (to a new name), along with making minor platform-specific adjustments as
% necessary.
%
% Run this script _in the desired target platform_, on a computer with a
% screen large enough to open the BrowserGUI at its 'native' pixel
% resolution.

f = open(figFile);
assert(strcmp(get(f,'Type'),'figure'));

assert(nnz(strcmp(targetPlatform,{'pc' 'mac' 'unix'}))==1);

% Set "outer" figure to have units of pixels, so that all screens that are
% large enough to accomodate the "native" size of the GUI open the GUI with
% that same native size. (If the outer figure had eg units of normalized,
% then its size would be scaled to the screensize.)
set(f,'Units','pixels');

allH = findobj(f);
allH = setdiff(allH,f);

for h = allH(:)'

    hType = get(h,'Type');

    if strcmp(targetPlatform,'mac')
        % Mac: tweak popup menus
        if strcmp(hType,'uicontrol') && strcmp(get(h,'Style'),'popupmenu')
            pos = get(h,'Position');
            pos(1) = pos(1)-2; % shift PUM to the left by 2 pixels (we assert that original units are ''pixels'', see below)
            set(h,'Position',pos);
        end
    end

    % Set units to normalized for all children uicontrols/uipanels etc.
    % This enables proportional rescaling of the GUI. 
    if isprop(h,'Units')
        hUnits = get(h,'Units');
        assert(strcmp(hUnits,'pixels'),'''Units'' in original fig file should be ''pixels''.');    
        set(h,'Units','Normalized');        
    end    
    
    % Set FontUnits to normalized     
    if isprop(h,'FontUnits')
        hFontUnits = get(h,'FontUnits');
        assert(strcmp(hFontUnits,'pixels'),'''FontUnits'' in original fig should be ''pixels''.');
        
        switch targetPlatform
            case 'pc'
                % TMW note: Setting a uipanel's FontUnits to normalized on
                % a PC causes it to hang. Similarly, using a PC to open a
                % .fig file that was created on a Mac, in which a uipanel's
                % FontUnits has been saved as 'Normalized', causes the PC
                % to hang.
                if ~strcmp(hType,'uipanel')
                    set(h,'FontUnits','Normalized');
                end
            otherwise
                set(h,'FontUnits','Normalized');
        end
    end
end

hgsave(f,newFigFile);

% close opened GUI
delete(f);
