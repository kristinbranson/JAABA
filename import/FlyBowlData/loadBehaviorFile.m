
function behaviorSet = loadBehaviorFile(filename)
    fid = fopen(filename);
    if fid < 0
        error('Unable to open behavior file ''%s''\n', filename);
    end
    
    subset = '';
    while(~feof(fid))
        tline = fgetl(fid);
        if isempty(tline)
            continue;
        end
        if(tline(1)=='*')
            tokens = textscan(tline, '%s', 'delimiter', ' ');
            subset = tokens{1}{1}(2:end);
            bi = 1;
            continue;
        end
        if isempty(subset)
            continue;
        end

        tokens = textscan(tline, '%s', 'delimiter', '\t');
        if length(tokens{1}) < 4
            continue;
        end
        longDesc = tokens{1}{1};
        pos = regexp(longDesc, '<|>');
        if ~isempty(pos)
            longDesc = longDesc(1 : min(pos)-1);
        end
        [startIdx endIdx] = regexp(longDesc, 'Left |Right ');
        if ~isempty(startIdx)
            longDescM = {longDesc(1 : endIdx(1))};
            longDesc = longDesc(endIdx(1)+1 : end);
            [~, endIdx] = regexp(longDesc, 'Abdominal |Thoracic |Tip ');
            if ~isempty(endIdx)
                longDescM = {longDescM{1}; longDesc(1 : endIdx(1)); longDesc(endIdx(1)+1 : end)};
            end
        else
            longDescM = {longDesc};
        end
        behaviorSet.(subset)(bi).longDesc = tokens{1}{1};
        behaviorSet.(subset)(bi).longDesc_display = longDescM;
        behaviorSet.(subset)(bi).short = tokens{1}{2};
        behaviorSet.(subset)(bi).value1 = sscanf(tokens{1}{3}, '%x');
        behaviorSet.(subset)(bi).value2 = sscanf(tokens{1}{4}, '%x');     
        bi = bi + 1;
    end
    fclose(fid);
end