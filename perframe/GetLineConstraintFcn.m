function fcn = GetLineConstraintFcn(xlim,ylim,origPosition,xLocs)
            
fcn = @constrainLineToRect;

%Store previous position matrix for use in constrainLineToRect.
line_pos_last = origPosition;
    %-----------------------------------------
    function new_pos = constrainLineToRect(pos)
                                
        previous_position_cached = ~isempty(line_pos_last);
        
        is_end_point_drag = previous_position_cached &&...
                            (any(pos(:,1) == line_pos_last(:,1)) &&...
                            any(pos(:,2) == line_pos_last(:,2)));
                                    
        if is_end_point_drag
            new_pos = line_pos_last;
%             new_pos = [constrainPointToRect(pos(1,:));constrainPointToRect(pos(2,:))];
        else
            %Apply correction made to first end point to both end points
            constrained_p1 = constrainPointToRect(pos(1,:));
            v1 = constrained_p1 - pos(1,:);
            temp_pos = pos + [v1; v1];
            
            % Now reconstrain both end points according to correction made
            % to second endpoint of partially constrained line.
            constrained_p2 = constrainPointToRect(temp_pos(2,:));
            v2 = constrained_p2 - temp_pos(2,:);
            new_pos = temp_pos + [v2; v2];
            
            [~,closestXlocs] = min( abs(temp_pos(1)-xLocs));
            new_pos(:,1) = xLocs(closestXlocs);
        end
        
        line_pos_last = new_pos;
    end

    function new_pos = constrainPointToRect(pos)

        x_candidate = pos(1);
        y_candidate = pos(2);

        x_new = min( xlim(2), max(x_candidate, xlim(1)) );
        y_new = min( ylim(2), max(y_candidate, ylim(1)) );

        new_pos = [x_new y_new];

    end

end
