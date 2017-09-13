classdef ComparisonPlot < OlyDat.Plot
% ComparisonPlot Compare data from two days.

    properties (Constant)
      SUCCESSTYPE_SUBSETS = {
        'all trials'
        'successful_1grab'
        'successful_2plusgrab'
        'unsuccessful'
        'any success'
        }
      AXISTITLE_FONTSIZE = 12; 
    end
    
    properties
        Name = 'Trial Plot';  
        AllowsMultiStats = false;
        UsesAuxVar = false;
        
        UsesGroupVar2 = true;
        GroupVarLabel = 'Group Var';
        GroupVar2Label = 'Restrict Var';
        
        AutoReplotAfterStatChange = true;
        AutoReplotAfterGroupChange = true;
    end
    
    methods
        
        function descStr = doPlot(obj,ax,data,bstat,grp,plotCfg,expCoordinator) %#ok<INUSL>
          % This is a plot of bstat vs trial, with 2ish grouping vars
          % * primary grouping var is grp/plotCfg.grp. Most typical usage
          % is date. Trials in different groups are represented by
          % different markers.
          % * secondary grouping var is successtype. Trials with different 
          % successtype are represented by different colors. auxVar enables
          % restriction to canned subsets of successtypes.
          
          MARKER = struct('sym',{'.' '^' '+' 'd'},'size',{30 20 20 20},'legsize',{16 6 6 6});
          NGRPSUPPORTED = numel(MARKER);
          COLOR = struct();
          COLOR.SUCCESS1 = ExpPP.TRIALSUCCESSCOLOR.successful_1grab;
          COLOR.SUCCESS2PLUS = ExpPP.TRIALSUCCESSCOLOR.successful_2plusgrab;
          COLOR.FAILURE = [0 0 0];
          XFLD = 'trial';
         
          assert(iscell(grp)&&numel(grp)==2);
          grp = grp{1}; % grp2 for restriction only
          
          % collect all the pieces 
          eid = [data.experiment_id]';
          xAll = [data.(XFLD)]';
          if isnumeric(grp) || islogical(grp)
            grp = arrayfun(@num2str,grp,'uni',0);
          end          
          [tfKeep,tfFailure,tfSuccess1,tfSuccess2] = ComparisonPlot.successSubset(data);
          assert(isequal(size(bstat),size(grp),size(eid),size(xAll),size(tfKeep),size(tfFailure),size(tfSuccess1),size(tfSuccess2))); 
          
%           % Restrict everything based on auxvar
%           eid = eid(tfKeep);
%           bstat = bstat(tfKeep);
%           grp = grp(tfKeep);
%           xAll = xAll(tfKeep);
%           tfFailure = tfFailure(tfKeep);
%           tfSuccess1 = tfSuccess1(tfKeep);
%           tfSuccess2 = tfSuccess2(tfKeep);
          
          % Determine number of groups to include in plot
          ungrp = unique(grp);
          Ngrp = numel(ungrp);
          if Ngrp > NGRPSUPPORTED
            warningNoTrace('ComparisonPlot:tooManyGroups',...
              'Too many groups selected. Only showing data for %d groups: %s',...
            NGRPSUPPORTED,civilizedStringFromCellArrayOfStrings(ungrp(1:NGRPSUPPORTED)));
          end
          Ngrpplot = min(Ngrp,NGRPSUPPORTED);

          hold(ax,'on');

          hPlotLeg = nan(Ngrpplot,1); % for legend
          legLbl = cell(Ngrpplot,1); % for legend
          allX = zeros(0,1);
          allY = zeros(0,1);
          allEid = zeros(0,1);
          for iGrp = 1:Ngrpplot
            g = ungrp{iGrp};
            tfGrp = strcmp(g,grp);
            tfGrpF = tfGrp & tfFailure;
            tfGrpS1 = tfGrp & tfSuccess1;            
            tfGrpS2 = tfGrp & tfSuccess2;            
            
            xGrpF = xAll(tfGrpF);
            yGrpF = bstat(tfGrpF);
            xGrpS1 = xAll(tfGrpS1);
            yGrpS1 = bstat(tfGrpS1);
            xGrpS2 = xAll(tfGrpS2);
            yGrpS2 = bstat(tfGrpS2);   
            marker = MARKER(iGrp);
            plot(ax,xGrpF,yGrpF,marker.sym,'markersize',marker.size,'Color',COLOR.FAILURE);
            plot(ax,xGrpS1,yGrpS1,marker.sym,'markersize',marker.size,'Color',COLOR.SUCCESS1);
            plot(ax,xGrpS2,yGrpS2,marker.sym,'markersize',marker.size,'Color',COLOR.SUCCESS2PLUS);            
            
            allX = [allX; xAll(tfGrp)]; %#ok<AGROW>
            allY = [allY; bstat(tfGrp)]; %#ok<AGROW>
            allEid = [allEid; eid(tfGrp)]; %#ok<AGROW>
              
            % plot a (nan,nan) point in black for leg
            hPlotLeg(iGrp) = plot(ax,nan,nan,marker.sym,'markersize',marker.legsize,'Color',COLOR.FAILURE);
            legLbl{iGrp} = g;                                   
          end
          % legend; plot (nan,nan) points for colors
          tmpFlds = fieldnames(COLOR);
          for f = tmpFlds(:)',f=f{1}; %#ok<FXSET>
            hPlotLeg(end+1,1) = plot(ax,nan,nan,'s','Color',COLOR.(f),'MarkerFaceColor',COLOR.(f)); %#ok<AGROW>
            legLbl{end+1,1} = lower(f); %#ok<AGROW>
          end
  
          xl = xlim(ax);
          xl(1) = 0;
          xlim(ax,xl);
          legend(ax,hPlotLeg,legLbl,'location','best','interpreter','none');
          hold(ax,'off');
          grid(ax,'on');                    
          ylabel(ax,plotCfg.stat.PrettyName,'interpreter','none','fontsize',ComparisonPlot.AXISTITLE_FONTSIZE);
          xlabel(ax,XFLD,'interpreter','none','fontsize',ComparisonPlot.AXISTITLE_FONTSIZE);          
          titlestr = BoxPlotPlot.titleStr(plotCfg,numel(bstat));
          title(ax,titlestr,'interpreter','none','fontsize',ComparisonPlot.AXISTITLE_FONTSIZE);
          
          OlyDat.XYPlotClickHandler(ax,allX,allY,allEid,expCoordinator,true);
          
          descStr = {'We can put stats here etc'};
        end
        
    end
    
    methods (Static)
      
      function [tfKeep,tfFailure,tfSuccess1,tfSuccess2] = successSubset(data)
        
        SUCCESSFIELD = 'auto_Grab_success';
        SUCCESSTYPEFIELD = 'auto_Grab_successtype';
        
        successEnum = {data.(SUCCESSTYPEFIELD)}';
        tfFailure = strcmp(successEnum,'unsuccessful');
        tfSuccess1 = strcmp(successEnum,'successful_1grab');
        tfSuccess2 = strcmp(successEnum,'successful_2plusgrab');
        assert(all(tfFailure | tfSuccess1 | tfSuccess2));

        % sideline: check SUCCESSFIELD
        tfSuccessTmp = [data.(SUCCESSFIELD)]';
        assert(isequal(tfSuccessTmp,tfSuccess1 | tfSuccess2));

        tfKeep = true(size(data));
%         switch plotCfgAux.Name
%           case 'all trials',           tfKeep = true(size(data));
%           case 'successful_1grab',     tfKeep = tfSuccess1;
%           case 'successful_2plusgrab', tfKeep = tfSuccess2;
%           case 'unsuccessful',         tfKeep = tfFailure;
%           case 'any success',          tfKeep = tfSuccess1 | tfSuccess2;
%           otherwise,                   assert(false,'Unknown successtype subset.');
%         end
        
      end
      
    end    
    
end