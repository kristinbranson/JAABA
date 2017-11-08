classdef BoxPlotPlot < OlyDat.Plot
  % BoxPlotPlot Plot distributions of a single stat using a boxplot.
  %
  %   This plot uses the boxplot command from the MATLAB stats toolbox to
  %   plot the distributions of a single, grouped variable.
  %
  %   At the moment this plot does not have interactivity (it is not
  %   clickable, nor does it respond to experiment selection).
  %
  %   The Kruskal-Wallis statistic is used to probe for differences between
  %   groups.
  %
  %   See OlyDat.Plot, boxplot.
  
  properties (Constant)
    STATTEXT_ARGS = {'HorizontalAlignment','center','FontUnits','pixels','FontSize',12,'FontWeight','normal'};
  end
  
  properties
    Name = 'Box Plot';
    AllowsMultiStats = false;
    UsesAuxVar = false;
    
    UsesGroupVar2 = true;    
    GroupVarLabel = 'Group Var';
    GroupVar2Label = 'Group Var';
    
    AutoReplotAfterStatChange = true;
    AutoReplotAfterGroupChange = true;

    BoxPlotPVArgs = {};
  end
  
  methods
    
    % plotObj = BoxPlotPlot()
    % plotObj = BoxPlotPlot(pvArgs)
    %
    % pvArgs is an optional cell array of PV arguments to be passed to
    % boxplot. See 'help boxplot' for more information.
    function obj = BoxPlotPlot(pvArgs)
      if nargin >= 1
        assert(iscell(pvArgs) && mod(numel(pvArgs),2)==0 && ...
          iscellstr(pvArgs(1:2:end)),...
          'pvArgs must be a cell array of additional PV-arguments for boxplot.');
        obj.BoxPlotPVArgs = pvArgs;
      end
    end
    
    function descStr = doPlot(obj,ax,data,bstat,grp,plotCfg,~)       
                      
      assert(iscell(grp)&&numel(grp)==2);
      grp1 = grp{1};
      grp2 = grp{2};
      grp1C = categorical(grp1);
      grp2C = categorical(grp2);    
      grp1cats = categories(grp1C);
      grp2cats = categories(grp2C);
      
%       % restrict to subset
%       tfKeep = ComparisonPlot.successSubset(data,plotCfg.aux);
%       assert(isequal(size(bstat),size(grp1C),size(grp2C),size(tfKeep)));
%       bstat = bstat(tfKeep);
%       grp1C = grp1C(tfKeep);
%       grp2C = grp2C(tfKeep);      
      
      grpBig = {grp1C grp2C};      
      
      % Currently using two separate grouping vars for boxplot, in order
      % to support 'fullfactors'=true. Below we create a single combined
      % grouping var for computing stats
      boxplot(ax,bstat,grpBig,'labelorientation','inline','fullfactors',true);
                 
      % label stats on plot
      [ns,mds,p25,p75,gname] = grpstats(bstat,grpBig,{'numel',@median,...
        @(x)prctile(x,25),@(x)prctile(x,75),'gname'});
      yl = ylim(ax);
      dy = diff(yl);
      yl(2) = yl(2) + .15*dy;
      ylim(ax,yl);
      dytxt = 0.035*dy;
      ypos = yl(2)-dytxt;

      xLbl = 0;
      for iGrp1 = 1:numel(grp1cats)
      for iGrp2 = 1:numel(grp2cats)
        
        gstr1 = grp1cats{iGrp1};
        gstr2 = grp2cats{iGrp2};
        xLbl = xLbl+1;        
        % We iterate over the groups as shown/ordered in the boxplot. grpstats()
        % does not have a 'fullfactors' option, so it will not contain
        % groups for which there is no data.
        
        tfStat = strcmp(gstr1,gname(:,1)) & strcmp(gstr2,gname(:,2));
        if any(tfStat)        
          assert(nnz(tfStat)==1);
          nstr = sprintf('N: %d',ns(tfStat));
          mdstr = sprintf('mdn: %.2f',mds(tfStat));
          p25str = sprintf('p25: %.2f',p25(tfStat));
          p75str = sprintf('p75: %.2f',p75(tfStat));

          text(xLbl,ypos,nstr,obj.STATTEXT_ARGS{:});
          text(xLbl,ypos-dytxt,mdstr,obj.STATTEXT_ARGS{:});
          text(xLbl,ypos-2*dytxt,p25str,obj.STATTEXT_ARGS{:});
          text(xLbl,ypos-3*dytxt,p75str,obj.STATTEXT_ARGS{:});
        end
      end            
      end
      
      titlestr = obj.titleStr(plotCfg,numel(bstat));
      title(ax,titlestr,'interpreter','none','fontsize',ComparisonPlot.AXISTITLE_FONTSIZE);
      ylabel(ax,plotCfg.stat.PrettyName,'interpreter','none','fontweight','normal','fontsize',ComparisonPlot.AXISTITLE_FONTSIZE);
      xlabelstr = sprintf('%s/%s',plotCfg.grp.PrettyName,plotCfg.grp2.PrettyName);
      xlabel(ax,xlabelstr,'interpreter','none','fontweight','normal','fontsize',ComparisonPlot.AXISTITLE_FONTSIZE);
      
      % create single 
      gComb = strcat(cellstr(grp1C),'#',cellstr(grp2C));
      pkw = kruskalwallis(bstat,gComb,'off');
      %pan = anova1(bstat,gComb,'off');
      
      [ptmat,ptmatlbl] = BoxPlotPlot.ttestMatrix(bstat,gComb);
      gOrder = strcat(gname(:,1),'#',gname(:,2));
      [tf,loc] = ismember(gOrder,ptmatlbl);
      assert(all(tf));
      ptmat = ptmat(loc,loc);
      ptmatlbl = ptmatlbl(loc);
      descStr = [...
        {sprintf('KruskalWallis p-value: %.4g',pkw)};... %{sprintf('1D Anova p-value: %.4g',pan)};...
        {''};...
        {'Pairwise t-test (uncorrected):'};...
        sprintf('  %s\n',ptmatlbl{:});...
        num2str(ptmat,3)];
      
      zoom reset;
    end
    
  end
  
  methods (Static)
    
    function tstr = titleStr(plotCfg,nbstat)
      grpvarstr = BoxPlotPlot.groupVarStr(plotCfg.grpVars);
      grpvarstr2 = BoxPlotPlot.groupVarStr(plotCfg.grpVars2);
      tstr1 = sprintf('%s, N=%d. Group1: %s {%s}.',...
        plotCfg.stat.PrettyName,nbstat,... 
        plotCfg.grp.PrettyName,grpvarstr);
      tstr2 = sprintf('Group2: %s {%s}.',...
        plotCfg.grp2.PrettyName,grpvarstr2);
      tstr = {tstr1;tstr2};
    end
    
    function str = groupVarStr(grpvars)
      if isnumeric(grpvars)
        grpvars = arrayfun(@num2str,grpvars,'uni',0);
      end
      str = sprintf('%s,',grpvars{:});
      str = str(1:end-1);
      
      MAXSTR = 46;
      if numel(str)>MAXSTR
        str = sprintf('<%d values>',numel(grpvars));
      end
    end
    
    function [pmat,glbl] = ttestMatrix(y,g)
      assert(isvector(y)&&isnumeric(y)&&isvector(g)&&numel(y)==numel(g));
      
      [gi,glbl] = collapsegroup(g);
      if isnumeric(glbl)
        glbl = arrayfun(@num2str,glbl,'uni',0);
      end
      ngrp = max(gi);
      ycell = arrayfun(@(zzI)y(gi==zzI),1:ngrp,'uni',0); % ycell{i} contains all y data for gi==i
      
      pmat = nan(ngrp,ngrp);
      for i = 1:ngrp
      for j = 1:ngrp
        [~,pmat(i,j)] = ttest2(ycell{i},ycell{j},'vartype','unequal'); 
      end
      end      
    end
    
  end
  
end
