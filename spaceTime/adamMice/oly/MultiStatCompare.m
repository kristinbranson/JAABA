classdef MultiStatCompare < MultiStatCompareBase
  % MultiStatCompare Compare predictive strength of stats on a binary
  % grouping var
  
  properties
    Name = 'Multi-Stat Comparison';
    AllowsMultiStats = true;
    UsesAuxVar = false;
    
    UsesGroupVar2 = true;
    GroupVarLabel = 'Observed Var';
    GroupVar2Label = 'Restrict Var';
    
    AutoReplotAfterStatChange = true;
    AutoReplotAfterGroupChange = true;
  end
  
  methods
    
    function descStr = doPlot(obj,ax,data,bstat,grp,plotCfg,~)
      % * primary grouping var is grp/plotCfg.grp. Must be a binary
      % grouping vector, ie grp can only take on two values.
      % * auxVar enables restriction to canned subsets of successtypes.
                      
      statnames = {plotCfg.stat.Name}';
      assert(isequal(size(bstat),[numel(data) numel(statnames)]));      
      
      assert(iscell(grp) && numel(grp)==2,'Uses groupvar2 for data restriction');
      grp1 = grp{1}; %grp2 just used for restriction
      
      % check/massage grp
      if islogical(grp1)
        grp1 = double(grp1);
      end
      [grp1,grpnames] = collapsegroup(grp1);
      if max(grp1)==1
        OlyDat.Plot.messageAxisCenter(ax,'Grouping variable only takes on one value for available data.');
        descStr = '';
        return;
      elseif max(grp1)>2
        OlyDat.Plot.messageAxisCenter(ax,'Grouping variable takes on more than two distinct values.');
        descStr = '';
        return;
      end
      grp1 = grp1-1; % needs to be 0-based
      if isnumeric(grpnames)
        grpnames = arrayfun(@num2str,grpnames,'uni',0);
      end
      
%       % restrict to subset
%       tfKeep = ComparisonPlot.successSubset(data,plotCfg.aux);
%       assert(isequal([size(bstat,1) 1],size(grp1),size(tfKeep)));      
%       bstat = bstat(tfKeep,:);
%       grp1 = grp1(tfKeep);
      
      newax = obj.setupAxes(ax);      
      
      grpvarstr = BoxPlotPlot.groupVarStr(plotCfg.grpVars);
      grpvarstr2 = BoxPlotPlot.groupVarStr(plotCfg.grpVars2);
      titlestr1 = sprintf('%d stats, N=%d. Observed: %s {%s}.',...
        size(bstat,2),size(bstat,1),...
        plotCfg.grp.PrettyName,grpvarstr);
      titlestr2 = sprintf('Restriction: %s {%s}. ',...
        plotCfg.grp2.PrettyName,grpvarstr2);
      titlestr = {titlestr1;titlestr2};
      CompareStatsAcrossGroups(bstat,grp1,statnames,grpnames,'axes',newax,'title',titlestr);
      
      descStr = {'We can put stats here etc'};
    end
    
  end
      
end