classdef RepeatMultiStatCompare < MultiStatCompareBase & RepeatPlotBase
  
  properties
    Name = 'Repeat Multi-Stat Comparison';
    AllowsMultiStats = true;
    UsesAuxVar = false;
    
    UsesGroupVar2 = true;
    GroupVarLabel = 'Fixed Effect';
    GroupVar2Label = 'Restrict Var';
    
    AutoReplotAfterStatChange = false;
    AutoReplotAfterGroupChange = false;
  end
  
  methods        
            
    function obj = RepeatMultiStatCompare(hObj)
      assert(isa(hObj,'HandleObj'));
      obj = obj@RepeatPlotBase(hObj);
    end
    
    function descStr = doPlot(obj,ax,data,bstat,grp,plotCfg,~)
      
      [ok,~,tfKeep,datesAllC,gRptC,gCtlExpWashC] = obj.collectData(data);
      if ~ok
        descStr = obj.FailedRepeatSelection(ax);
        return;
      end
      assert(isequal(size(tfKeep),[size(bstat,1) 1]));
      bstat = bstat(tfKeep,:);
      
      gFixedGrpName = plotCfg.grp.Name;
      gFixed = grp{1};
      assert(isequal(size(tfKeep),size(gFixed)));
      gFixed = gFixed(tfKeep);
      
      if isempty(gFixedGrpName)
        % No group1var selected. Default fixed effect: Ctl-Exp compare
        gFixed = gCtlExpWashC;
        gFixedVals = {'ctl' 'exp'};
      else
        % User-specified fixed effect
        gFixed = categorical(gFixed);
        gFixedVals = categories(gFixed);
        if numel(gFixedVals)~=2
          OlyDat.Plot.messageAxisCenter(ax,...
            sprintf('Fixed effect ''%s'' must take on precisely two values.',gFixedGrpName));
          descStr = '';
          return;
        end
      end
         
      newax = obj.setupAxes(ax);
      bstatnames = {plotCfg.stat.Name}';
      descStr = obj.doPlotCore(newax,bstat,bstatnames,gRptC,gFixed,gFixedVals,datesAllC,plotCfg);     
    end
    
  end
  
  methods (Static)
    
    function descStr = doPlotCore(ax,X,Xnames,gRptC,gFixed,gFixedCats,datesC,plotCfg)
      % gFixedCats: 2-el cellstr, categories of gFixed to use in analysis
      % datesC, plotCfg: used only for legend/labeling
      
      assert(numel(ax)==2 && all(ishandle(ax)));
      [ntrials,nstats] = size(X);     
      assert(iscellstr(Xnames) && isvector(Xnames) && numel(Xnames)==nstats);
      assert(iscategorical(gRptC) && isvector(gRptC) && numel(gRptC)==ntrials);
      assert(iscategorical(gFixed) && isvector(gFixed) && numel(gFixed)==ntrials);
      assert(iscellstr(gFixedCats) && numel(gFixedCats)==2 && all(ismember(gFixedCats,categories(gFixed))));
      assert(iscategorical(datesC) && isvector(datesC) && numel(datesC)==ntrials);
      nrpt = numel(unique(cellstr(gRptC)));
      
      tf0 = gFixed==gFixedCats{1};
      tf1 = gFixed==gFixedCats{2};
      tf01 = tf0 | tf1;
      
      % Look for stats with at least 2 non-nan values in Ctls and Exps 
      n0 = sum(~isnan(X(tf0,:)),1);
      n1 = sum(~isnan(X(tf1,:)),1);
      tfGoodStats = n0 > 1 & n1 > 1;
      tfBadStats = ~tfGoodStats;
      if all(tfBadStats)
        warningNoTrace('CompareStats:noData',...
          'There are no stats that have at least two control and experimental datapoints.');
        descStr = 'No good stats';
        return;
      elseif any(tfBadStats)
        badStatNames = Xnames(tfBadStats);
        badStatNamesStr = civilizedStringFromCellArrayOfStrings(badStatNames);
        warningNoTrace('CompareStats:noData',...
          'Skipping stats that don''t have at least two control and experimental datapoints: %s',...
          badStatNamesStr);
      end
      % Restrict to goodstats
      X = X(:,tfGoodStats);
      Xnames = Xnames(tfGoodStats);
      nstats = size(X,2);
      
      % For each stat, compute p-value for C vs W, ignoring 
      gtmp = nan(size(tf0));
      gtmp(tf0) = 0;
      gtmp(tf1) = 1;      
      pval.kbttest = kbttest(X(tf01,:),gtmp(tf01));
      [~,pval.ttest] = arrayfun(@(i)ttest2(X(tf0,i),X(tf1,i),'vartype','unequal'),1:nstats);
      pval.ranksum = arrayfun(@(i)ranksum(X(tf0,i),X(tf1,i)),1:nstats);
      
      % Now compute pvals accounting for repeats
      pval.aovInt = nan(nstats,1);
      pval.aovCEMain = nan(nstats,1);
      pval.aovLmeCE = nan(nstats,1);
      pval.aovLme2CE = nan(nstats,1);
      pval.aov1CE = nan(nstats,1);
      hWB = waitbar(0,'Computing p-values');
      for i = 1:nstats
        try
          [pval.aovInt(i),pval.aovCEMain(i),pval.aovLmeCE(i),pval.aovLme2CE(i),pval.aov1CE(i)] = ...
            Repeat.anls(X(tf01,i),cellstr(gRptC(tf01)),cellstr(gFixed(tf01)));
        catch ME
          warningNoTrace('RepeatMultiStatCompare:anlsFail',...
            'Repeat-analysis had trouble with stat %d, ''%s'': %s',...
            i,Xnames{i},ME.message);
          % all pvals will be nan
        end
        waitbar(i/nstats);
      end
      delete(hWB);
      
      % sort by aovLmeCE
      tfNanP = isnan(pval.aovLmeCE);
      if any(tfNanP)
        fprintf(2,'%d nan LME p-vals\n',nnz(tfNanP));
      end
      [~,idxPval] = sort(pval.aovLmeCE);
      
      % sorted p-value plot
      axes(ax(1));
      hold off;
      plot([0,nstats+1],log10([.05,.05]),'c--');
      hold on;
      h.lmeCE = plot(1:nstats,log10(pval.aovLmeCE(idxPval)),'-x','Color',[0.5 0 0]);
      h.lme2CE = plot(1:nstats,log10(pval.aovLme2CE(idxPval)),'-o','Color',[0.5 0 0]);
      h.aovCE = plot(1:nstats,log10(pval.aovCEMain(idxPval)),'b-o');      
      h.kbttest = plot(1:nstats,log10(pval.kbttest(idxPval)),'k-+');
      %h.ttest = plot(1:nstats,pval.ttest(idxPval),'x-','Color',[.5,0,0]);
      h.ranksum = plot(1:nstats,log10(pval.ranksum(idxPval)),'k--');
      set(gca,'XTick',1:nstats,'XTickLabel',{},'XLim',[0,nstats+1]);
%       ylim = get(gca,'YLim');
%       ylim(1) = 0;
%       set(gca,'YLim',ylim);
      ylabel('log10 p-value','fontweight','bold');
      box off;
      legend([h.lmeCE h.lme2CE h.aovCE h.kbttest,h.ranksum],...
        {'LME' 'LME2' 'Anova2' 't-test' 'ranksum'},'location','best','interpreter','none');
      grpvarstr = BoxPlotPlot.groupVarStr(plotCfg.grpVars);
      grpvarstr2 = BoxPlotPlot.groupVarStr(plotCfg.grpVars2);
      titlestr1 = sprintf('Stats over %d Repeats (%d days), N=%d.',...
        numel(unique(cellstr(gRptC(tf01)))),numel(unique(cellstr(datesC(tf01)))),ntrials);
      titlestr2 = sprintf('Restrict1: %s {%s}. Restrict2: %s {%s}.',...
        plotCfg.grp.PrettyName,grpvarstr,plotCfg.grp2.PrettyName,grpvarstr2);
      titlestr = {titlestr1;titlestr2};
      title(titlestr,'interpreter','none','fontweight','bold','fontsize',12);
  
      % comparison plot
      % Stat i has bars in [i,i+1)
      % All stats normalized 
      % First 2 bars: Overall ctl, overall exp
      % Next 2*N bars: rpt1 ctl, exp1 ctl, rpt2 ctl, exp2 ctl
      
      % AL: in the following code, fixedeffect0~ctl and fixedeffect1~exp
      % aggregate data first so can make calls to errorbar all at once
      % (errorbar whisker length issue)
      ebYctl = nan(nstats,nrpt+1); % row i contains y-centers for errorbar for: [overallctl rpt1ctl rpt2ctl ...]
      ebSctl = nan(nstats,nrpt+1); 
      ebNctl = nan(nstats,nrpt+1); 
      ebYexp = nan(nstats,nrpt+1); % row i contains y-centers for errorbar for: [overallexp rpt1exp rpt2exp ...]
      ebSexp = nan(nstats,nrpt+1); 
      ebNexp = nan(nstats,nrpt+1); 
      for iStat = 1:nstats
        y = X(tf01,idxPval(iStat)); % all data for this stat
        ymu = nanmean(y);
        ysd = nanstd(y);
        
        yCtlrs = (X(tf0,idxPval(iStat)) - ymu)/ysd; % all ctl data for this stat
        yExprs = (X(tf1,idxPval(iStat)) - ymu)/ysd; % all exp data "
        ebYctl(iStat,1) = nanmean(yCtlrs);
        ebSctl(iStat,1) = nanstd(yCtlrs);
        ebNctl(iStat,1) = nnz(~isnan(yCtlrs));
        ebYexp(iStat,1) = nanmean(yExprs);
        ebSexp(iStat,1) = nanstd(yExprs);
        ebNexp(iStat,1) = nnz(~isnan(yExprs));
        
        for iRpt = 1:nrpt
          tfRpt = double(gRptC)==iRpt;
          yCtlRptrs = (X(tf0&tfRpt,idxPval(iStat)) - ymu)/ysd;
          yExpRptrs = (X(tf1&tfRpt,idxPval(iStat)) - ymu)/ysd;
          
          ebYctl(iStat,1+iRpt) = nanmean(yCtlRptrs);
          ebSctl(iStat,1+iRpt) = nanstd(yCtlRptrs);
          ebNctl(iStat,1+iRpt) = nnz(~isnan(yCtlRptrs));
          ebYexp(iStat,1+iRpt) = nanmean(yExpRptrs);
          ebSexp(iStat,1+iRpt) = nanstd(yExpRptrs);  
          ebNexp(iStat,1+iRpt) = nnz(~isnan(yExpRptrs));
        end        
      end
      
      axes(ax(2));
      hold on;      
      dx = 1/((nrpt+1)*3+2); % each of {Overall,Rpt1,Rpt2,...} gets 3 slots, first 2 have bars; add 2 for extra padding between stats
      LINEWIDTH = struct();
      LINEWIDTH.OVERALL = 2;
      LINEWIDTH.INDIV = 1;
      for iCol = 1:nrpt+1
        if iCol==1
          linewidth = LINEWIDTH.OVERALL;
        else
          linewidth = LINEWIDTH.INDIV;
        end
        hCtlDot = lclerrorbar((1:nstats)+(iCol-1)*3*dx,ebYctl(:,iCol)',ebSctl(:,iCol)'./sqrt(ebNctl(:,iCol)'),...
          linewidth,'o','Color',[.5,0,.5],'MarkerFaceColor',[.5,0,.5]);
        hExpDot = lclerrorbar((1:nstats)+(iCol-1)*3*dx+dx,ebYexp(:,iCol)',ebSexp(:,iCol)'./sqrt(ebNexp(:,iCol)'),...
          linewidth,'o','Color',[0,.5,.5],'MarkerFaceColor',[0,.5,.5]);
      end
      % make legend
      hLeg = nan(1,4);
      hLeg(1:2) = [hCtlDot(1) hExpDot(1)];
      hLeg(3) = plot(nan,nan,'k','LineWidth',LINEWIDTH.OVERALL);
      hLeg(4) = plot(nan,nan,'k','LineWidth',LINEWIDTH.INDIV);
      legend(hLeg,{gFixedCats{1} gFixedCats{2} 'overall' 'singleday'},'interpreter','none');

      set(gca,'XTick',1:nstats,'XTickLabel',Xnames(idxPval),'XLim',[0 nstats+1],'YTickMode','auto','YTickLabelMode','auto');
      rotateticklabel(gca);
      ylabel('z-score','fontweight','bold');
      box off;      
      
      linkaxes(ax,'x');
      
      descStr = 'Check crosstab, ctl-wash comparisons etc';
    end
    
  end
  
end

function [hDot,hLine] = lclerrorbar(x,y,s,linewidth,varargin)
assert(isequal(size(x),size(y),size(s)));
iClr = find(strcmp(varargin,'Color'));
N = numel(s);
hDot = nan(N,1);
hLine = nan(N,1);
for i = 1:N
  hDot(i) = plot(x(i),y(i),varargin{:});
  hLine(i) = plot([x(i) x(i)],[y(i)-s(i) y(i)+s(i)],'Color',varargin{iClr+1},'LineWidth',linewidth);
end
end
