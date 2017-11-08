classdef RepeatSuccess < MultiStatCompareBase & RepeatPlotBase
  
  properties (Constant)
    MEDGRAY = [0.5 0.5 0.5];
    MEDBLUE = [0 0 0.75];
    LINESRED = [0.85 0.325 0.098];
    LINESORG = [0.929 0.694 0.125];
    LINESPRP = [0.635 0.078 0.184];
    LINESGRN = [0.466 0.674 0.188];
  end
  
  properties
    Name = 'SuccessRate Analysis';
    AllowsMultiStats = false;
    UsesAuxVar = false;
    
    UsesGroupVar2 = true;
    GroupVarLabel = 'Restrict Var';
    GroupVar2Label = 'Restrict Var';
    
    AutoReplotAfterStatChange = false;
    AutoReplotAfterGroupChange = false;
  end
  
  methods        
            
    function obj = RepeatSuccess(hObj)
      assert(isa(hObj,'HandleObj'));
      obj = obj@RepeatPlotBase(hObj);
    end
    
    function descStr = doPlot(obj,ax,data,bstat,~,plotCfg,~)
      [ok,~,tfKeep,datesAllC,gRptC,gCtlExpWashC] = obj.collectData(data);
      if ~ok
        descStr = obj.FailedRepeatSelection(ax);
        return;
      end
      assert(isequal(size(tfKeep),[size(bstat,1) 1]));
      bstat = bstat(tfKeep,:);
      if numel(unique(bstat))~=2
        descStr = obj.Failed(ax,'Selected data is not a dichotomous variable, or does not contain two distinct values.');
        return;
      end
      
      newax = obj.setupAxes(ax);
      
      bstatnames = {plotCfg.stat.Name}';
      descStr = obj.doPlotCore(newax(1),newax(2),bstat,bstatnames,gRptC,gCtlExpWashC,datesAllC,plotCfg);     
    end
    
  end
  
  methods (Static)
    
    function descStr = doPlotCore(axTop,axBot,X,Xnames,gRptC,gCEWC,datesC,plotCfg)
      
      [ntrials,nstats] = size(X);
      assert(iscellstr(Xnames) && isvector(Xnames) && numel(Xnames)==nstats);
      assert(iscategorical(gRptC) && isvector(gRptC) && numel(gRptC)==ntrials);
      assert(iscategorical(gCEWC) && isvector(gCEWC) && numel(gCEWC)==ntrials);
      assert(iscategorical(datesC) && isvector(datesC) && numel(datesC)==ntrials);
      nrpt = numel(unique(cellstr(gRptC)));
      assert(nrpt==numel(categories(gRptC))); % Repeat.anlsDichotomous works off of numel(categories(gRptC))
       
      CTLEXPWSH = {'ctl' 'exp' 'wsh'};
      CTLEXP = CTLEXPWSH(1:2);
      CTLWSH = CTLEXPWSH([1 3]);
      SPREXP = {'spr' 'exp'}; % "supercontrol" vs exp
      ZEROONE = {'0' '1'};
      
      % Ctl-Exp analysis
      tfCE = ismember(gCEWC,CTLEXP);
      gCEWC_CE = setcats(gCEWC(tfCE),CTLEXP);
      [lctCE,rctCE] = Repeat.anlsDichotomous(X(tfCE,:),gRptC(tfCE),gCEWC_CE);
      
      % Ctl-Wsh analysis
      tfCW = ismember(gCEWC,CTLWSH);
      gCEWC_CW = setcats(gCEWC(tfCW),CTLWSH);
      [lctCW,rctCW] = Repeat.anlsDichotomous(X(tfCW,:),gRptC(tfCW),gCEWC_CW);
      
      % Spr-Exp analysis
      gSE = cell(size(gCEWC));
      gSE(gCEWC=='ctl') = {'spr'};
      gSE(gCEWC=='wsh') = {'spr'};
      gSE(gCEWC=='exp') = {'exp'};
      assert(~any(cellfun(@isempty,gSE)));
      gSE = categorical(gSE,SPREXP);
      [lctSE,rctSE] = Repeat.anlsDichotomous(X,gRptC,gSE);
      
      ALL = struct();
      [ALL.tbl,ALL.chi2,ALL.chi2_p,ALL.lbls] = crosstab(X,gCEWC,gRptC);
      ALL.tblFlat = sum(ALL.tbl,3);
      [~,ALL.chi2_flat,ALL.chi2_p_flat] = crosstab(X,gCEWC); % for chi2 p-val on 2x3 contingency table (0/1 x C/E/W)
      
      % Checks on ALL, LCT, RCT
      assert(isequal(size(lctCE.m),[2 2]));
      assert(isequal(lctCE.mLbl,[CTLEXP(:) ZEROONE(:)]));
      assert(isequal(lctCE.m,ALL.tblFlat(:,1:2)'));
      
      [tmp1,tmp2,tmp3] = size(rctCE.m);
      assert(isequal([tmp1 tmp2 tmp3],[2 2 nrpt]));
      assert(isequal(rctCE.mLbl{1}',CTLEXP) && isequal(rctCE.mLbl{2}',ZEROONE)); % rctCE.mLbl{3} checked below
      assert(isequal(permute(ALL.tbl(:,1:2,:),[2 1 3]),rctCE.m));
      assert(numel(rctCE.stats_fisher_rpts)==nrpt);
      
      assert(isequal(lctCW.mLbl,[CTLWSH(:) ZEROONE(:)]));
      assert(isequal(lctCW.m,ALL.tblFlat(:,[1 3])'));
      
      assert(isequal(rctCW.mLbl{1}',CTLWSH) && isequal(rctCW.mLbl{2}',ZEROONE));
      assert(isequal(permute(ALL.tbl(:,[1 3],:),[2 1 3]),rctCW.m));
      %assert(numel(rctCW.stats_fisher_rpts)==nrpt);
      
      assert(isequal(lctSE.mLbl,[SPREXP(:) ZEROONE(:)]));
      tmp = [];
      tmp(:,1) = sum(ALL.tblFlat(:,[1 3]),2);
      tmp(:,2) = ALL.tblFlat(:,2);
      assert(isequal(lctSE.m,tmp'));
      
      assert(isequal(rctSE.mLbl{1}',SPREXP) && isequal(rctSE.mLbl{2}',ZEROONE));
      tmp = [];
      tmp(:,1,:) = sum(ALL.tbl(:,[1 3],:),2);
      tmp(:,2,:) = ALL.tbl(:,2,:);
      assert(isequal(permute(tmp,[2 1 3]),rctSE.m));
      %assert(numel(rctSE.stats_fisher_rpts)==nrpt);

      assert(isequal(ALL.lbls(1:2,1)',ZEROONE) && isequal(ALL.lbls(1:3,2)',CTLEXPWSH));
      tmp = ALL.lbls(:,3);
      tmp = tmp(~cellfun(@isempty,tmp));
      assert(isequal(tmp,rctCE.mLbl{3},rctCW.mLbl{3},rctSE.mLbl{3}));
      
      %%% Visualize      
      vizM = nan(nrpt+1,3); % age successful. Rows: overall, rpt1,rpt2,...rptK. Cols: ctl,exp,wsh
      vizMCnt = nan(nrpt+1,3); % trial counts. Rows/Cols same as vizM.
      vizOR_CE = nan(nrpt+1,3); % Ctl-Exp OR. Rows: etc. Cols: Log-OR,Log-OR-lowerCI,Log-OR-upperCI
      vizOR_CW = nan(nrpt+1,3);
      vizOR_SE = nan(nrpt+1,3); 
      
      % Overall
      vizM(1,:) = ALL.tblFlat(2,:)./sum(ALL.tblFlat,1);
      vizMCnt(1,:) = sum(ALL.tblFlat,1);
      vizOR_CE(1,:) = hlpOR(lctCE.stats_fisher);
      vizOR_CW(1,:) = hlpOR(lctCW.stats_fisher);
      vizOR_SE(1,:) = hlpOR(lctSE.stats_fisher);
      
      % Rpts
      for k = 1:nrpt
        vizM(k+1,:) = ALL.tbl(2,:,k)./sum(ALL.tbl(:,:,k),1);
        vizMCnt(k+1,:) = sum(ALL.tbl(:,:,k),1);
        vizOR_CE(k+1,:) = hlpOR(rctCE.stats_fisher_rpts{k});
        vizOR_CW(k+1,:) = hlpOR(rctCW.stats_fisher_rpts{k});
        vizOR_SE(k+1,:) = hlpOR(rctSE.stats_fisher_rpts{k});
      end
      
      vizMRowLbls = [{'Total'}; rctCE.mLbl{3}];
      % Append counts of ctl/exp to vizMRowLbls
      for k = 1:numel(vizMRowLbls)
        vizMRowLbls{k} = sprintf('%s (%d/%d/%d)',vizMRowLbls{k},vizMCnt(k,1),vizMCnt(k,2),vizMCnt(k,3));
      end

      barX = [0.9 2:size(vizM,1)];
      hBar = bar(axBot,barX,vizM,'grouped');
      set(hBar(1),'FaceColor',RepeatSuccess.LINESGRN);
      set(hBar(2),'FaceColor',RepeatSuccess.MEDBLUE);
      set(hBar(3),'FaceColor',RepeatSuccess.MEDGRAY);

      set(axBot,'XTick',barX,'XTickLabel',vizMRowLbls);
      ylim(axBot,[0 1.05]);
      set(axBot,'YTick',0:.2:1,'YTickLabelMode','auto');
      ylabel(axBot,'Proportion Grab Successful');
      xlabel(axBot,'Repeat (nos. C/E/W)');
      grid(axBot,'on');
      legend(axBot,hBar,CTLEXPWSH,'Location','Best');
      
      % OR
      orX = barX;
      hOR = bar(axTop,orX,[vizOR_CW(:,1) vizOR_CE(:,1) vizOR_SE(:,1)]);
      hold(axTop,'on');
      set(hOR(1),'FaceColor',RepeatSuccess.LINESORG,'EdgeColor','none');
      set(hOR(2),'FaceColor',RepeatSuccess.LINESRED,'EdgeColor','none');
      set(hOR(3),'FaceColor',RepeatSuccess.LINESPRP,'EdgeColor','none');
      hORErr_CW = errorbar(axTop,orX+hOR(1).XOffset,vizOR_CW(:,1),vizOR_CW(:,1)-vizOR_CW(:,2),vizOR_CW(:,3)-vizOR_CW(:,1));
      hORErr_CE = errorbar(axTop,orX+hOR(2).XOffset,vizOR_CE(:,1),vizOR_CE(:,1)-vizOR_CE(:,2),vizOR_CE(:,3)-vizOR_CE(:,1));
      hORErr_SE = errorbar(axTop,orX+hOR(3).XOffset,vizOR_SE(:,1),vizOR_SE(:,1)-vizOR_SE(:,2),vizOR_SE(:,3)-vizOR_SE(:,1));
      set(hORErr_CW,'LineStyle','none','LineWidth',2);
      tmpclr = get(hORErr_CW,'Color');
      set(hORErr_CE,'LineStyle','none','LineWidth',2,'Color',tmpclr);
      set(hORErr_SE,'LineStyle','none','LineWidth',2,'Color',tmpclr);
      
      xlim(axTop,xlim(axBot));
      set(axTop,'XTick',barX,'XTickLabel',[]);
      grid(axTop,'on');
      vizORGoodVals = [vizOR_CW(~isnan(vizOR_CW) & ~isinf(vizOR_CW)); ...
                       vizOR_CE(~isnan(vizOR_CE) & ~isinf(vizOR_CE)); ...
                       vizOR_SE(~isnan(vizOR_SE) & ~isinf(vizOR_SE))];
      if ~isempty(vizORGoodVals)
        tmp = max(abs(vizORGoodVals(:)));
        ylim(axTop,[-1.1*tmp 1.1*tmp]);
      end
      set(axTop,'YTickMode','auto','YTickLabelMode','auto');
      ylabel(axTop,'Log Odds Ratio');
      legend(axTop,hOR,{'ctl-wsh' 'ctl-exp' 'spr-exp'},'Location','Best');
      
      grpvarstr = BoxPlotPlot.groupVarStr(plotCfg.grpVars);
      grpvarstr2 = BoxPlotPlot.groupVarStr(plotCfg.grpVars2);
      titlestr1 = sprintf('Success Rate Analysis: %s, %d Repeats (%d days), N=%d.',...
        Xnames{1},nrpt,numel(unique(cellstr(datesC))),ntrials);
      titlestr2 = sprintf('Restrict1: %s {%s}. Restrict2: %s {%s}.',...
        plotCfg.grp.PrettyName,grpvarstr,plotCfg.grp2.PrettyName,grpvarstr2);
      titlestr = {titlestr1;titlestr2};
      title(titlestr,'interpreter','none','fontweight','bold','fontsize',12);      
      
      descStr = {...
        sprintf('p_chi2, lumped, 3-way C-E-W: %.3g',ALL.chi2_p_flat);
        sprintf('p_chi2, lumped (C-W,C-E,S-E):');
        sprintf(' %7.3g %7.3g %7.3g',lctCW.p_chi2,lctCE.p_chi2,lctSE.p_chi2); ...
        sprintf('p_fish, lumped (C-W,C-E,S-E):');
        sprintf(' %7.3g %7.3g %7.3g',lctCW.p_fisher,lctCE.p_fisher,lctSE.p_fisher); ...
        sprintf('p_cmh (C-W,C-E,S-E):');
        sprintf(' %7.3g %7.3g %7.3g',rctCW.p_cmh,rctCE.p_cmh,rctSE.p_cmh); ...
        sprintf('fisher OddsRatio w/95%% CI, lumped:'); ...
        sprintf(' C-W: %7.3g [%7.3g %7.3g]',...
          lctCW.stats_fisher.OddsRatio,...
          lctCW.stats_fisher.ConfidenceInterval(1),...
          lctCW.stats_fisher.ConfidenceInterval(2));
        sprintf(' C-E: %7.3g [%7.3g %7.3g]',...
          lctCE.stats_fisher.OddsRatio,...
          lctCE.stats_fisher.ConfidenceInterval(1),...
          lctCE.stats_fisher.ConfidenceInterval(2));
        sprintf(' S-E: %7.3g [%7.3g %7.3g]',...
          lctSE.stats_fisher.OddsRatio,...
          lctSE.stats_fisher.ConfidenceInterval(1),...
          lctSE.stats_fisher.ConfidenceInterval(2));
        sprintf('C-E by rpts: p_fisher OddsRatio [OR CIs]');
        };      
      for i = 1:numel(rctCE.p_fisher_rpts)
        descStr{end+1,1} = ...
          sprintf(' RPT%d: %7.3g %7.3g [%7.3g %7.3g]',i,...
          rctCE.p_fisher_rpts(i),...
          rctCE.stats_fisher_rpts{i}.OddsRatio,...
          rctCE.stats_fisher_rpts{i}.ConfidenceInterval(1),...
          rctCE.stats_fisher_rpts{i}.ConfidenceInterval(2)); %#ok<AGROW>
      end
    end
    
    function descStr = Failed(ax,msg)
      xl = xlim(ax);
      yl = ylim(ax);
      hTxt = text(mean(xl),mean(yl),msg);
      set(hTxt,'horizontalalignment','center');
      
      descStr = {msg};
    end
    
  end
  
end

function orvec = hlpOR(stats)
or = hlpLogsafe(stats.OddsRatio);
orlo = hlpLogsafe(stats.ConfidenceInterval(1));
orhi = hlpLogsafe(stats.ConfidenceInterval(2));
orvec = [or orlo orhi];
end

function y = hlpLogsafe(x)
% negative values of x were generating complex values
if x<0
  y = nan;
else
  y = log(x);
end
end

% function [hDot,hLine] = lclerrorbar(x,y,s,linewidth,varargin)
% assert(isequal(size(x),size(y),size(s)));
% iClr = find(strcmp(varargin,'Color'));
% N = numel(s);
% hDot = nan(N,1);
% hLine = nan(N,1);
% for i = 1:N
%   hDot(i) = plot(x(i),y(i),varargin{:});
%   hLine(i) = plot([x(i) x(i)],[y(i)-s(i) y(i)+s(i)],'Color',varargin{iClr+1},'LineWidth',linewidth);
% end
% end
