classdef RepeatBoxPlot < OlyDat.Plot & RepeatPlotBase
  
  properties
    Name = 'Repeat Box Plot';
    AllowsMultiStats = false;
    UsesAuxVar = false;
    
    UsesGroupVar2 = true;
    GroupVarLabel = 'Fixed Effect';
    GroupVar2Label = 'Restrict Var';
    
    AutoReplotAfterStatChange = true;
    AutoReplotAfterGroupChange = false;
  end
  
  methods    
        
    function obj = RepeatBoxPlot(hObj)
      assert(isa(hObj,'HandleObj'));
      obj = obj@RepeatPlotBase(hObj);
    end
    
    function descStr = doPlot(obj,ax,data,bstat,grp,plotCfg,~)
      
      [ok,Nrpt,tfKeep,datesAllC,gRptC,gCtlExpWashC] = obj.collectData(data);
      if ~ok
        descStr = obj.FailedRepeatSelection(ax);
        return;
      end
      
      assert(isequal(size(tfKeep),[size(bstat,1) 1]));
      bstat = bstat(tfKeep,:);
        
      gFixed = plotCfg.grp.Name;
      gFixedVals = grp{1};
      assert(isequal(size(tfKeep),size(gFixedVals)));
      gFixedVals = gFixedVals(tfKeep);
      gFixedVals = ExpPP.ensureStringCategorical(gFixedVals,gFixed);

      if isempty(gFixed) 
        % <no group>; use Default C/E/W fixed effect
      
        boxplot(ax,bstat,{gRptC gCtlExpWashC datesAllC},...
          'colorgroup',gCtlExpWashC,...
          'colors',[0 0.5 0;0 0 1;0.5 0.5 0.5],'factorgap',[10 0 0],...
          'factorseparator',1,'boxstyle','outline','widths',.8);      
        grid on;

        % compute stats
        assert(isequal(size(bstat),size(gRptC),size(gCtlExpWashC)));
        tfCE = ismember(gCtlExpWashC,{'ctl' 'exp'});
        tfCW = ismember(gCtlExpWashC,{'ctl' 'wsh'});

        [anlsCE.pAnovaInt,anlsCE.pAnovaCEMain,anlsCE.pLmeCE,anlsCE.pLme2CE,~,anlsCE.lmeStats] = ...
          Repeat.anls(bstat(tfCE),cellstr(gRptC(tfCE)),cellstr(gCtlExpWashC(tfCE)));
        [anlsCW.pAnovaInt,anlsCW.pAnovaCEMain,anlsCW.pLmeCE,anlsCW.pLme2CE,~,anlsCW.lmeStats] = ...
          Repeat.anls(bstat(tfCW),cellstr(gRptC(tfCW)),cellstr(gCtlExpWashC(tfCW)));
        
        descStr = {...
          'Ctl-Wash pInteract/pMain/pLME/pLME2:';...
          sprintf('%0.3g %0.3g %0.3g %0.3g',anlsCW.pAnovaInt,anlsCW.pAnovaCEMain,anlsCW.pLmeCE,anlsCW.pLme2CE); ...
          'Ctl-Exp pInteract/pMain/pLME/pLME2:';...
          sprintf('%0.3g %0.3g %0.3g %0.3g',anlsCE.pAnovaInt,anlsCE.pAnovaCEMain,anlsCE.pLmeCE,anlsCE.pLme2CE); ...
          '';
          'LME notes:';
          ' C-E LME vs LME2 compare:';
          sprintf(' p=%.3g',anlsCE.lmeStats.p_lmecmp);
          ' C-E, LME Rpt intrcpt [std CIlo CIhi]:';
          sprintf(' %s',num2str(anlsCE.lmeStats.lme.RptIntercept_SDEst));
          ' C-E, LME2 Rpt intrcpt [std CIlo CIhi]:';
          sprintf(' %s',num2str(anlsCE.lmeStats.lme2.RptIntercept_SDEst));
          ' C-E, LME2 Rpt CE [std CIlo CIhi]:';
          sprintf(' %s',num2str(anlsCE.lmeStats.lme2.RptCE_SDEst));
          };
        
        fixedEffectStr = 'Ctl vs Exp';
        fixedEffectXLabel = '';
      else
        % User specified group for fixed effect
        
        gFixedValsUn = unique(gFixedVals);
%         if numel(gFixedValsUn)>2
%           errstr = sprintf('Fixed effect ''%s'' is not a dichotomous variable.',gFixed);
%           OlyDat.Plot.messageAxisCenter(ax,errstr);
%           descStr = '';
%           return;
%         end
        
        boxplot(ax,bstat,{gRptC gFixedVals},...
          'colorgroup',gFixedVals,...
          'colors',[0 0.5 0;0 0 1],'factorgap',[10 0],...
          'factorseparator',1,'boxstyle','outline','widths',.8);
        grid on;

        % compute stats
        assert(isequal(size(bstat),size(gRptC),size(gFixedVals)));
        [anls.pAnovaInt,anls.pAnovaFxdMain,anls.pLme,anls.pLme2,~,anls.lmeStats] = ... 
          Repeat.anls(bstat,cellstr(gRptC),gFixedVals);

        fixedEffectStr = sprintf('%s vs ',gFixedValsUn{:});
        fixedEffectStr = fixedEffectStr(1:end-4);
        fixedEffectXLabel = gFixed;

        descStr = {...
          fixedEffectStr; ...
          sprintf('pInteract/pMain/pLME/pLME2:'); ...
          sprintf('%0.3g %0.3g %0.3g %0.3g',...
            anls.pAnovaInt,anls.pAnovaFxdMain,anls.pLme,anls.pLme2); ...
          'LME notes:'; ...
          ' LME vs LME2 compare:'; ...
          sprintf(' p=%.3g',anls.lmeStats.p_lmecmp); ...
          ' LME Rpt intrcpt [std CIlo CIhi]:'; ...
          sprintf(' %s',num2str(anls.lmeStats.lme.RptIntercept_SDEst)); ...
          ' LME2 Rpt intrcpt [std CIlo CIhi]:'; ...
          sprintf(' %s',num2str(anls.lmeStats.lme2.RptIntercept_SDEst)); };
        nFixedValsUn = numel(gFixedValsUn);
        if nFixedValsUn>2
          descStr{end+1,1} = '';
          descStr{end+1,1} = sprintf('NOTE: %d-way Fixed effect comparison',nFixedValsUn);
        end
      end
        
      %grpvarstr = BoxPlotPlot.groupVarStr(plotCfg.grpVars);
      grpvarstr2 = BoxPlotPlot.groupVarStr(plotCfg.grpVars2);
      titlestr1 = sprintf('%s over %d Repeats, N=%d.',...
        plotCfg.stat.PrettyName,Nrpt,numel(bstat));
      titlestr2 = sprintf('Fixed effect: %s.',fixedEffectStr);
      titlestr3 = sprintf('Restrict: %s {%s}.',...
        plotCfg.grp2.PrettyName,grpvarstr2);
      titlestr = {titlestr1;titlestr2;titlestr3};      
      title(ax,titlestr,'interpreter','none','fontsize',10);
      ylabel(ax,plotCfg.stat.PrettyName,'interpreter','none','fontweight','normal','fontsize',ComparisonPlot.AXISTITLE_FONTSIZE);
      xlabelstr = {'Repeat number';fixedEffectXLabel};
      xlabel(ax,xlabelstr,'interpreter','none','fontweight','bold','fontsize',10);
      
      switch version('-release')
        case '2014b'
          hTxt = findobj(ax,'Type','Text');
          set(hTxt,'Rotation',20);
      end
    end
    
  end
  
end
