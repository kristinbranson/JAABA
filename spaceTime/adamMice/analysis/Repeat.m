classdef Repeat
  % A Repeat is a set of three days, Ctl-Exp-Wash.
  
  methods (Static)

    function [pAnovaInt,pAnovaFxdMain,pLmeFxd,pLme2Fxd,pAnova1Fxd,lmeStats] = ...
        anls(y,gRpt,gFixed)
      % y: data vector
      % gRpt: cellstr, grouping vector for Repeat
      % gFixed: cellstr, binary grouping vector for fixed effect. Typically 
      %   gFixed labels eg Ctl vs Exp or Ctl vs Wash, but in general can 
      %   be any dichotomous quantity that will be treated as a fixed 
      %   effect in a mixed model with gRpt.
      
      assert(isvector(y));
      assert(isvector(gRpt)&&iscellstr(gRpt)&&numel(gRpt)==numel(y));
      assert(isvector(gFixed)&&iscellstr(gFixed)&&numel(gFixed)==numel(y)); %&&numel(unique(gFixed))<=2);
      
      % Sanity check: 1-way anova, CE fixed effect only
      pAnova1Fxd = anova1(y,gFixed,'off');
      
      % 2-way anova with interactions, to determine strength of Rpt/CE interaction      
      [aov.int.p,aov.int.t,aov.int.stats,aov.int.terms] = anovan(y,{gRpt,gFixed},...
        'display','off','random',1,'varnames',{'rpt' 'CE'},'model','interaction');

      % 2-way anova, main effects only
      [aov.main.p,aov.main.t,aov.main.stats,aov.main.terms] = anovan(y,{gRpt,gFixed},...
        'display','off','random',1,'varnames',{'rpt' 'CE'},'model','linear');

      % LME, main effects 
      tbl = table(y,gRpt,gFixed,'variablenames',{'y' 'rpt' 'ce'});
      lme = fitlme(tbl,'y ~ ce + (1|rpt)');
      aov.lme = anova(lme);            
      % LME, random slope for ce
      lme2 = fitlme(tbl,'y ~ ce + (1+ce|rpt)');
      aov.lme2 = anova(lme2);
      
      %%% LME checks/diagnostics %%%
      assert(lme.VariableInfo('rpt',:).IsCategorical);
      assert(lme.VariableInfo('ce',:).IsCategorical);
      assert(lme2.VariableInfo('rpt',:).IsCategorical);
      assert(lme2.VariableInfo('ce',:).IsCategorical);
      
      [~,~,stats] = lme.covarianceParameters();
      stats = stats{1};
      assert(strcmp(stats.Group,'rpt'));
      assert(strcmp(stats.Name1,'(Intercept)'));
      assert(strcmp(stats.Name2,'(Intercept)'));
      lmeStats = struct();      
      lmeStats.lme.RptIntercept_SDEst = [stats.Estimate stats.Lower stats.Upper];
      
      [~,~,stats2] = lme2.covarianceParameters();
      stats2 = stats2{1};
      assert(strcmp(strtrim(stats2.Group(1,:)),'rpt'));
      assert(strcmp(stats2.Name1{1},'(Intercept)'));
      assert(strcmp(stats2.Name2{1},'(Intercept)'));
      lmeStats.lme2.RptIntercept_SDEst = [stats2.Estimate(1) stats2.Lower(1) stats2.Upper(1)];
      nFixed = numel(unique(gFixed));
      if nFixed==2
        % Following will not hold if nFixed>2; (will be multiple stds)
        assert(strcmp(strtrim(stats2.Group(3,:)),'rpt'));
        assert(strncmp(stats2.Name1{3},'ce_',3)); % could be ce_exp or ce_ctl etc
        assert(strncmp(stats2.Name2{3},'ce_',3));
        lmeStats.lme2.RptCE_SDEst = [stats2.Estimate(3) stats2.Lower(3) stats2.Upper(3)];
      end
      
      % Compare LMEs
      lmeCmp = compare(lme,lme2,'CheckNesting',true);
      assert(isequal(cellstr(lmeCmp.Model),{'lme';'lme2'}));
      lmeStats.p_lmecmp = lmeCmp.pValue(2);
            
      assert(isequal(aov.int.terms(end,:),[1 1]));
      pAnovaInt = aov.int.p(end);
      
      assert(isequal(aov.main.terms(2,:),[0 1]));
      pAnovaFxdMain = aov.main.p(2);
      
      assert(strcmp(aov.lme.Term{2},'ce'));
      pLmeFxd = aov.lme.pValue(2);
      assert(strcmp(aov.lme2.Term{2},'ce'));
      pLme2Fxd = aov.lme2.pValue(2);
    end
    
    function [LCT,RCT] = anlsDichotomous(y01,gRpt,gCE)
      % Grab successrate analysis
      % y01: dichotomous outcome variable (0/1)
      % gRpt: Repeat grouping vec
      % gCE: Ctl-Exp grouping vec
      
      assert(isvector(y01) && isnumeric(y01) && all(y01==0 | y01==1));
      assert(isvector(gRpt) && numel(gRpt)==numel(y01) && iscategorical(gRpt));
      assert(isvector(gCE) && numel(gCE)==numel(y01) && iscategorical(gCE));
      
      N = numel(y01);
      gCEcats = categories(gCE);
      gRptcats = categories(gRpt);
      nCEcats = numel(gCEcats);
      nRptcats = numel(gRptcats);
      assert(nCEcats==2,'Expect 2 categories for Ctl-Exp');
            
      % lumped crosstab
      LCT = struct();
      [LCT.m,LCT.chi2,LCT.p_chi2,LCT.mLbl] = crosstab(gCE,y01);
      assert(sum(LCT.m(:))==N);
      % permute rows/cols of LCT.m so it matches RCT.m below
      [tf1,loc1] = ismember(gCEcats,LCT.mLbl(:,1));
      [tf2,loc2] = ismember({'0' '1'},LCT.mLbl(:,2));
      assert(all(tf1) && all(tf2));
      LCT.mLbl(:,1) = LCT.mLbl(loc1,1);
      LCT.mLbl(:,2) = LCT.mLbl(loc2,2);
      LCT.m = LCT.m(loc1,loc2);

      % Fischer's exact (lumped)
      tmp = which('fishertest');
      tfFisherAvail = ~isempty(tmp);
      if ~tfFisherAvail
        warning('No fishertest() function available. Do you have Stats+MATLAB R2014b?');
      end
      if tfFisherAvail
        [~,LCT.p_fisher,LCT.stats_fisher] = fishertest(LCT.m);
      else
        LCT.p_fisher = nan;
        LCT.stats_fisher = struct('OddsRatio',nan,'ConfidenceInterval',[nan nan]);
      end
      
      % repeated crosstabs
      RCT = struct();
      RCT.m = nan(nCEcats,2,nRptcats);
      for iRpt = 1:nRptcats
        for iCE = 1:nCEcats
          for iY = 1:2
            RCT.m(iCE,iY,iRpt) = nnz(gRpt==gRptcats(iRpt) & gCE==gCEcats(iCE) & y01==iY-1);
          end
        end
      end
      RCT.mLbl{1} = gCEcats;
      RCT.mLbl{2} = {'0';'1'};
      RCT.mLbl{3} = gRptcats;
      assert(isequal(sum(RCT.m,3),LCT.m));
            
      % CMH
      [RCT.chiMH,RCT.p_cmh] = CMH(RCT.m);
      RCT.p_fisher_rpts = nan(nRptcats,1);
      RCT.stats_fisher_rpts = cell(nRptcats,1);
      for i = 1:nRptcats
        if tfFisherAvail
          [~,RCT.p_fisher_rpts(i),RCT.stats_fisher_rpts{i}] = fishertest(RCT.m(:,:,i));
        else
          RCT.p_fisher_rpts(i) = nan;
          RCT.stats_fisher_rpts{i} = struct('OddsRatio',nan,'ConfidenceInterval',[nan nan]);
        end
      end      
    end     
    
    function [repeats,ok] = SelectRepeats(datesUn)
      % [repeats,ok] = SelectRepeats(datesUn)
      %
      % repeats: struct array with fields ctl, exp, wsh. Field values are
      % elements of datesUn. Field values are not guaranteed to be distinct
      % values of datesUn. The .wsh field is optional and may be [] for
      % any/all elements of repeats.
      % ok: scalar logical. If true, repeats are good and explicitly
      % verified by user.
      
      REPEATSEMPTY = struct('ctl',cell(0,1),'exp',[],'wsh',[]);
      repeats = REPEATSEMPTY;
      
      while true
        [idx,ok] = listdlg('ListString',datesUn,'PromptString','Select Repeat',...
          'CancelString','Done with Repeats','listsize',[320 300]);
        if ~ok
          % Done or killed
          break;
        end
        if numel(idx)~=2 && numel(idx)~=3
          uiwait(warndlg('Please select 2 or 3 dates for Control, Experiment, and (optional) Wash.',...
            'Invalid Repeat','modal'));
          continue;
        end
        
        incWash = numel(idx)==3;
        tfSuccessfulRpt = false;
        while true
          if incWash
            % exit from loop:
            % User cancels, tfSuccessfulTrip remains false
            % User says the ordering is correct, tfSuccessfulTrip is true
            str = sprintf('Is this correct?\nCtl: %s\nExp: %s\nWsh: %s\n',datesUn{idx(1)},datesUn{idx(2)},datesUn{idx(3)});
            btn = questdlg(str,'Order repeat','Yes','No, rotate up','No, switch last two','Yes');
            if isempty(btn)
              break;
            end
            switch btn
              case 'Yes', tfSuccessfulRpt = true; break;
              case 'No, rotate up', idx = idx([2 3 1]);
              case 'No, switch last two', idx = idx([1 3 2]);
            end
          else
            str = sprintf('Is this correct?\nCtl: %s\nExp: %s\n',datesUn{idx(1)},datesUn{idx(2)});
            btn = questdlg(str,'Order repeat','Yes','No, switch','Yes');
            if isempty(btn)
              break;
            end
            switch btn
              case 'Yes', tfSuccessfulRpt = true; break;
              case 'No, switch', idx = idx([2 1]);
            end
          end
        end
        
        if tfSuccessfulRpt
          repeats(end+1).ctl = datesUn{idx(1)}; %#ok<AGROW>
          repeats(end).exp = datesUn{idx(2)};
          if incWash
            repeats(end).wsh = datesUn{idx(3)};
          end
        end
      end
      
      % repeats could be empty, or user could have canceled etc; need to
      % verify
      ok = Repeat.uiVerify(repeats);
    end
    
    function ok = uiVerify(repeats,varargin)
      prestr = myparse(varargin,...
        'prestr','Are these repeats correct?');
      
      if isempty(repeats)
        ok = false;
        return;
      end
      
      str = arrayfun(@(x)sprintf('%s -> %s -> %s',x.ctl,x.exp,x.wsh),repeats,'uni',0);
      % check washes
      washes = {repeats.wsh};
      tfEmpty = cellfun(@isempty,washes);
      if all(tfEmpty) || ~any(tfEmpty)
        str = [{prestr};str(:)];
      else
        warnstr = 'WARNING: some repeats have washes specified and others do not.';
        str = [{prestr;warnstr};str(:)];
      end        
      
      btn = questdlg(str,'Confirm repeats','Yes','No','Yes');
      if isempty(btn) || strcmp(btn,'No')
        ok = false;
      else
        ok = true;
      end      
    end
    
  end
  
end
