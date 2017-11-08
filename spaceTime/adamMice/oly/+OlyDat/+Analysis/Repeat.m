classdef Repeat
   
    properties (Constant)
        LINE_FLD = 'line_name';
        CONTROL_LINE = 'pBDPGAL4U';
    end
    
    methods (Static) % Applications
        
        % Distros + repeat superplot
        function mtgFig(hFig,rd,fld)
            clf(hFig);
            figure(hFig);
            
            ginfo.g = {rd.gmrVsCtrl}';
            ginfo.gName = 'gVc';

            hAx = subplot(3,1,1);
            y = cat(1,rd.(fld));
            assert(isequal(size(y),[numel(rd) 1]));
            BoxData.Variation.plotDensitiesCore(hAx,y,ginfo,...
                'showStats',true);
            axis1 = axis;
            set(hAx,'Box','off');
            set(hAx,'FontSize',30);
            set(hAx,'XTickLabel',[]);
            
            hAx = subplot(3,1,2:3);
            line2Idx = OlyDat.Analysis.Repeat.line2DataIdx(rd);
            alllinenames = {rd.line_name}';
            tfRpted = cellfun(@(x)numel(line2Idx(x))>1,alllinenames);
            
            rdRpted = rd(tfRpted);
            line2IdxRpted = OlyDat.Analysis.Repeat.line2DataIdx(rdRpted); % prob unnnecessary
            grptedInfo.g = {rdRpted.gmrVsCtrl}';
            grptedInfo.gName = 'gVc';
            
            BoxData.Variation.experimentLineSuperPlot(hAx,rdRpted,...
                fld,line2IdxRpted,grptedInfo);
            axis2 = axis;
            axis2(1:2) = axis1(1:2);
            axis(axis2);
            set(hAx,'FontSize',30);
            ylabel('Rank');
            xlabel('Stat');
            
            % append globalLineSD to title
            [~,globalLineSD] = OlyDat.Analysis.Repeat.repeatabilityAllLines(rd,fld);
            hT = get(hAx,'Title');
            str = sprintf('%s globalLineSD=%.4g.',get(hT,'String'),globalLineSD);
            set(hT,'String',str);
        end
        
        function cmpRepeatability(a0rd,a1rd,a0fld,a1fld,a0RptInfo,a1RptInfo,varargin)
            % cmpRepeatability(a0rd,a1rd,a0RptInfo,a1RptInfo,varargin)
            % "Before/after" comparison of repeatability.
            % a0rd: before dataset
            % a1rd: after dataset
            % a0RptInfo: before repeatInfo, struct with fields 'dm',
            % 'sdAll','sdCtl', (outputs of repeatabilityAllLines)
            % a1RptInfo: etc
            % varargin: optional PV pairs:
            %   * 'doPlot': if true, make a plot. Default=true.

            assert(isstruct(a0rd));
            assert(isstruct(a1rd));
            assert(ischar(a0fld)&&isfield(a0rd,a0fld));
            assert(ischar(a1fld)&&isfield(a1rd,a1fld));
            assert(isstruct(a0RptInfo)&&all(isfield(a0RptInfo,{'dm' 'sdAll' 'sdCtl'})));
            assert(isstruct(a1RptInfo)&&all(isfield(a1RptInfo,{'dm' 'sdAll' 'sdCtl'})));
            
            opts.doPlot = true;
            opts = pvparse({},opts,'error',varargin);
                                    
%             fprintf(1,'a0SdAll a1SdAll: %.5g %.5g\n',a0RptInfo.sdAll,a1RptInfo.sdAll);
%             fprintf(1,'a0SdCtl a1SdCtl: %.5g %.5g\n',a0RptInfo.sdCtl,a1RptInfo.sdCtl);
%             
            [a0RptSD a0RptLns] = OlyDat.Analysis.Repeat.rptInfoRValsForRepeatedLines(a0RptInfo);
            [a1RptSD a1RptLns] = OlyDat.Analysis.Repeat.rptInfoRValsForRepeatedLines(a1RptInfo);
            
            assert(isequal(a0RptLns,a1RptLns));
            linecount.lines = a0RptLns;
            linecount.tfImproved = a1RptSD < a0RptSD;
            linecount.tfWorsened = a1RptSD > a0RptSD;
            linecount.tfUnchanged = a1RptSD==a0RptSD;
            linecount.bars = zeros(numel(linecount.lines),1);
            linecount.bars(linecount.tfImproved) = 1;
            linecount.bars(linecount.tfWorsened) = -1;
            linecount.Nimproved = nnz(linecount.tfImproved);
            linecount.Nworsened = nnz(linecount.tfWorsened);
            linecount.Nunchanged = nnz(linecount.tfUnchanged);
            linecount.str = sprintf('Num lines improved/worsened/unchanged: %d %d %d',...
                linecount.Nimproved,linecount.Nworsened,linecount.Nunchanged);
            disp(linecount.str);
            
            % print a single line for the wiki
            a0RelRptly = mean(a0RptSD)/a0RptInfo.sdAll;
            a1RelRptly = mean(a1RptSD)/a1RptInfo.sdAll;
            fprintf(1,'|| %s/%s || Global GMR SD || CTL SD || Mean GMR Rpt SD || Mdn GMR Rpt SD || Rel. rptbly || Num lines improved/worsened/unchanged || Notes ||\n',a0fld,a1fld);
            fprintf(1,'| No model | %.4g | %.4g | %.4g | %.4g | %.3g | n/a | |\n',...
                a0RptInfo.sdAll,a0RptInfo.sdCtl,mean(a0RptSD),median(a0RptSD),a0RelRptly);
            fprintf(1,'| <model> | %.4g | %.4g | %.4g | %.4g | %.3g | %d/%d/%d | |\n',...
                a1RptInfo.sdAll,a1RptInfo.sdCtl,mean(a1RptSD),median(a1RptSD),a1RelRptly,...
                linecount.Nimproved,linecount.Nworsened,linecount.Nunchanged);
            fprintf(1,'percent rptably change: %.3g, %.3g, %.3g, %.3g\n',...
                (mean(a1RptSD)-mean(a0RptSD))/mean(a0RptSD),...
                (median(a1RptSD)-median(a0RptSD))/median(a0RptSD),...
                (a1RelRptly-a0RelRptly)/a0RelRptly,...
                linecount.Nimproved/((linecount.Nimproved+linecount.Nworsened)/2)-1);
            
            if opts.doPlot
                
                % hist of lineSDs, before+after
                figure('units','normalized','outerposition',[0 0 1 1])
                ax = subplot(2,1,1);
                zlclRValHist(ax,a0RptSD,a0fld);
                axs = axis;
                ax = subplot(2,1,2);
                zlclRValHist(ax,a1RptSD,a1fld);
                axis(axs);
%                 
%                 % superplot2 and temporalplot
%                 figmax;
%                 ax1 = subplot(1,2,1);
%                 tfRepeated = zlclGetTFRepeated(a0rd,a0RptInfo);
%                 OlyDat.Analysis.Repeat.visualizeRepeatability(ax1,a0rd(tfRepeated),a0fld,[],'includeLineAvSig',2);
%                 ax2 = subplot(1,2,2);
%                 tfRepeated = zlclGetTFRepeated(a1rd,a1RptInfo);
%                 OlyDat.Analysis.Repeat.visualizeRepeatability(ax2,a1rd(tfRepeated),a1fld,[],'includeLineAvSig',2);
%                 zlclMakeSuperplotAxesComparable(ax1,ax2);
%                 
                % superplot2, lines in same order
                figmax;
                ax1 = subplot(1,5,1:2);
                tfRepeated = zlclGetTFRepeated(a0rd,a0RptInfo);
                a0lines = OlyDat.Analysis.Repeat.visualizeRepeatability(ax1,a0rd(tfRepeated),a0fld,[],'includeLineAvSig',2);
                ax2 = subplot(1,5,3:4);
                tfRepeated = zlclGetTFRepeated(a1rd,a1RptInfo);
                OlyDat.Analysis.Repeat.visualizeRepeatability(ax2,a1rd(tfRepeated),a1fld,[],'lineOrdering',a0lines,'includeLineAvSig',1,'lineAvSigColor',[0.8 0.8 0.8]);
                zlclMakeSuperplotAxesComparable(ax1,ax2);
                
                % plot a horizontal bar chart for each line, showing whether it
                % improved or worsened.
                ax3 = subplot(1,5,5);
                linecountTmp = linecount;
                linecountTmp.lines{end+1,1} = 'pBDPGAL4U'; % add control line back in
                linecountTmp.bars(end+1,1) = nan; % control line almost always improves, but doesn't matter here
                assert(isequal(sort(a0lines),sort(linecountTmp.lines)));
                [~,loc] = ismember(a0lines,linecountTmp.lines);
                assert(isequal(a0lines,linecountTmp.lines(loc)));
                barh(ax3,1:numel(a0lines),linecountTmp.bars(loc));         
%                 
%                 % superplot2, lines in order of SD
%                 figmax;
%                 ax1 = subplot(1,2,1);
%                 tfRepeated = zlclGetTFRepeated(a0rd,a0RptInfo);
%                 OlyDat.Analysis.Repeat.visualizeRepeatability(ax1,a0rd(tfRepeated),a0fld,[],'includeLineAvSig',0,'lineOrdering','rptSD','shiftToZeroMean',false);
%                 ax2 = subplot(1,2,2);
%                 tfRepeated = zlclGetTFRepeated(a1rd,a1RptInfo);
%                 OlyDat.Analysis.Repeat.visualizeRepeatability(ax2,a1rd(tfRepeated),a1fld,[],'includeLineAvSig',0,'lineOrdering','rptSD','shiftToZeroMean',false);
%                 zlclMakeSuperplotAxesComparable(ax1,ax2);
% 
%                 % superplot2, lines in order of SD, second chart in same
%                 % order as first
%                 figmax;
%                 ax1 = subplot(1,5,1:2);
%                 tfRepeated = zlclGetTFRepeated(a0rd,a0RptInfo);
%                 a0lines = OlyDat.Analysis.Repeat.visualizeRepeatability(ax1,a0rd(tfRepeated),a0fld,[],'includeLineAvSig',2,'lineOrdering','rptSD','shiftToZeroMean',true);
%                 ax2 = subplot(1,5,3:4);
%                 tfRepeated = zlclGetTFRepeated(a1rd,a1RptInfo);
%                 OlyDat.Analysis.Repeat.visualizeRepeatability(ax2,a1rd(tfRepeated),a1fld,[],'includeLineAvSig',2,'lineOrdering',a0lines,'lineAvSigColor',[0.8 0.8 0.8],'shiftToZeroMean',true);
%                 zlclMakeSuperplotAxesComparable(ax1,ax2);
%                 % plot a horizontal bar chart for each line, showing whether it
%                 % improved or worsened.
%                 ax3 = subplot(1,5,5);
%                 assert(isequal(sort(a0lines),sort(linecount.lines)));
%                 [~,loc] = ismember(a0lines,linecount.lines);
%                 assert(isequal(a0lines,linecount.lines(loc)));
%                 barh(ax3,1:numel(a0lines),linecount.bars(loc));      
            end
            
        end
        
    end
    
    methods (Static) % Core API
        
        function lines = visualizeRepeatability(ax,data,fld,expCoordinator,varargin)
            if nargin < 4
                expCoordinator = [];
            end
            
            fprintf(1,'Relying on BoxData plot.\n');
            line2Idx = OlyDat.Analysis.Repeat.line2DataIdx(data);
            ginfo.g = {data.gmrVsCtrl}';
            ginfo.gName = 'gVc';
            lines = BoxData.Variation.experimentLineSuperPlot(ax,data,fld,line2Idx,ginfo,expCoordinator,varargin{:});
        end
        
        % batch processing of repeatability for cut+paste to excel
        function [m rowLbl colLbl] = repeatabilityBatch(rd,flds)
            colLbl = {'globalSD' 'avLineSdRptedLns' 'ctlLineSd'};
                  
            Nfld = numel(flds);
            m = cell(0,3);
            for i = 1:Nfld
                f = flds{i};
                rowLbl{i,1} = f; %#ok<AGROW>
                
                [dmat sdAll sdCtl] = OlyDat.Analysis.Repeat.repeatabilityAllLines(rd,f);
                nrpt = cell2mat(dmat(:,2));
                repeatSd = cell2mat(dmat(:,4));
                avLineSdRptedLns = mean(repeatSd(nrpt>1));
                m(i,:) = {sdAll avLineSdRptedLns sdCtl};
            end
            
        end
        
        % Calculate repeatability for all lines.
        % dmat: Nline x 5 cell matrix. 
        %   Col 1: line name. 
        %   Col 2: nrpt.
        %   Col 3: r (dimensionless repeatability factor). 
        %   Col 4: repeat-sd.
        %   Col 5: repeat-mu.
        % 
        % sdAll: sd of distribution of repeat-means across all lines. THIS
        % INCLUDES LINES NOT REPEATED (Nrpt==1)
        % sdCtl: sd of distribution of control data.
        function [dmat sdAll sdCtl] = repeatabilityAllLines(data,fld)
            line2Idx = OlyDat.Analysis.Repeat.line2DataIdx(data);
            [lines nrpt mu sd] = OlyDat.Analysis.Repeat.lineDistribution(data,fld,line2Idx);
            
            % compute sdAll, sdCtl, avLineSd
            sdAll = std(mu);
            tfCtl = strcmp(lines,OlyDat.Analysis.Repeat.CONTROL_LINE);
            assert(nnz(tfCtl)==1);
            sdCtl = sd(tfCtl);
            
            % generate dmat            
            Nline = numel(lines);            
            dmat = cell(Nline,5);
            dmat(:,1) = lines(:);
            for lIdx = 1:Nline
                %fprintf(1,'%d\n',lIdx);
                dmat(lIdx,2:5) = {nrpt(lIdx) sd(lIdx)/sdAll sd(lIdx) mu(lIdx)};
            end
        end
        
        % fld: a scalar stat that is either >0 for positive hit, ==0 for
        % nonhits, or <0 for negative hits.        
        function [dm dmHdr] = repeatabilityDiscrete(rd,fld)
            line2Idx = OlyDat.Analysis.Repeat.line2DataIdx(rd);
            lines = line2Idx.keys;
            if ismember('pBDPGAL4U',lines)
                warning('OlyDat:Analysis:Repeat:ignoringCtlLine',...
                    'Ignoring control line.');
                lines = setdiff(lines,'pBDPGAL4U');
            end
            Nline = numel(lines);
                        
            dmHdr = {'line' 'Nrpt' 'NhitLow' 'NhitHi' 'Nnonhit' 'hitStatStr' 'hitStat'};
            dm = cell(Nline,7);
            for lnIdx = 1:Nline
                ln = lines{lnIdx};
                idx = line2Idx(ln);
                
                hs = [rd(idx).(fld)];
                assert(numel(hs)==numel(idx));
                
                dm{lnIdx,1} = ln;
                dm{lnIdx,2} = numel(idx); % nrpt
                dm{lnIdx,3} = nnz(hs<0); % Nhit low
                dm{lnIdx,4} = nnz(hs>0); % Nhit hi
                dm{lnIdx,5} = nnz(hs==0); % Nnonhit
                dm{lnIdx,6} = OlyDat.Analysis.Repeat.score2SignStr(hs);
                dm{lnIdx,7} = hs;
            end
            
            % Compute stats
            nrpt = cell2mat(dm(:,2));
            nhitlo = cell2mat(dm(:,3));
            nhithi = cell2mat(dm(:,4));
            nnonhit = cell2mat(dm(:,5));
            stats.Ntot = sum(nrpt); % these are counts of experiments, not lines
            stats.Nhitlo = sum(nhitlo); 
            stats.Nhithi = sum(nhithi);
            stats.Nnonhit = sum(nnonhit);
            stats.NhitloFrac = stats.Nhitlo/stats.Ntot;
            stats.NhithiFrac = stats.Nhithi/stats.Ntot;
            stats.NnonhitFrac = stats.Nnonhit/stats.Ntot;
            assert(stats.Ntot==stats.Nhitlo+stats.Nhithi+stats.Nnonhit);
            fprintf(1,'%d lines.\n',Nline);
            fprintf(1,'EXPERIMENT counts:\n');
            fprintf(1,'Nhit-, Nhit+, Nnonhits, Ntot: %d, %d, %d, %d\n',...
                stats.Nhitlo,stats.Nhithi,stats.Nnonhit,stats.Ntot);
            fprintf(1,'Nhit-/Ntot, Nhit+/Ntot, Nnonhit/Ntot: %.3g, %.3g, %.3g\n',...
                stats.NhitloFrac,stats.NhithiFrac,stats.NnonhitFrac);
            
            % Compare experimentally realized hit distribution vs
            % theoretical/random distro with "zero" repeatability
            for nr = 2:max(nrpt)
                tfThisNRpt = nr==nrpt;
                OlyDat.Analysis.Repeat.repeatabilityDiscreteAssess(...
                    nr,...
                    nhitlo(tfThisNRpt),...
                    nhithi(tfThisNRpt),...
                    nnonhit(tfThisNRpt));
            end
        end
        
        function repeatabilityDiscreteAssess(nrpt,nhitlo,nhithi,nnonhit)
            assert(isequal(size(nhitlo),size(nhithi),size(nnonhit)));
            assert(all( sum([nhitlo(:) nhithi(:) nnonhit(:)],2)==nrpt ));

            stats.Nline = numel(nhitlo);
            stats.Ntot = nrpt*stats.Nline; % experiment
            stats.NhitloFrac = sum(nhitlo)/stats.Ntot;
            stats.NhithiFrac = sum(nhithi)/stats.Ntot;
            stats.NnonhitFrac = sum(nnonhit)/stats.Ntot;
            
            fprintf('## Nrpt==%d, %d experiments, %d lines:\n',nrpt,stats.Ntot,stats.Nline);
            fprintf(1,'Experiment fracs: Nhit-/Ntot, Nhit+/Ntot, Nnonhit/Ntot: %.3g, %.3g, %.3g\n',...
                stats.NhitloFrac,stats.NhithiFrac,stats.NnonhitFrac);            
            
            lineCountActual = [];
            lineCountTheory = [];
            for n0 = 0:nrpt
            for np = 0:nrpt-n0
                nm = nrpt-n0-np;
                str = [repmat('o',1,n0) repmat('+',1,np) repmat('-',1,nm)];
                lineCountActual(end+1) = nnz(nhitlo==nm & nhithi==np & nnonhit==n0); %#ok<AGROW>
                lineCountTheory(end+1) = stats.Nline * ...
                                   stats.NnonhitFrac^n0 * ...
                                   stats.NhitloFrac^nm * ...
                                   stats.NhithiFrac^np * ...
                                   factorial(nrpt) / (factorial(n0)*factorial(np)*factorial(nm)); %#ok<AGROW>
                %fprintf(1,' %s: Exp: %d. Random: %.1f.\n',str,lineCountActual(end),lineCountTheory(end));
                fprintf(1,'| %s | %d | %.1f |\n',str,lineCountActual(end),lineCountTheory(end));
            end
            end
            
            fprintf(1,' Sum Exp: %d. Sum Random: %.1f.\n',sum(lineCountActual),sum(lineCountTheory));
        end
        
    end
    
    % helpers
    methods (Static,Hidden)
        
        function str = score2SignStr(z)
            validateattributes(z,{'numeric'},{'vector'});
            str = char(zeros(size(z)));
            for c = 1:numel(z)
                if z(c)<0
                    str(c) = '-';
                elseif z(c)==0
                    str(c) = 'o';
                else
                    str(c) = '+';
                end
            end 
        end
        
        % rptInfo: struct with fields 'dm', 'sdAll', 'sdCtl', outputs of
        % OlyDat.Analysis.Repeat.repeatabilityAllLines.
        %
        % rsRpted: repeat-sd for repeated lines
        % lns: linenames for repeated lines
        function [rsRpted lns] = rptInfoRValsForRepeatedLines(rptInfo)
            dm = rptInfo.dm;
            assert(strcmp(dm{end,1},'pBDPGAL4U'));
            dm = dm(1:end-1,:);
            
            rs = cell2mat(dm(:,4)); % repeat-sd
            nrpts = cell2mat(dm(:,2));
            assert(isequal(nrpts==1,isnan(rs)));
            lns = dm(:,1);
            
            tfRpted = nrpts>1;
            rsRpted = rs(tfRpted);
            lns = lns(tfRpted);
        end
                
        % For each line, calculate the mean, sd, etc over its repeats.
        %
        % line2Idx: optional, map from lines to data idxs.
        %
        % lines: col cellstr. unique lines in data.
        % nrpt: col vec. number of repeats for each line.
        % mu: col vec. repeat-average for each line.
        % sd: col vec. repeat-sd for each line.
        function [lines nrpt mu sd] = lineDistribution(data,fld,line2Idx)
            assert(ischar(fld));
            assert(isstruct(data) && isfield(data,fld));
            if nargin < 3            
                line2Idx = OlyDat.Analysis.Repeat.line2DataIdx(data);
            else
                assert(isa(line2Idx,'containers.Map'));
            end
                
            lines = line2Idx.keys;
            seen = zeros(size(data)); % for assertion
            for lIdx = numel(lines):-1:1
                ln = lines{lIdx};
                idx = line2Idx(ln);
                assert(~isempty(idx));
                stat = cat(1,data(idx).(fld));
                assert(isequal(size(stat),[numel(idx) 1]));
                
                nrpt(lIdx) = numel(stat);
                mu(lIdx) = nanmean(stat);
                sd(lIdx) = nanstdBiasCorrected(stat);
                
                seen(idx) = seen(idx)+1;
            end
            
            assert(all(seen==1));
        end
        
        function line2Idx = line2DataIdx(data)
            line2Idx = OlyDat.Analysis.fieldVal2DataIdx(data,OlyDat.Analysis.Repeat.LINE_FLD);
        end
        
    end
    
end

function zlclRValHist(ax,rs,fld)
hist(ax,rs,30);
titlestr = sprintf('%s. %d lines. mean std %.5g %.5g',fld,numel(rs),mean(rs),std(rs));
title(ax,titlestr,'interpreter','none','FontSize',18);
end

function tfKeep = zlclGetTFRepeated(rd,ri)
ulines = ri.dm(:,1);
ulinesnrpts = cell2mat(ri.dm(:,2));
ulinesrpted = ulines(ulinesnrpts>1);

allLines = {rd.line_name}';
tfKeep = ismember(allLines,ulinesrpted);
end

function zlclMakeSuperplotAxesComparable(ax1,ax2)
axs1 = axis(ax1);
axs2 = axis(ax2);
delx1 = axs1(2)-axs1(1);
delx2 = axs2(2)-axs2(1);
delx = max(delx1,delx2);

axs1(1:2) = [mean(axs1(1:2))-delx/2 mean(axs1(1:2))+delx/2];
axis(ax1,axs1);
axs2(1:2) = [mean(axs2(1:2))-delx/2 mean(axs2(1:2))+delx/2];
axis(ax2,axs2);
end