classdef AverageEthogram
  
  methods (Static)  
    
    function [hFig,s] = make(data,plottype,pred)
      if ~exist('pred','var')
        tf = true(size(data));
      elseif ischar(pred)
        assert(isfield(data,pred),'No field ''%s'' in data.',pred);
        tf = [data.(pred)]';
        fprintf(1,'Including data based on field ''%s'': including %d/%d experiments.\n',pred,nnz(tf),numel(tf));
      else
        assert(false,'Unknown type for predicate.');
      end
      d = data(tf);
      
      % get groups
      stype = {d.auto_Grab_successtype}';
      g = struct();
      g.suc1 = strcmp(stype,'successful_1grab');
      g.suc2 = strcmp(stype,'successful_2plusgrab');
      g.suc = [d.auto_Grab_success]';
      g.fail = ~g.suc;      
      grps = fieldnames(g);
      Ngrp = numel(grps);
      
      % setup axes
      hFig = figure;
      TOPMARGIN = .1;
      BOTMARGIN = .1;
      LEFTPOS = .13;
      WIDTH = .8;
      YFUDGE = .005;
      height = (1-TOPMARGIN-BOTMARGIN)/Ngrp-YFUDGE;
      axs = nan(Ngrp,1);
      for i = 1:Ngrp
       axs(Ngrp-i+1,1) = subplot('position',[LEFTPOS BOTMARGIN+(i-1)*(height+YFUDGE) WIDTH height]);
      end
      
      XMAX = 5000;
      COLORS = [
        0         0    1.0000; ...
        0    0.5000         0; ...
        1.0000         0         0; ...
        0    0.7500    0.7500; ...
        0.7500         0    0.7500; ...
        0.7500    0.6000         0; ...
        0.2500    0.2500    0.2500];
      s = struct();

      for iGrp = 1:Ngrp
        group = grps{iGrp};
        tf = g.(group);
        NInc = nnz(tf);
        dInc = d(tf);
        
        % collect data        
        maxf = -inf;
        bigI = zeros(0,XMAX);
        for beh = ExpPP.BASICBEHAVIORS(1:end),beh=beh{1}; %#ok<FXSET>
          switch plottype
%             case 'GSSS'
%               fldstart = sprintf('auto_GSSS_%s_0',beh);
%               fldbl = sprintf('auto_GSSS_%s_bl',beh);
%               s.(beh).start = [d(tf).(fldstart)]';
%               s.(beh).bl = [d(tf).(fldbl)]';
%               s.(beh).end = s.(beh).start+s.(beh).bl;        
%               assert(isequal(numel(s.(beh).start),numel(s.(beh).bl),numel(s.(beh).end),NInc));
%               assert(all(s.(beh).end<XMAX));
% 
%               s.(beh).I = zeros(NInc,XMAX);
%               for i = 1:NInc
%                 s.(beh).I(i,s.(beh).start(i):s.(beh).end(i)) = 1;          
%               end
            case 'first'
              titlestr = 'Average Ethogram with First Bouts';
              fldstart = sprintf('auto_%s_t0s',beh);
              fldend = sprintf('auto_%s_t1s',beh);
              fldb0 = sprintf('auto_%s_0',beh);
              fldb0len = sprintf('auto_%s_0_bl',beh);

              s.(group).(beh).I = zeros(NInc,XMAX); % I(dataIdx,xIdx)=1 if dataIdx'th experiment has a bout of behavior at t=xIdx
%               s.(group).(beh).Ibl = nan(NInc,XMAX); % If dataIdx'th experiment has a bout of behavior at t=xIdx, I(dataIdx,xIdx)=bl; otherwise i(dataIdx,xIdx)=nan
%               s.(group).(beh).bls = cell(NInc,1);
              s.(group).(beh).Nbout = nan(NInc,1);
%               s.(group).(beh).Istart = zeros(NInc,XMAX);              
              s.(group).(beh).b0 = nan(NInc,1);
              s.(group).(beh).b0len = nan(NInc,1);
              for iDat = 1:NInc
                s.(group).(beh).b0(iDat) = dInc(iDat).(fldb0);
                s.(group).(beh).b0len(iDat) = dInc(iDat).(fldb0len);
                
                strs = dInc(iDat).(fldstart);
                ends = dInc(iDat).(fldend);
                assert(isequal(size(strs),size(ends)));
                Nbout = numel(strs);
                s.(group).(beh).Nbout(iDat) = Nbout;
                for iBout = 1:Nbout
                  idx = strs(iBout):ends(iBout);
%                   bl = ends(iBout)-strs(iBout);
                  assert(all(s.(group).(beh).I(iDat,idx)==0));
                  s.(group).(beh).I(iDat,idx) = 1;         
%                   s.(group).(beh).Ibl(iDat,idx) = bl;
%                   s.(group).(beh).bls{iDat}(end+1) = bl;
% 
%                   idxStart = strs(iBout):strs(iBout)+4;
%                   s.(group).(beh).Istart(iDat,idxStart) = 1;
                end
              end

              s.(group).(beh).f = sum(s.(group).(beh).I,1)/NInc;
              maxf = max(maxf,max(s.(group).(beh).f));
              s.(group).(beh).w = ones(size(s.(group).(beh).f));              
            case 'grabseq'
              titlestr = 'Average Ethogram with Grab Sequences';
              fldstart = sprintf('auto_%s_t0s',beh);
              fldend = sprintf('auto_%s_t1s',beh);
              fldgs00b0 = sprintf('auto_GS00_%s_0',beh);
              fldgs00b0len = sprintf('auto_GS00_%s_bl',beh);
              fldgsssb0 = sprintf('auto_GSSS_%s_0',beh);
              fldgsssb0len = sprintf('auto_GSSS_%s_bl',beh);              

              s.(group).(beh).I = zeros(NInc,XMAX); % I(dataIdx,xIdx)=1 if dataIdx'th experiment has a bout of behavior at t=xIdx
%               s.(group).(beh).Ibl = nan(NInc,XMAX); % If dataIdx'th experiment has a bout of behavior at t=xIdx, I(dataIdx,xIdx)=bl; otherwise i(dataIdx,xIdx)=nan
%               s.(group).(beh).bls = cell(NInc,1);
              s.(group).(beh).Nbout = nan(NInc,1);
%               s.(group).(beh).Istart = zeros(NInc,XMAX);
              
              s.(group).(beh).gs00b0 = nan(NInc,1);
              s.(group).(beh).gs00b0len = nan(NInc,1);
              s.(group).(beh).gsssb0 = nan(NInc,1);
              s.(group).(beh).gsssb0len = nan(NInc,1);
              for iDat = 1:NInc
                if isfield(dInc,fldgs00b0)
                  s.(group).(beh).gs00b0(iDat) = dInc(iDat).(fldgs00b0);
                  s.(group).(beh).gs00b0len(iDat) = dInc(iDat).(fldgs00b0len);
                end
                if isfield(dInc,fldgsssb0)
                  s.(group).(beh).gsssb0(iDat) = dInc(iDat).(fldgsssb0);
                  s.(group).(beh).gsssb0len(iDat) = dInc(iDat).(fldgsssb0len);
                end
                strs = dInc(iDat).(fldstart);
                ends = dInc(iDat).(fldend);
                assert(isequal(size(strs),size(ends)));
                Nbout = numel(strs);
                s.(group).(beh).Nbout(iDat) = Nbout;
                for iBout = 1:Nbout
                  idx = strs(iBout):ends(iBout);
%                   bl = ends(iBout)-strs(iBout);
                  assert(all(s.(group).(beh).I(iDat,idx)==0));
                  s.(group).(beh).I(iDat,idx) = 1; 
%                   s.(group).(beh).Ibl(iDat,idx) = bl;
%                   s.(group).(beh).bls{iDat}(end+1) = bl;
% 
%                   idxStart = strs(iBout):strs(iBout)+4;
%                   s.(group).(beh).Istart(iDat,idxStart) = 1;
                end
              end

              s.(group).(beh).f = sum(s.(group).(beh).I,1)/NInc;
              maxf = max(maxf,max(s.(group).(beh).f));
              s.(group).(beh).w = ones(size(s.(group).(beh).f));              
          end

          bigI = [bigI;s.(group).(beh).I]; %#ok<AGROW>
        end
        s.(group).NInc = NInc;
        s.(group).maxf = maxf;
      
%         hFig = figure;
%         ax = axes;
%         spy(bigI);
%         titlestr = sprintf('%d/%d exps',NInc,numel(d));
%         title(titlestr,'fontweight','bold');
%         grid on;
        
%         hFig2 = figure;
%         ax = axes;
        ax = axs(iGrp);
        hold(ax,'on');
        MINNUM_STARTANDLEN = 5; % require this many experiments for bout start/len
        for iBeh = 1:6
          beh = ExpPP.BASICBEHAVIORS{iBeh};
          AverageEthogram.plotBout(ax,-iBeh,1:XMAX,s.(group).(beh).f,s.(group).(beh).w,COLORS(iBeh,:));
          switch plottype
            case 'first'              
              boutStarts = s.(group).(beh).b0;
              boutLens = s.(group).(beh).b0len;
              assert(nnz(~isnan(boutStarts))==nnz(~isnan(boutLens)));
              if nnz(~isnan(boutStarts))>=MINNUM_STARTANDLEN
                AverageEthogram.plotBoutStartAndLen(ax,-iBeh,boutStarts,boutLens,[0 0 0]);
              end
            case 'grabseq'
              
              boutStarts = s.(group).(beh).gs00b0;
              boutLens = s.(group).(beh).gs00b0len;
              assert(nnz(~isnan(boutStarts))==nnz(~isnan(boutLens)));
              if nnz(~isnan(boutStarts))>=MINNUM_STARTANDLEN
                b075 = AverageEthogram.plotBoutStartAndLen(ax,-iBeh+.165,boutStarts,boutLens,[0 0 0]);
                tfGrabSeqLbl = iGrp==1 && iBeh==1;
                if tfGrabSeqLbl
                  text(b075+40,-iBeh+.25,'GS00','fontsize',7,'fontangle','italic');
                end
              end               
          
              boutStarts = s.(group).(beh).gsssb0;
              boutLens = s.(group).(beh).gsssb0len;
              assert(nnz(~isnan(boutStarts))==nnz(~isnan(boutLens)));
              if nnz(~isnan(boutStarts))>=MINNUM_STARTANDLEN
                b075 = AverageEthogram.plotBoutStartAndLen(ax,-iBeh-.3,boutStarts,boutLens,[0 0 0]);
                tfGrabSeqLbl = iGrp==1 && iBeh==1;
                if tfGrabSeqLbl
                  text(b075+40,-iBeh-.25,'GSSS','fontsize',7,'fontangle','italic');
                end
              end
          end
        end
        lblstr = sprintf('%s (n=%d)',group,NInc);
        ylabel(ax,lblstr,'fontweight','bold');
        grid(ax,'on');
        set(ax,'YTick',[],'YTickLabel',[]);
        if iGrp==Ngrp
          %set(ax,'XTickLabel',500);
        else
          set(ax,'XTickLabel',[]);
        end
        if iGrp==1
          title(ax,titlestr,'fontweight','bold');
        end
        
      end
      linkaxes(axs);
      ylim(axs(1),[-7 0]);
    end
      
    function b075 = plotBoutStartAndLen(ax,y,b0,bl,c)
      % plot a bout length bar of length bl starting at b0, indicating
      % dispersion in bls.
      %
      % bl: all boutlenghts
      % b0: all boutstarts
      
      b0md = nanmedian(b0);
      b025 = prctile(b0,25);
      b075 = prctile(b0,75);      
      blmd = nanmedian(bl);
      
      LINEWIDTH = 4;
      axes(ax);
      AverageEthogram.herrbar(b0md,y,b0md-b025,b075-b0md,c);
      plot([b0md b0md+blmd],[y y],'-','Color',c,'LineWidth',LINEWIDTH);
    end
    
    function herrbar(x,y,l,u,c)
      YHEIGHT = .28;
      h = errorbar_x(x,y,l,u,'.');
      h(1).Color = c;
      h(2).Color = c;
      h(1).YData([4 7]) = y-YHEIGHT/30;%-YHEIGHT;
      h(1).YData([5 8]) = y+YHEIGHT;
      h(1).Tag = 'boutstartEB';
      h(2).Tag = 'boutstartEBMarker';
    end
      
    function plotBout(ax,y,x,f,w,c,varargin)
      assert(isscalar(y));
      assert(isequal(size(x),size(f)));
%       [fmin,fmax,varwidth] = myparse(varargin,...
%         'fmin',0,...
%         'fmax',1,...
%         'varwidth',true);
      
      a = f;
      
      % verts
      Nx = numel(x);
      Nv = 2*Nx;
      verts = nan(Nv,2);
      verts(1:2:end,1) = x;
      verts(2:2:end,1) = x;
      verts(1:2:end,2) = y-w/2;
      verts(2:2:end,2) = y+w/2;
      
      % faces
      Nf = Nx-1;
      faces = nan(Nf,4); % trapezoids
      for i = 1:Nf
        vIdx = 2*i-1;
        faces(i,:) = [vIdx vIdx+1 vIdx+3 vIdx+2];
      end
      
      fvalpha = nan(Nv,1);
      fvalpha(1:2:end) = a;
      fvalpha(2:2:end) = a;
      
      patch('Parent',ax,'Vertices',verts,'Faces',faces,...
        'FaceVertexAlphaData',fvalpha,'FaceAlpha','interp','EdgeColor','none','FaceColor',c);      
    end
    
  end  
  
end
