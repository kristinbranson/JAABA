classdef MultiStatCompareBase < OlyDat.Plot
  
  properties
    % map from "original" ax to 2-vectors for subplot-axes. 
    % MultiStatCompare plots need two sub-axes rather than the single large axis
    ax2AxMap;     
  end
  
  methods
    
      function obj = MultiStatCompareBase
        obj.ax2AxMap = containers.Map('KeyType','double','ValueType','any');
      end
    
      function newax = setupAxes(obj,ax)
        % Call this to hide the original big axis and return two subaxes for plotting.
        
        % Browser comes with one axis only; hide it and replace with 2
        axpos = get(ax,'Position');
        axfig = get(ax,'Parent');
        set(ax,'Visible','off','deletefcn',@(zSrc,zEvt)obj.cleanupPlot(zSrc));
        
        XLABEL_SPACEFRAC = 0.25;
        BETWEEN_SPACEFRAC = 0.02;
        newbottom = axpos(2)+XLABEL_SPACEFRAC*axpos(4);
        newheight = axpos(4)*(1-XLABEL_SPACEFRAC-BETWEEN_SPACEFRAC)/2;
        spacing = axpos(4)*BETWEEN_SPACEFRAC;
        newpos1 = axpos;
        newpos1(2) = newbottom;
        newpos1(4) = newheight;
        newpos2 = axpos;
        newpos2(2) = newbottom+newheight+spacing;
        newpos2(4) = newheight;
        
        if obj.ax2AxMap.isKey(double(ax))
          newax = obj.ax2AxMap(double(ax));
          delete(newax(ishandle(newax)));
        end
        newax = nan(2,1);
        newax(2) = axes('Parent',axfig,'box','on','Position',newpos1,'XTick',[],'XTickLabel',[],'YTick',[],'YTickLabel',[]);
        newax(1) = axes('Parent',axfig,'box','on','Position',newpos2,'XTick',[],'XTickLabel',[],'YTick',[],'YTickLabel',[]);
        obj.ax2AxMap(double(ax)) = newax;
      end
      
      function cleanupPlot(obj,ax)
        % To support two subplots rather than one, MultiStatCompare needs
        % to clean up after every plot. Browser calls the axis deleteFcn
        % after every plot.
        
        if obj.ax2AxMap.isKey(double(ax))
          newax = obj.ax2AxMap(double(ax));
          delete(newax(ishandle(newax)));
          obj.ax2AxMap.remove(double(ax));
        end
        
        set(ax,'Visible','on','DeleteFcn',[]);
      end
      
  end
end
