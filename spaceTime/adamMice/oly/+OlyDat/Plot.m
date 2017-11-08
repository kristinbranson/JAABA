classdef Plot < handle
%Plot OlyDat Browser Plot
%   A Plot object represents a plot/analysis in the OlyDat Browser. A plot
%   has various responsibilities (all except for the actual plotting are
%   "as appropriate"):
%
%   * Custom data preprocessing 
%   * The actual plotting
%   * Interactivity: sending and responding to experiment-selection signals
%   * Math/statistical analysis
%
%   See also OlyDat.Plot/preprocessData, OlyDat.Plot/doPlot,
%   OlyDat.Browser, OlyDat.BrowserStat, OlyDat.ExperimentCoordinator.
    
    
    properties (Abstract)
        Name;
        AllowsMultiStats; % (logical) If true, this plot handles multiple statistics
        UsesAuxVar; % Either true, false, or 'custom'. 
                    % If true: this plot uses the auxiliary pulldown, as configured at Browser startup.
                    % If false: this plot does not use the aux pulldown. (It will be inactivated.)
                    % If 'custom': This plot configures the auxiliary pulldown itself (using groupVarChanged()). 
		    %	           pltCfg.auxVar will be empty when doPlot() is called. This class can use
		    % 		   auxVarChanged() to cache the current auxvar setting.
		    %		   TODO Requiring this class to cache the the current auxvar is unnecessary.
        
        UsesGroupVar2; % Either true or false
        GroupVarLabel; % Text label for Group Var 1
        GroupVar2Label;
        
        AutoReplotAfterStatChange; % (logical) If true, stat selection causes replot
        AutoReplotAfterGroupChange;
    end
    
    methods
    
        function data = preprocessData(obj,data,plotCfg) %#ok<INUSD,MANU>
            %data = preprocessData(obj,data,plotCfg)
            %   Plots should override this method if they need custom data
            %   manipulation/flattening for their plot. Typical uses
            %   include flattening statistics in the .Scores field to the
            %   top level, averaging scores over tubes, normalizing stats,
            %   etc.
            %
            %   plotCfg is a struct with fields .stat, .grp, .grpVals,
            %   .aux.

            % default implementation: no-op
        end        
       
        function groupChanged(obj,grpStat,allData,hCtrls) %#ok<INUSD,MANU>
            % Called when user selects a new grouping variable. (When the
            % plot is first selected, the current grouping variable is
            % considered to be selected.)
            %
            % grpStat: OlyDat.BrowserStat for new grouping var.
            % allData: struct array, all data loaded in browser.
            % hCtrls: scalar struct:
            %         * hCtrls.pumAuxVar is the auxiliary var pulldown.
            
            % This method is used when a plot needs to eg customize what
            % auxiliary vars are available, depending on the grouping var.
            % Default implementation: no-op
        end
        
        function auxVarChanged(obj,auxVar,allData,hCtrls) %#ok<MANU,INUSD>
            % Called when a user selects a new auxiliary variable. (When
            % the plot is first selected, the current aux variable is
            % considered to be selected.)
            %
            % auxVar: OlyDat.BrowserStat for new aux var. if
            % UsesAuxVar=='custom', auxVar will be [].
            % allData: struct array, all data loaded in browser.
            % hCtrls: scalar struct:
            %           * hCtrls.pumAuxVar is the auxiliary var pulldown.
            
            % Default impl: no-op
        end            
        
    end
    
    methods (Abstract)
        %detailStr = doPlot(obj,ax,data,stat,grp,plotCfg,expCoordinator)
        %   * ax: axes object
        %   * data: preprocessed data struct for current dataset
        %   * stat: Nrows x Nstats numeric matrix. As a convenience the
        %   Browser puts the statistic(s) of interest into a numeric array
        %   for plotting. Nrows is the number of elements in the
        %   preprocessed data struct.
        %   * grp: Nrows x 1 column vector, either numeric or cellstr. As a
        %   convenience the Browser puts the value of the grouping variable
        %   (if any) in a column vector for plotting.
        %   * plotCfg: The plot configuration specified by the user within
        %   the Browser. This is a struct with fields .stat, .grp,
        %   .grpVals, and .aux.
        %   * expCoordinator: ExperimentCoordinator object used by Browser.
        %   Plots that want to participate in plot interactivity should
        %   register appropriately with the coordinator.
        %   * detailStr: Plots should return a string/cellstr that will be
        %   displayed in the Browser. Typical uses include a description of
        %   the plot, or math details of any statistics performed, etc.
        %
        % OlyDat calls this method when drawing a plot in an axis. This can
        % occur either in the main Browser window, or in standalone figures
        % (eg in the case of 'Compare plots'). If doPlot() performs any
        % actions that require cleanup, that cleanup should be put in ax's
        % DeleteFcn. The Browser will ensure this method is called when a
        % plot is going away. Unregistering from the expCoordinator is an
        % example of such a cleanup action.
        %
        % TODO: unregistering from the expCoordinator is a common action,
        % we should make this easier for users.
        str = doPlot(obj,ax,data,bstat,grp,pltCfg,expCoordinator)
    end

    %% Utilities
    
    methods (Static)
      
      function messageAxisCenter(ax,txt)
        % Put a message in the center of axis ax.
        
        axes(ax);
        xl = xlim(ax);
        yl = ylim(ax);
        hTxt = text(mean(xl),mean(yl),txt);
        set(hTxt,'horizontalalignment','center','interpreter','none');
      end
      
    end
    
end