classdef Browser < handle
%Browser Interactively view Olympiad data
%   Given a data set, the Browser lets you choose variables of interest,
%   plot those variables in various ways under various groupings, perform
%   various analyses, etc. You can select experiments of interest and dig
%   into the details for that experiment.
%
%   The Browser is extensible with respect to new plots/analyses. To create
%   a new plot/analysis, subclass OlyDat.Plot.
%
%   See also OlyDat.Browser/Browser, OlyDat.Plot, OlyDat.BrowserStat,
%   OlyDat.ExperimentDetailHandler, OlyDat.DataSelector,
%   ExperimentHighlighter.
 
    properties (Constant)
        %FO_FILE_SYSTEM_ROOT_NFS = '/groups/sciserv/flyolympiad/';
        CURATION_FIELD_NOT_IN_DATA_WID = 'OlyDat:Browser:curationFieldNotInData';
        ACCEPTED_PV_CTOR_ARGS = {'ExperimentFileSystemFolder';'PassAllSelectedExpsToDetailHandler'};
    end

    properties (Constant)
        %FO_FILE_SYSTEM_ROOT_DEFAULT = zlclInitFO_FILE_SYSTEM_ROOT_DEFAULT();
    end
    
	properties (Hidden,SetAccess=private)
        fUsername; % username supplied at login
        
		fBrowserH; % handle to BrowserGUI
		fBrowserGData; % handle to Browser guidata
		fAxes; % handle to axes in BrowserGUI
        
        fAssayName; % char. name of assay, eg 'box', 'bowl', etc
        fData; % Nx1 col struct array. dataset. 
        fDataTF; % Nx1 logical array. Indicators for whether an experiment is currently in consideration (true) or excluded (false)
        fEid2Idx; % map from experiment_id -> idx into fData.
        fExpDateFormat; % string; date format for exp_datetime field        
		
        fPlots; % col array of OlyDat.Plot objects
        fSelectedPlot; % idx into fPlots for currently selected plot
        fExpCoordinator; % scalar ExperimentCoordinator object
        fComparePlots; % col array of handles to currently open "compare plots"
        
        fBStats; % col vec of BrowserStats
        fGroupStats; % col vec of BrowserStats
        fAuxStats; % col vec of BrowserStats
        fExpDetailStats; % col vec of BrowserStats

        fSelectedExpIdxs; % idx into fData for the currently selected experiment.
        fSelectedExpIdxsIdx; % idx into fSelectedExpIdxs for currently viewed experiment in expDetail pane.
        
        fCustomGroupInfo = struct('Name',cell(0,1),'MemberIdx',cell(0,1)); % struct array with fields 'Name' and 'MemberIdx'.
        
        fCurationFieldObjs; % col vec of OlyDat.CurationFields
        fCurationFieldObjsTblOrder; % col cellstr, names of fCurationFieldObjs to be used in Curation Manager Table
        fCurationFieldObjsTSVOrder; % col cellstr, names of fCurationFieldObjs to be used in Curation TSV file
        fCurationInfo; % map from eid to curation info struct
        fCurationMarkedThisSession; % map from eid to logical, curated this session or not
        fCurationUnsavedCurationInfo; % scalar logical; if true, there are edits to fCurationInfo that have not been writen to a curation file.
        fCurationFile; % char. Currently loaded curation filename.
        fCurationWriter; % scalar handle to CurationFileWriter object
        fCurationExpName2ExpIDMap; % map from experiment_name to experiment_id; used for reading/writing curation files
        fCurationH; % handle to BrowserCurationGUI
        fCurationGData; % handle to BrowserCuration guidata
        fCurationInfoEntryGUIName; % guiname for curation info entry dialog
        
        fExpDetailHandler; % scalar OlyDat.ExperimentDetailHandler object
        fPassAllSelectedExpsToDetailHandler = false; % scalar logical. If true, pass all selected exps to detailHandler, not just selected+viewed experiment. 
        
        fExperimentFileSystemFolder = []; % File system location containing experiment folders
    end
    
    properties (Dependent)
        fDataIsLoaded; % scalar logical. Is there currently data loaded
        fNumExps; % scalar int
        fNumCurrentExps; % scalar int
        fNumExpDetailStats; % scalar int
        fCurrentData; % col struct array. Current data set (full data set minus exclusions)
        fAllSelectedData; % struct array. Data for currently selected experiments.
        fAllSelectedEids; % col vec of selected experimentIDs
        fSelectedDisplayedData; % scalar struct. Currently selected, displayed, data.
        fSelectedDisplayedIdx; % idx into fData for currently selected, displayed, exp.
        fBStatUIControl; % scalar handle to uicontrol. Either fBrowserGData.lbBStats or fBrowserGData.pumBStats, depending on which is visible.
        fAllCustomGroupNames; % cell array of current custom group names
        fCurrentPlot; % currently selected OlyDat.Plot object
        fCurrentAuxStat; % currently selected aux BrowserStat object. If fCurrentPlot.UsesAuxVar=='custom', this returns [].
        % TODO The handling of fCurrentAuxStat for 'custom' places an
        % unnecessary burden on OlyDat.Plot. In the case of 'custom' we can
        % still pass the current auxvar string selection to doPlot(),
        % auxVarChanged(), etc, rather than passing [].
    end
    
    properties (Hidden,Dependent)
        fCurationFields; % col cellstr of all curation field names
        fCurationWritableFields; % col cellstr of writable curation fields
        fCurationReadonlyFields; % col cellstr of readonly curation fields
    end
        
    % Dependent property accessors
    methods
        function tf = get.fDataIsLoaded(obj)
            tf = ~isempty(obj.fData);
        end
        function v = get.fNumExps(obj)
            v = numel(obj.fData);
        end
        function v = get.fNumCurrentExps(obj)
            v = nnz(obj.fDataTF);            
        end
        function v = get.fNumExpDetailStats(obj)
            v = numel(obj.fExpDetailStats);
        end
        function v = get.fCurrentData(obj)
            v = obj.fData(obj.fDataTF);
        end
        function v = get.fAllSelectedData(obj)
            if isempty(obj.fSelectedExpIdxs)
                v = [];
            else
                v = obj.fData(obj.fSelectedExpIdxs);
            end
        end
        function v = get.fAllSelectedEids(obj)
            d = obj.fAllSelectedData;
            if isempty(d)
                v = zeros(0,1);
            else
                v = [d.experiment_id]';
            end
        end
        function v = get.fSelectedDisplayedData(obj)
            if isempty(obj.fSelectedExpIdxs)
                v = [];
            else
                v = obj.fData(obj.fSelectedExpIdxs(obj.fSelectedExpIdxsIdx));
            end
        end
        function v = get.fSelectedDisplayedIdx(obj)
            if isempty(obj.fSelectedExpIdxs)
                v  = [];
            else
                v = obj.fSelectedExpIdxs(obj.fSelectedExpIdxsIdx);
            end
        end
        function v = get.fAllCustomGroupNames(obj)
            v = {obj.fCustomGroupInfo.Name}';
        end
        function h = get.fBStatUIControl(obj)
            if strcmp(get(obj.fBrowserGData.lbBStats,'Visible'),'on')
                h = obj.fBrowserGData.lbBStats;
            else
                assert(strcmp(get(obj.fBrowserGData.pumBStats,'Visible'),'on'));
                h = obj.fBrowserGData.pumBStats;
            end
        end
        function v = get.fCurationFields(obj)
            v = {obj.fCurationFieldObjs.Name}';
        end
        function v = get.fCurationWritableFields(obj)
            tf = [obj.fCurationFieldObjs.TFEditable]';
            flds = obj.fCurationFieldObjs(tf);
            v = {flds.Name}';
        end
        function v = get.fCurationReadonlyFields(obj)
            tf = [obj.fCurationFieldObjs.TFEditable]';
            flds = obj.fCurationFieldObjs(~tf);
            v = {flds.Name}';
        end
        function v = get.fCurrentPlot(obj)
            v = obj.fPlots{obj.fSelectedPlot};
        end
        function v = get.fCurrentAuxStat(obj)
            crntPlot = obj.fCurrentPlot;
            if isequal(crntPlot.UsesAuxVar,'custom')
                v = [];
            else
                idx = get(obj.fBrowserGData.pumAuxVar,'Value');
                v = obj.fAuxStats(idx);
            end
        end
    end
    
    methods (Access=private)
        function delete(obj)
            obj.fExpDetailHandler.close();
            delete(obj.fCurationH);
            delete(obj.fComparePlots); % maybe not strictly necessary
            delete(obj.fBrowserH);
        end
    end
	
    methods (Static)
        
        % data: struct data array 
        % fldname: name of field to extract
        % bstat: numel(data)-by-dataWidth array with nans filled in
        % appropriately when data(idx).fldname is empty.
        function bstat = extractSingleBStat(data,dataFieldName)
            bStatData = {data.(dataFieldName)}';
            tfMissingData = cellfun(@isempty,bStatData);
            if all(tfMissingData)
                % All data fields are empty; assume unit data width
                bstat = nan(numel(data),1);
            else
                nonEmptyData = cat(1,bStatData{~tfMissingData}); % will throw if data have different widths
                assert(size(nonEmptyData,1)==nnz(~tfMissingData));
                dataWidth = size(nonEmptyData,2);
                bstat = nan(numel(data),dataWidth);
                bstat(~tfMissingData,:) = nonEmptyData;
            end            
        end
        
    end
    
	methods
		
        function setUsername(obj,uname)
            assert(ischar(uname));
            obj.fUsername = uname;            
        end
        
		function obj = Browser(varargin)
        %Browser Construct a Browser and bring up the Browser GUI.
        %   obj = Browser(assayName,data,plots,bstats,groupstats,auxstats,
        %   expdetailstats,expdetailhandler,curationfldinfo,curationinfoentryguiname,varargin) 
        %   constructs a browser object.
        %   assayName is eg 'box' or 'bowl' or 'gap', etc. data is a data
        %   set, typically pulled from SAGE using the OlyDat.DataSelector.
        %   Plots is a cell array of OlyDat.Plot objects. bstats,
        %   groupstats, auxstats, expdetailstats are arrays of
        %   OlyDat.BrowserStat objects. Expdetailhandler is an
        %   OlyDat.ExperimentDetailHandler object. Curationfldinfo is a
        %   struct with three fields: .objs, .tblOrder, and .tsvOrder:
        %       * curationfldinfo.objs is an array of OlyDat.CurationField
        %       objects.
        %       * curationfldinfo.tblOrder is a cellstr of the names of
        %       curationfldinfo.objs to be used in the curation manager
        %       table.
        %       * curationfldinfo.tsvOrder is a cellstr of the names of the
        %       curationfldinfo.objs to be used in the curation TSV file.        
        %   Curationinfoentryguiname is the gui to use for the
        %   curationInfoEntry dialog.
        %
        %  varargin accepted PV pairs:
        %  'ExperimentFileSystemFolder': (local) path containing experiment folders
        %  'PassAllSelectedExpsToDetailHandler': scalar logical
        
            error(nargchk(10,inf,nargin,'struct'));
                        
            assayName = varargin{1};
            data = varargin{2};
            plots = varargin{3};
            bstats = varargin{4};
            groupstats = varargin{5};
            auxstats = varargin{6};
            expdetailstats = varargin{7};
            expdetailhandler = varargin{8};
            curationfldinfo = varargin{9};
            curationinfoentryguiname = varargin{10};
            varargin = varargin(11:end);
            
            assert(ischar(assayName),'''assayName'' argument must be a string.');
            assert(isstruct(data),'''data'' argument must be a struct array.');            
            assert(isstruct(curationfldinfo) && ...
                   all(isfield(curationfldinfo,{'objs';'tblOrder';'tsvOrder'})) && ...
                   isa(curationfldinfo.objs,'OlyDat.CurationField') && ...
                   iscellstr(curationfldinfo.tblOrder) && ...
                   iscellstr(curationfldinfo.tsvOrder));
            
%             LoginGUI(@(x)obj.setUsername(x));
%             if isempty(obj.fUsername)
%                 % for now, continue with empty username
%             end
        
            obj.fBrowserH = BrowserGUI(obj);
            obj.fBrowserGData = guidata(obj.fBrowserH);
			obj.fAxes = obj.fBrowserGData.ax;
                        
            obj.fAssayName = assayName;
            obj.fPlots = plots;
            plotNames = cellfun(@(x)x.Name,obj.fPlots,'UniformOutput',false);
			set(obj.fBrowserGData.pumPlotType,'String',plotNames);
            set(obj.fBrowserGData.pumPlotType,'Value',1);
            obj.fSelectedPlot = 1;
            obj.fExpCoordinator = OlyDat.ExperimentCoordinator();
            obj.fExpCoordinator.registerObject(obj,'handleExpSelectionLimited');
            obj.fComparePlots = zeros(0,1);
            
            obj.fBStats = bstats(:);
            obj.fGroupStats = groupstats(:); %[OlyDat.BrowserStat('','<no group>');groupstats(:)];
            obj.fAuxStats = auxstats(:);
            obj.fExpDetailStats = expdetailstats(:);
                        
            % BStat lb/pum
            set(obj.fBrowserGData.lbBStats,'Visible','off'); % cosmetic, minimize time of overlapping uicontrols
			set(obj.fBrowserGData.lbBStats,'String',{obj.fBStats.PrettyName}');
            set(obj.fBrowserGData.pumBStats,'String',{obj.fBStats.PrettyName}');
            set(obj.fBrowserGData.lbBStats,'Value',1);
            set(obj.fBrowserGData.pumBStats,'Value',1);
            set(obj.fBrowserGData.lbBStats,'Max',2); % enable multi-select
            
            set(obj.fBrowserGData.pumGroup,'String',{obj.fGroupStats.PrettyName}');
            set(obj.fBrowserGData.pumGroup,'Value',1);
            set(obj.fBrowserGData.pumGroup2,'String',{obj.fGroupStats.PrettyName}');
            set(obj.fBrowserGData.pumGroup2,'Value',1);
            obj.refreshGroupsList();
            
            % aux vars
            obj.fAuxStats = auxstats;
            set(obj.fBrowserGData.pumAuxVar,'String',{obj.fAuxStats.PrettyName}');
            			
			% expfields            
			obj.fExpDetailStats = expdetailstats;

            expDetailData = {expdetailstats.PrettyName}';
            expDetailData(:,2) = {''};
            set(obj.fBrowserGData.tblExpDetail,'Data',expDetailData);
            
            obj.fExpDetailHandler = expdetailhandler;

            obj.loadData(data);
            
            obj.handleExpSelectionLimited(zeros(0,1)); % start with no selection
            
            % Curation: browser curation state
            obj.fCurationFieldObjs = curationfldinfo.objs;
            obj.fCurationFieldObjsTblOrder = curationfldinfo.tblOrder;
            obj.fCurationFieldObjsTSVOrder = curationfldinfo.tsvOrder;
            obj.initCurationInfoAndMarkedThisSessionAndUnsavedCurationInfo();
            
            % Curation: curation file/writer
            obj.fCurationFile = '';
            obj.initCurationWriter();
            obj.initCurationExpName2ExpIDMap();
            
            % Curation: Browser Curation guis
            obj.initCurationTable();            
            obj.fCurationInfoEntryGUIName = curationinfoentryguiname;
            
            % optional PVs
            assert(mod(numel(varargin),2)==0);
            for c = 1:2:numel(varargin)
                pname = varargin{c};
                val = varargin{c+1};
                assert(ischar(pname) && ismember(pname,obj.ACCEPTED_PV_CTOR_ARGS));
                objpropname = ['f' pname];
                obj.(objpropname) = val;
            end
            
            % FOR ADAM 
            set(obj.fBrowserGData.pbEnterCurationInfo,'Enable','off');
            set(obj.fBrowserGData.pbCreateCustomGroup,'Enable','off');
            set(obj.fBrowserGData.pbLineReport,'Enable','off');
            set(obj.fBrowserGData.pbFileSystem,'Enable','off');
            
            obj.newPlot(1);
        end
        
        
        % Check that data has the necessary fields, etc to work with the
        % Browser. This may throw.        
        function data = checkDataForBrowserRequiredFields(obj,data)
            
            % Note: Browser hardcodes against experiment_id,
            % experiment_name, exp_datetime, and fields in
            % fCurationExperimentIdentifyingFields
            
            % experiment_id-- used in curation, and in experiment-selection
            % signaling. mandatory field
            assert(isfield(data,'experiment_id'),'Data must have field ''experiment_id''.');
            cellEids = {data.experiment_id}';
            assert(all(cellfun(@(x)isnumeric(x)&&isscalar(x),cellEids)),'''experiment_id'' must have scalar numeric value.');
            matEids = cellfun(@double,cellEids);
            assert(numel(unique(matEids))==numel(matEids),'Data has duplicate values of ''experiment_id''.');
            
            % experiment_name. mandatory field
            assert(isfield(data,'experiment_name'),'Data must have field ''experiment_name''.');
            eNames = {data.experiment_name}';
            assert(iscellstr(eNames),'''experiment_name'' must have a string value.');
            assert(numel(unique(eNames))==numel(eNames),'Data has duplicate values of ''experiment_name''.');
            
            % exp_datetime. optional field
            obj.fExpDateFormat = '';
            if isfield(data,'exp_datetime')
                dateStrs = {data.exp_datetime}';
                assert(iscellstr(dateStrs),'''exp_datetime'' must take string values.');
                numChars = cellfun(@numel,dateStrs);
                uniqueNumChars = unique(numChars);
                assert(numel(uniqueNumChars)==1,'Invalid format for ''exp_datetime'' field.');
                switch uniqueNumChars
                    case 8
                        obj.fExpDateFormat = 'yyyymmdd';
                    case 15
                        obj.fExpDateFormat = 'yyyymmddTHHMMSS';
                    otherwise
                        error('Invalid format for ''exp_datetime'' field.');
                end
                try
                    datenum(dateStrs,obj.fExpDateFormat);
                catch %#ok<CTCH>
                    error('Invalid format for ''exp_datetime'' field.');
                end
            end
            
            % line_name. optional field
            if isfield(data,'line_name')
                lNames = {data.line_name}';
                assert(iscellstr(lNames),'''line_name'' must take string values.');
            else
                warning('OlyDat:Browser:noLineName','Data has no ''line_name'' field.');
                [data.line_name] = deal('Unknown line');                
            end
            
            %%% BrowserStats
            
            % BStats, GroupStats, AuxStats: check fields for existence and
            % warn. If one of these are missing nothing bad will happen 
            % unless the user actually chooses to use the missing stat.
            warnst = warning('backtrace','off');
            zlclThrowWarningForMissingFields(data,obj.fBStats);
            zlclThrowWarningForMissingFields(data,obj.fGroupStats(2:end)); % XXX hardcoded, first element of fGroupStats is the "no group" stat BUT NOT FOR AGGRESSION.
            %zlclThrowWarningForMissingFields(data,obj.fAuxStats);
            warning(warnst);
            
            % ExpDetailStats: check fields for existence and add fields if
            % one is missing. If one of these is missing,
            % experiment-selection will break because selection updates the
            % expdetail pane.
            data = zlclAddMissingFields(data,obj.fExpDetailStats);            
        end
        
        function exportPlot(obj)
            % make tmp plot
            h = figure;
            ax = axes;
            obj.plotInAxes(ax,false);
            
            % get a location/file
            exportpath = ExpPP.loadConfigVal('exportpath');
            if isempty(exportpath)
              exportpath = pwd;
            end
            filterspec = fullfile(exportpath,'*.fig');
            [fname,pth] = uiputfile(filterspec,'Export figure');
            
            if isequal(fname,0) || isequal(pth,0)
              % none
            else
              figname = fullfile(pth,fname);
              hgsave(h,figname);            
              ExpPP.saveConfigVal('exportpath',pth);
            end
            delete(h);
        end
    
        % TODO loadData button
        function loadData(obj,data)
        %loadData(obj,data)
        %   Load a dataset into the Browser.
        
            if isempty(data)
                error('OlyDat:Browser:emptyData','Data is empty.');
            end
            
            data = data(:);
            
            % Data integrity check for Browser. This will throw if
            % something is wrong with data.
            data = obj.checkDataForBrowserRequiredFields(data);
            
            obj.fData = data;
            obj.fDataTF = true(numel(data),1);
            eid = [data.experiment_id]';
            obj.fEid2Idx = containers.Map(eid,1:numel(eid));            			

            numExpStr = sprintf('%d',obj.fNumExps);
            set(obj.fBrowserGData.txNumExp,'String',numExpStr);
            
            if isfield(obj.fData,'exp_datetime')
                allDates = sort({obj.fData.exp_datetime});             
                startDate = datenum(allDates{1},obj.fExpDateFormat);
                endDate = datenum(allDates{end},obj.fExpDateFormat);
                fromToStr = sprintf('%s to %s',datestr(startDate,2),datestr(endDate,2));
                set(obj.fBrowserGData.txFromTo,'String',fromToStr);
            else
                set(obj.fBrowserGData.txFromTo,'String','Unknown dates');
            end
            
            obj.refreshNumCurrentExps();
        end
        
        function plot(obj)
        % plot(obj)
        %   Apply the currently selected plot to the currently selected
        %   data in the main browser axis, using the currently selected
        %   statistic/grouping var/etc.
            
            ax = obj.fAxes;
            
            % Whatever plot used to be in this axis is going away, so we
            % have to call the axis DeleteFcn per the agreement with
            % OlyDat.Plot.
            fcn = get(ax,'DeleteFcn');
            if ~isempty(fcn)
                feval(fcn,ax,[]); % no eventdata
            end
        
            hold(ax,'off');
            cla(ax,'reset');
            axis(ax,'auto');
             
            obj.plotInAxes(ax,true);
            
            % re-select currently selected experiments to update selection
            % in main window plot
            obj.selectExperiment(obj.fAllSelectedEids);
            
            obj.plotRefreshed();
        end
        
        % The load operation a) reset the Browser curation info and then b)
        % loads curation info from the specified file for the currently
        % loaded experiments. There will almost always be experiments
        % loaded in the Browser that are not in the file; there may be
        % experiments in the file that are not loaded in the Browser.
        function curationLoadCurationFile(obj,fname)
            if nargin < 2 || isempty(fname)
                fname = obj.curationCurationFileHelper(@uigetfile);
                if isempty(fname)
                    return;
                end
            end
            assert(ischar(fname),'Filename must be a string.');
            if exist(fname,'file')~=2
                error('OlyDat:Browser:fileNotFound',...
                    'File ''%s'' not found.',fname);
            end
            
            tfOkToContinue = obj.curationCheckForUnsavedCurationInfo();
            if ~tfOkToContinue
                return;
            end
            
            % When loading a curation file, we do not want to throw the
            % warning that says "Curation field 'blah' not present in data.
            warnst = warning('off',obj.CURATION_FIELD_NOT_IN_DATA_WID);
            obj.initCurationInfoAndMarkedThisSessionAndUnsavedCurationInfo();
            warning(warnst);
            
            expName2CurationInfo = obj.fCurationWriter.read(fname);
            
            % incorporate expName2CurationInfo into fCurationInfo
            fileExpNames = expName2CurationInfo.keys;
            for c = 1:numel(fileExpNames)
                expName = fileExpNames{c};
                if obj.fCurationExpName2ExpIDMap.isKey(expName)
                    % expName is in the current dataSet
                    eid = obj.fCurationExpName2ExpIDMap(expName);
                    existingInfo = obj.fCurationInfo(eid);
                    fileInfo = expName2CurationInfo(expName);
                    existingInfo = zlclStructMerge(existingInfo,fileInfo);
                    obj.fCurationInfo(eid) = existingInfo;
                end
            end
            
            obj.curationRefreshCurationTableFromInfo;
            obj.fCurationFile = fname;
        end
        
        % Take all curation info *marked this session* and append to an
        % existing curation file. A warning is thrown if the file already
        % contains an experiment_id that is to be written. There may be
        % unmarked experiments in the Browser for which no data is written
        % to the file; there may be curation data in the file for
        % experiments not currently loaded in the Browser.
        function curationAppendCurationDataToFile(obj,fname)
            if nargin < 2 || isempty(fname)
                fname = obj.curationCurationFileHelper(@uigetfile);
                if isempty(fname)
                    return;
                end
            end
            assert(ischar(fname),'Filename must be a string.');
            if exist(fname,'file')~=2
                error('OlyDat:Browser:fileNotFound',...
                    'File ''%s'' not found.',fname);
            end
            
            infoToSave = obj.curationGetInfoMarkedThisSession();
            if infoToSave.isempty()
                msgbox('No curation marks have been made in this session. No action will be taken.',...
                    'No action taken');
                assert(~obj.fCurationUnsavedCurationInfo);
                obj.fCurationFile = fname;
            else            
                [tfSuccess rowsWritten] = obj.fCurationWriter.append(fname,infoToSave,true);
                if tfSuccess
                    msg = sprintf('%d row(s) appended to TSV file, %d row(s) overwritten.',rowsWritten.append,rowsWritten.overwrite);
                    msgbox(msg,'Curation Marks Written');
                    obj.fCurationUnsavedCurationInfo = false;
                    obj.fCurationFile = fname;
                end
            end
        end
        
        % Save all "marked" curation info to a new curation file.
        function curationWriteCurationDataToNewFile(obj,fname)
            if nargin < 2 || isempty(fname)
                fname = obj.curationCurationFileHelper(@uiputfile);
                if isempty(fname)
                    return;
                end
            end
            assert(ischar(fname),'Filename must be a string.');
            pathstr = fileparts(fname);            
            if exist(pathstr,'dir')~=7
                error('OlyDat:Browser:dirNotFound',...
                    'Directory ''%s'' not found.',pathstr);
            end
            
            infoToSave = obj.curationGetInfoMarkedThisSession();
            [tfSuccess rowsWritten] = obj.fCurationWriter.writeNewFile(fname,infoToSave);
            if tfSuccess
                msg = sprintf('%d row(s) written to TSV file.',rowsWritten);
                msgbox(msg,'Curation Marks Written');
                obj.fCurationUnsavedCurationInfo = false;
                obj.fCurationFile = fname;
            end
        end
        
    end
    
    methods (Hidden)
        
        function addComparePlot(obj)
            obj.fComparePlots(end+1,1) = figure('CloseRequestFcn',@(src,evt)obj.rmComparePlot(src));
            obj.fComparePlots(end,2) = axes;
            obj.plotInAxes(obj.fComparePlots(end,2),false);
        end
        
        function rmComparePlot(obj,hFig)
            [tf loc] = ismember(hFig,obj.fComparePlots(:,1));
            assert(tf);
            obj.fComparePlots(loc,:) = [];
            delete(hFig);
        end
        
        function clearComparePlots(obj)
            while ~isempty(obj.fComparePlots)
                obj.rmComparePlot(obj.fComparePlots(1,1));
            end
        end
        
        % cfg.stat: Nx1 array of BrowserStat objs
        % cfg.grp: scalar BrowserStat obj
        % cfg.grpVars: col int array or cellstr of values for grouping var
        % cfg.aux: scalar BrowserStat obj
        function cfg = getPlotConfig(obj)
            gdata = obj.fBrowserGData;

            % stat
            idx = get(obj.fBStatUIControl,'Value');
            cfg.stat = obj.fBStats(idx);
            
            % grp
            idx = get(gdata.pumGroup,'Value');
            cfg.grp = obj.fGroupStats(idx);
            cfg.grpVars = obj.hlpGroupVar(gdata.lbGroupVars,cfg.grp.ValidValues);
            % grp2
            if strcmp(get(gdata.pumGroup2,'Enable'),'on')              
                idx = get(gdata.pumGroup2,'Value');
                cfg.grp2 = obj.fGroupStats(idx);
                cfg.grpVars2 = obj.hlpGroupVar(gdata.lbGroupVars2,cfg.grp2.ValidValues);
            else
                cfg.grp2 = [];
                cfg.grpVars2 = [];
            end
            
            cfg.aux = obj.fCurrentAuxStat;
        end        
    end    
    methods (Static)
        function grpvars = hlpGroupVar(lb,validvals)
            strs = get(lb,'String');
            idx = get(lb,'Value');
            if ~isempty(strs)                
                grpvars = strs(idx);
                assert(~isempty(validvals));
                if isnumeric(validvals)
                    grpvars = cellfun(@str2double,grpvars);
                end
            else
                grpvars = [];
            end          
        end
    end    
    methods        
        function plotInAxes(obj,ax,tfUpdatePlotDesc)
            if ~obj.fDataIsLoaded
                return;
            end           
            
            plotCfg = obj.getPlotConfig();
            plt = obj.fCurrentPlot;
            
            % if necessary, do custom data flattening
            data = plt.preprocessData(obj.fCurrentData,plotCfg);
            assert(numel(data)==numel(obj.fCurrentData));

            % collect stats of interest in bstat matrix; this is Nexps x
            % numStats
            if plt.AllowsMultiStats
                % multi-stat data handling: for now, each stat must have a
                % scalar value. bstat is a Nexp x Nstat matrix.                
                Nstat = numel(plotCfg.stat);
                bstat = nan(numel(data),Nstat);
                for c = 1:Nstat
                    bstObj = plotCfg.stat(c);
                    tmpfld = bstObj.getDataFieldName();
                    tmpdata = {data.(tmpfld)};
                    ismissingdata = cellfun(@isempty,tmpdata);
                    assert(all(cellfun(@(x)isscalar(x)&&(isnumeric(x)||islogical(x)),tmpdata(~ismissingdata))),...
                      'Found nonscalar/nonnumeric data in field %s.',tmpfld);
                    bstat(~ismissingdata,c) = cell2mat(tmpdata(~ismissingdata))';
                end
            else
                % single-stat data handling: each stat may be a row value
                % of arbitrary width. bstat is a Nexp x dataWidth matrix.
                assert(numel(plotCfg.stat)==1);
                
                dataFieldName = plotCfg.stat.getDataFieldName(); 
                bstat = OlyDat.Browser.extractSingleBStat(data,dataFieldName);
            end
            
            % get rid of rows that are all nans
            tfNonNanData = ~all(isnan(bstat),2);

            % grp, grpVars
            [tfGrp1,grp,tfInGroups] = obj.hlpGroupCfg(data,plotCfg,'grp','grpVars');
            [tfGrp2,grp2,tfInGroups2] = obj.hlpGroupCfg(data,plotCfg,'grp2','grpVars2');
            assert(isequal(size(data),size(grp),size(grp2),size(tfInGroups),size(tfInGroups2)));
                        
            tf = tfNonNanData & tfInGroups & tfInGroups2;
            if plt.UsesGroupVar2
              grpArg = {grp(tf) grp2(tf)};              
            else              
              grpArg = grp(tf);
            end
            
            % ViewDetail
            vd = struct();
            vd.nstat = size(bstat,2);
            vd.nobsshown = nnz(tf);
            vd.nallnanrows = nnz(~tfNonNanData);
            [vd.tbl,~,~,lbl] = crosstab(grp(tfNonNanData),grp2(tfNonNanData));            
            vd.lbl1 = lbl(1:size(vd.tbl,1),1);
            if ~tfGrp1
              assert(isscalar(vd.lbl1)&&strcmp(vd.lbl1{1},'1'));
              vd.lbl1{1} = '<all>';
            end
            vd.lbl2 = lbl(1:size(vd.tbl,2),2);
            if ~tfGrp2
              assert(isscalar(vd.lbl2)&&strcmp(vd.lbl2{1},'1'));
              vd.lbl2{1} = '<all>';
            end
            % Determine which entries in tbl are currently included in plot
            [tblshown,~,~,lblshown] = crosstab(grp(tf),grp2(tf)); 
            lbl1shown = lblshown(1:size(tblshown,1),1);
            lbl2shown = lblshown(1:size(tblshown,2),2);
            vd.lbl1tfshown = ismember(vd.lbl1,lbl1shown);
            vd.lbl2tfshown = ismember(vd.lbl2,lbl2shown);            
            
            if any(tf)
%                 try
                    str = plt.doPlot(ax,data(tf),bstat(tf,:),grpArg,plotCfg,obj.fExpCoordinator);
%                 catch ME
%                     str = '';
%                     warning('OlyDat:Browser:errDuringPlot',...
%                         'Error caught during plot: %s',ME.message);
%                 end
            else
                text(0.5,0.5,'No Data','HorizontalAlignment','center','Parent',ax);
                str = '';
            end
                        
            if tfUpdatePlotDesc
                obj.handlePlotDescription(str);
                obj.handleViewDetail(vd);
            end
            
        end
    end
    methods (Static) 
        function [tfGroup,grp,tfInGroups] = hlpGroupCfg(data,plotCfg,fldGrp,fldGrpVars)
          
          % get/verify the grouping var
          % output: grp, tfGroup
          grpStat = plotCfg.(fldGrp);
          if ~isempty(grpStat) && ~isempty(grpStat.Name)
            tfGroup = true;
            grp = {data.(grpStat.Name)}';
            if iscellstr(grp)
              % none
            else
              assert(all(cellfun(@(x)isscalar(x)&&(isnumeric(x)||islogical(x)),grp)));
              grp = cell2mat(grp);
              if ~isempty(grpStat.ValidValues)
                assert(isnumeric(grpStat.ValidValues));
              end
            end
          else
            % hardcoded special case: an empty plotCfg.grp.Name
            % indicates "no group"
            tfGroup = false;
            grp = ones(size(data));
          end
          
          % get inclusion in selected vars (if applicable)
          grpVars = plotCfg.(fldGrpVars);
          if ~isempty(grpVars)
            assert(tfGroup);
            tfInGroups = ismemberwithequalnans(grp,grpVars);
          else
            tfInGroups = true(size(grp));
          end
        
          assert(isequal(size(grp),size(tfInGroups)));
        end
    end
    
    %% Experiment selection
    methods
        
        % Select the given experiment ids.
        function selectExperiment(obj,eids)
            obj.fExpCoordinator.sendSignal(obj,eids);
        end
                    
        % Browser-only response to experiment selection (expDetail,
        % expDetailPlots). To trigger all experiment-selection-related
        % callbacks, call selectExperiment.
        function handleExpSelectionLimited(obj,eids)
            idxs = nan(size(eids));
            for c = 1:numel(eids)
                idxs(c) = obj.fEid2Idx(eids(c));
            end
            obj.handleExpSelectionLimitedIdx(idxs);            
        end
        
        function handleExpSelectionLimitedIdx(obj,idxs)
            tfNoSelection = isempty(idxs);
            
            % Sset fSelectedExpIdxs, fSelectedExpIdxsIdx
            obj.fSelectedExpIdxs = idxs;
            if tfNoSelection
                obj.fSelectedExpIdxsIdx = 0;
            else
                obj.fSelectedExpIdxsIdx = 1; % By default, start by displayed first in selected group
            end            

            % Enable/disable controls on Experiment Detail panel
            if tfNoSelection
                enableEnum = 'off';
            else
                enableEnum = 'on';
            end            
            set(obj.fBrowserGData.cbRemoveExp,'Enable',enableEnum);
            %%% FOR ADAM
%             set(obj.fBrowserGData.pbEnterCurationInfo,'Enable',enableEnum);
%             set(obj.fBrowserGData.pbCreateCustomGroup,'Enable',enableEnum);
%             set(obj.fBrowserGData.pbLineReport,'Enable',enableEnum);
            set(obj.fBrowserGData.pbExpDetail,'Enable',enableEnum);
            set(obj.fBrowserGData.pbPrevExpDetail,'Enable',enableEnum);
            set(obj.fBrowserGData.pbNextExpDetail,'Enable',enableEnum);
%             set(obj.fBrowserGData.pbFileSystem,'Enable',enableEnum);
            
            obj.refreshExpDetailViewingStr();
            obj.refreshExpDetailTable();
            obj.refreshExperimentDetailPlots();            
        end
        
    end
        
    %% Experiment Detail Pane
    
    % Detail table, detail plots, viewstr, next, prev
    methods 
        
        function refreshExpDetailViewingStr(obj)
            set(obj.fBrowserGData.edExpDetailViewDisplayIdx,'String',num2str(obj.fSelectedExpIdxsIdx));
            str = sprintf('out of %d selected',numel(obj.fSelectedExpIdxs));
            set(obj.fBrowserGData.txExpDetailViewing,'String',str);
        end
        
        function refreshExpDetailTable(obj)
            detailData = get(obj.fBrowserGData.tblExpDetail,'Data'); 
            if ~isempty(obj.fSelectedExpIdxs)
                validateattributes(obj.fSelectedExpIdxsIdx,{'numeric'},...
                    {'integer' 'positive' '<=' numel(obj.fSelectedExpIdxs)});
                
                displayIdx = obj.fSelectedDisplayedIdx;
                assert(~isempty(displayIdx));
                dat = obj.fData(displayIdx);
                for c = 1:obj.fNumExpDetailStats
                    fname = obj.fExpDetailStats(c).Name;
                    detailData{c,2} = dat.(fname);
                end
                set(obj.fBrowserGData.tblExpDetail,'Data',detailData);
                set(obj.fBrowserGData.cbRemoveExp,'Value',~obj.fDataTF(displayIdx));
            else
                detailData(:,2) = {''};
                set(obj.fBrowserGData.tblExpDetail,'Data',detailData);
            end 
        end
        
        function openLineReport(obj)
            if isempty(obj.fAssayName)
                warndlg('Assay name is unspecified; cannot open line report.',...
                    'Unknown Assay Name');
                return;
            end
            
            if isempty(obj.fAllSelectedData)
                selectedLines = cell(0,1);
            else
                selectedLines = {obj.fAllSelectedData.line_name}';
                % selectedLines guaranteed to be cellstr
            end
            
            tfOkLineName = cellfun(@(x)~strcmp(x,'Unknown line'),selectedLines);
            selectedLines = selectedLines(tfOkLineName);
            if isempty(selectedLines)
                warndlg('No (valid) lines in current selection.','No Lines Selected.');
                return;
            end
            
            url = OlyDat.lineReportURL(obj.fAssayName,selectedLines);
            web(url,'-browser');
        end
        
        function openExperimentDetailPlots(obj)
            if obj.fPassAllSelectedExpsToDetailHandler
                obj.fExpDetailHandler.open(obj.fAllSelectedData);
            else
                obj.fExpDetailHandler.open(obj.fSelectedDisplayedData);
            end
        end
                
        function refreshExperimentDetailPlots(obj)
            if obj.fPassAllSelectedExpsToDetailHandler
                obj.fExpDetailHandler.refresh(obj.fAllSelectedData);
            else
                obj.fExpDetailHandler.refresh(obj.fSelectedDisplayedData); 
            end
        end
        
        function prevExpDetail(obj)
            if obj.fSelectedExpIdxsIdx==1
                % none
            else
                obj.fSelectedExpIdxsIdx = obj.fSelectedExpIdxsIdx - 1;
                obj.refreshExpDetailViewingStr();
                obj.refreshExpDetailTable();
                obj.refreshExperimentDetailPlots();
            end            
        end
        
        function nextExpDetail(obj)
            if obj.fSelectedExpIdxsIdx==numel(obj.fSelectedExpIdxs)
                % none
            else
                obj.fSelectedExpIdxsIdx = obj.fSelectedExpIdxsIdx + 1;
                obj.refreshExpDetailViewingStr();
                obj.refreshExpDetailTable();
                obj.refreshExperimentDetailPlots();
            end                        
        end
        
        % Create a custom group out of the currently selected experiments.
        function createCustomGroup(obj,grpname)
            
            if ~isvarname(grpname)
                newGrpName = genvarname(grpname);
                warning('Browser:nonVarNameGrpName',...
                    'Group name ''%s'' is not a valid MATLAB variable name. Using modified name ''%s''.',...
                    grpname,newGrpName);
                grpname = newGrpName;
            end
            
            % Validate fieldname for data struct
            allCustomGroupNames = obj.fAllCustomGroupNames;
            if ismember(grpname,allCustomGroupNames)
                error('Browser:cantCreateCustomGroup',...
                    'A custom group with name ''%s'' already exists.',grpname);
            end
            customGroupFieldName = obj.mangleGroupName(grpname);
            if isfield(obj.fData,customGroupFieldName)
                % one in a kazillion
                error('Browser:cantCreateCustomGroup',...
                    'The field ''%s'' already exists in the data.',grpname);
            end
            
            % Add new field to data
            tfInGroup = false(numel(obj.fData),1);
            tfInGroup(obj.fSelectedExpIdxs) = true;
            [obj.fData(tfInGroup).(customGroupFieldName)] = deal('in group');
            [obj.fData(~tfInGroup).(customGroupFieldName)] = deal('out of group');
            
            % Add BrowserStat
            groupStatClass = class(obj.fGroupStats); 
            
            % assume accepted constructor sig: ctor(name,prettyName)
            newGroupBrowserStat = feval(groupStatClass,customGroupFieldName,...
                sprintf('Custom group: %s',grpname));
            newGroupBrowserStat.ValidValues = {'in group';'out of group'};
            obj.fGroupStats(end+1,1) = newGroupBrowserStat;
            set(obj.fBrowserGData.pumGroup,'String',{obj.fGroupStats.PrettyName}');
            obj.refreshGroupsList();
        end
        
        % Open filesystem for currently selected experiment.
        function openExperimentFileSystem(obj)
            dat = obj.fSelectedDisplayedData;
            if ~isempty(dat) && isfield(dat,'file_system_path') && ~isempty(dat.file_system_path)
                
                if isempty(obj.fExperimentFileSystemFolder)
                    expPath = uigetdir([],'Specify folder containing experiments:');
                    if isequal(expPath,0)
                        return;
                    end
                    obj.fExperimentFileSystemFolder = expPath;
                end                                    
                
                % Strip off exp folder and add to ExperimentFileSysFolder
                fsp = dat.file_system_path;
                [~,expFlder] = fileparts(fsp);
                if isempty(expFlder)
                    warndlg(sprintf('Corrupt file system path ''%s''. Aborting.',fsp));
                    return;
                end
                fullfsp = fullfile(obj.fExperimentFileSystemFolder,expFlder);
                    
                while exist(fullfsp,'dir')==0
                    queststr = sprintf('Cannot find location ''%s''. Do you want to specify a new experiment folder?',...
                        fullfsp);
                    resp = questdlg(queststr,'Experiment folder inaccessible','Yes','Cancel','Yes');
                    switch resp
                        case 'Cancel'
                            return;
                    end
                    expPath = uigetdir([],'Specify folder containing experiments');
                    if isequal(expPath,0)
                        return;
                    end
                    obj.fExperimentFileSystemFolder = expPath;
                    fullfsp = fullfile(obj.fExperimentFileSystemFolder,expFlder);
                end
                
                if ispc
                    winopen(fullfsp);                    
                elseif ismac
                    system(['open ' fullfsp]);
                elseif isunix
                    web(fullfsp);                    
                end
            else
                warndlg('No experiment selected, or missing file_system_path.',...
                    'Can''t open experiment folder');
            end
        end
        
    end
    
    % Include/exclude
    methods
        
        % Exclude/include all currently selected experiments if tf is
        % true/false.
        function setRemoveExp(obj,tf,tfUpdateCurationMgr)
            if nargin < 3
                tfUpdateCurationMgr = true;
            end
            assert(~isempty(obj.fSelectedExpIdxs));
            
            % exclude state
            obj.fDataTF(obj.fSelectedExpIdxs) = ~tf;
            
            % DataSet panel
            obj.refreshNumCurrentExps;
            
            % Experiment Detail checkbox
            set(obj.fBrowserGData.cbRemoveExp,'Value',tf);
            
            % Curation Manager
            if tfUpdateCurationMgr
                eids = obj.fAllSelectedEids;
                tblData = get(obj.fCurationGData.tblCuration,'Data');
                tblEids = cell2mat(tblData(:,obj.fCurationGData.tableMetadataCurationField2Col.experiment_id));
                
                tfSelected = ismember(tblEids,eids);
                tblData(tfSelected,obj.fCurationGData.tableMetadataExcludedColIdx) = {logical(tf)};
                set(obj.fCurationGData.tblCuration,'Data',tblData);
            end
            
            % plot
            obj.plot();
            
            % resend experiment selection signal as appropriate
            % (not sure why I have to do this)
            if ~isempty(obj.fSelectedExpIdxs)
                eids = obj.fAllSelectedEids;
                obj.selectExperiment(eids);
            end
        end
        
        function unexcludeAllExps(obj)
            obj.fDataTF = true(obj.fNumExps,1);
            obj.refreshNumCurrentExps;
            obj.refreshExpDetailTable();
            curationData = get(obj.fCurationGData.tblCuration,'Data');
            curationData(:,obj.fCurationGData.tableMetadataExcludedColIdx) = {false};
            set(obj.fCurationGData.tblCuration,'Data',curationData);
            obj.plot;
        end
        
    end
    
    %% Data Set Pane
    methods
        
        function refreshNumCurrentExps(obj)
            str = sprintf('%d',obj.fNumCurrentExps);
            set(obj.fBrowserGData.txNumCurrentExps,'String',str);
        end
        
    end
    
    %% Misc
    methods
        
        function handlePlotDescription(obj,str)
          % Do some Yair voodoo to enable horizontal scrollbar. Note:
          % resizing the table (eg when figure is resized) will cause
          % horizontal scrollbar to go away. Refreshing plot will bring it
          % back.          
          pause(0.01);
          drawnow; % on occassion, getting "no method setHorizontalScrollBarPolicy on class handle.handle"
          try
            jScrollPane = findjobj(obj.fBrowserGData.etPlotDescription);         
            % Modify the scroll-pane's scrollbar policies
            %set(jScrollPane,'VerticalScrollBarPolicy',20);  % or: jScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED
            %jScrollPane.setHorizontalScrollBarPolicy(30);  % or: jScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED
            jViewPort = jScrollPane.getViewport;
            jEditbox = jViewPort.getComponent(0);
            jEditbox.setWrapping(false); % do *NOT* use set(...)!!!
            newPolicy = 30;% jScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED;
            set(jScrollPane,'HorizontalScrollBarPolicy',newPolicy);          
          catch ME
            warning('Browser:warning','Plot description horizontal scrollbar issues: %s.',ME.message);
          end            
          
          set(obj.fBrowserGData.etPlotDescription,'String',str);
        end        
        
        function handleViewDetail(obj,vdinfo)   
            T = vdinfo.tbl;
            T(:,end+1) = sum(T,2);
            T(end+1,:) = sum(T,1);
            T = num2str(T);
            
            lbl1str = vdinfo.lbl1;
            for i = 1:numel(lbl1str)
              if vdinfo.lbl1tfshown(i), lbl1str{i} = [lbl1str{i} '*']; end
            end
            lbl1str = [lbl1str;{'TOT'}];
            lbl1str = strcat(lbl1str,{' '});
            T = [char(lbl1str) T];
            
            lbl2str = vdinfo.lbl2;            
            for i = 1:numel(lbl2str)
              if vdinfo.lbl2tfshown(i), lbl2str{i} = [lbl2str{i} '*']; end
            end
            lbl2str{end+1} = 'TOT';
            lbl2str = sprintf('%s,',lbl2str{:});
            lbl2str = ['  ' lbl2str(1:end-1)];
            
            str = cell(0,1);
            str{end+1,1} = sprintf('Nstat: %d. Nallnan: %d. Nshown: %d.',...
              vdinfo.nstat,vdinfo.nallnanrows,vdinfo.nobsshown);
            str{end+1,1} = 'Contingency Table';
            str{end+1,1} = 'Columns:';
            str{end+1,1} = lbl2str;
            str = [str;cellstr(T)];
            str{end+1,1} = '(*) Group currently included';
            
            % See comments in handlePlotDescription
            pause(0.01);
            drawnow;
            try
              jScrollPane = findjobj(obj.fBrowserGData.etViewDetail);
              jViewPort = jScrollPane.getViewport;
              jEditbox = jViewPort.getComponent(0);
              jEditbox.setWrapping(false); 
              newPolicy = 30;
              set(jScrollPane,'HorizontalScrollBarPolicy',newPolicy);
            catch ME
              warning('Browser:warning','View detail horizontal scrollbar issues: %s.',ME.message);
            end
            
            set(obj.fBrowserGData.etViewDetail,'String',str);
        end
        
        function statChanged(obj)
            plt = obj.fCurrentPlot;
            if plt.AutoReplotAfterStatChange
              obj.plot();
            else
              obj.plotNeedsRefresh();
            end
        end
        
        function groupChanged(obj)
            obj.refreshGroupsList();
            obj.callPlotGroupChanged();
            plt = obj.fCurrentPlot;
            if plt.AutoReplotAfterGroupChange
              obj.plot();
            else
              obj.plotNeedsRefresh();
            end
        end
        function group2Changed(obj)
            obj.refreshGroupsList();
            plt = obj.fCurrentPlot;
            if plt.AutoReplotAfterGroupChange
              obj.plot();
            else
              obj.plotNeedsRefresh();
            end
        end       
        function groupOrGroup2SelectionChanged(obj)
            plt = obj.fCurrentPlot;
            if plt.AutoReplotAfterGroupChange
              obj.plot();
            else
              obj.plotNeedsRefresh();
            end          
        end
        
        function plotNeedsRefresh(obj)
          set(obj.fBrowserGData.pbPlot,'foregroundcolor',[1 0 0],'fontweight','bold');
        end
        function plotRefreshed(obj)
          set(obj.fBrowserGData.pbPlot,'foregroundcolor',[0 0 0],'fontweight','normal');
        end
        
        function callPlotGroupChanged(obj)
            currentPlot = obj.fCurrentPlot;
            grpIdx = get(obj.fBrowserGData.pumGroup,'Value');
            grpStat = obj.fGroupStats(grpIdx);                
            hCtrls = struct();
            hCtrls.pumAuxVar = obj.fBrowserGData.pumAuxVar;
            currentPlot.groupChanged(grpStat,obj.fData,hCtrls);            
        end
        
        function refreshGroupsList(obj)          
          gSets = struct('pum',{'pumGroup' 'pumGroup2'},'lb',{'lbGroupVars' 'lbGroupVars2'});
          for gIdx = [1 2]
            pumTag = gSets(gIdx).pum;
            lbTag = gSets(gIdx).lb;
            
            grpIdx = get(obj.fBrowserGData.(pumTag),'Value');
            grp = obj.fGroupStats(grpIdx);
            groupVals = grp.ValidValues;
            if ~isempty(groupVals) && isnumeric(groupVals);
                groupVals = num2cell(groupVals);
                groupVals = cellfun(@num2str,groupVals,'UniformOutput',false);
            end
            
            previousGroupVals = get(obj.fBrowserGData.(lbTag),'String');
            if isequal(previousGroupVals,groupVals)
                % none
            else
                set(obj.fBrowserGData.(lbTag),'String',groupVals,'Value',1:numel(groupVals));
            end
          end
        end
	
        function newPlot(obj,idx)
            if nargin < 2
                idx = get(obj.fBrowserGData.pumPlotType,'Value');
            end
            plot = obj.fPlots{idx};

            %% Manage BStat uicontrol
            tfLastPlotAllowedMultiStats = strcmp(get(obj.fBrowserGData.lbBStats,'Visible'),'on');
            if plot.AllowsMultiStats==tfLastPlotAllowedMultiStats
                % previous plot and current plot agree on allowing
                % multistats
            elseif plot.AllowsMultiStats
                % previous plot didn't allow multistat, new one does
                currentBStat = get(obj.fBrowserGData.pumBStats,'Value');
                set(obj.fBrowserGData.lbBStats,'Value',currentBStat);
                
                set(obj.fBrowserGData.lbBStats,'Visible','on');
                set(obj.fBrowserGData.pumBStats,'Visible','off');
            else
                % previous plot allowed multistat, new one doesn't
                currentBStat = get(obj.fBrowserGData.lbBStats,'Value');
                set(obj.fBrowserGData.pumBStats,'Value',currentBStat(1));
                        
                set(obj.fBrowserGData.pumBStats,'Visible','on');
                set(obj.fBrowserGData.lbBStats,'Visible','off');
            end
            
            if isequal(plot.UsesAuxVar,true) || isequal(plot.UsesAuxVar,'custom')
                set(obj.fBrowserGData.pumAuxVar,'Enable','on');
            elseif isequal(plot.UsesAuxVar,false)
                set(obj.fBrowserGData.pumAuxVar,'Enable','off');
            else
                assert(false,'Invalid value of UsesAuxVar.');
            end
            
            if plot.UsesGroupVar2==true
                tmpEnable = 'on';
            else
                tmpEnable = 'off';
            end            
            set(obj.fBrowserGData.pumGroup2,'Enable',tmpEnable);
            set(obj.fBrowserGData.lbGroupVars2,'Enable',tmpEnable);
            set(obj.fBrowserGData.pbUseAllGroupVals2,'Enable',tmpEnable);
            
            set(obj.fBrowserGData.txGroupVar,'String',plot.GroupVarLabel);
            set(obj.fBrowserGData.txGroupVar2,'String',plot.GroupVar2Label);
            
            obj.fSelectedPlot = idx;
            
            obj.callPlotGroupChanged();
            obj.callPlotAuxVarChanged();

            obj.plot();
        end
        
        function callPlotAuxVarChanged(obj)
            auxStat = obj.fCurrentAuxStat;            
            crntPlot = obj.fCurrentPlot;
            hCtrls = struct();
            hCtrls.pumAuxVar = obj.fBrowserGData.pumAuxVar;
            crntPlot.auxVarChanged(auxStat,obj.fData,hCtrls);            
        end
        
        function auxVarChanged(obj)
            obj.callPlotAuxVarChanged();            
            obj.plot();
        end
        
        function browserGUIWantsToClose(obj)
            tfOkToContinue = obj.curationCheckForUnsavedCurationInfo();
            if tfOkToContinue
                delete(obj);
            end
        end
        
    end

    %% Curation
    methods
        
        % Initialize fCurationInfo prop based on fData. Both read-only (eg
        % experiment-identifying) and writeable fields (eg curation marks)
        % are read from fData.
        function initCurationInfoAndMarkedThisSessionAndUnsavedCurationInfo(obj)
            
            % init fCurationInfo to contain empty structs with the right fields
            fields = obj.fCurationFields;
            Nfields = numel(fields);
            data = cell(numel(obj.fData),Nfields); 
            curationInfoStructs = cell2struct(data,fields,2);
            curationInfoStructs = num2cell(curationInfoStructs);
            eids = {obj.fData.experiment_id}';
            obj.fCurationInfo = containers.Map(eids,curationInfoStructs);
            
            % init fCurationInfo to match what is in fData
            for c = 1:Nfields
                fld = fields{c};
                if ~isfield(obj.fData,fld)
                  % AL20150127, Don't bother with this warning
%                     warnst = warning('backtrace','off');
%                     warning(obj.CURATION_FIELD_NOT_IN_DATA_WID,...
%                         'Curation field ''%s'' not present in dataset. Initializing this field to empty value.',fld);
%                     warning(warnst);
                end
            end
            obj.curationResetCurationInfoToData(cell2mat(eids));
            
            obj.fCurationMarkedThisSession = containers.Map(eids,repmat({false},size(eids)));
            obj.fCurationUnsavedCurationInfo = false;
        end
        
        function initCurationWriter(obj)
            [tf loc] = ismember(obj.fCurationFieldObjsTSVOrder,obj.fCurationFields);
            assert(all(tf));
            tsvObjs = obj.fCurationFieldObjs(loc);
            tsvNames = {tsvObjs.TSVName}';
            assert(all(cellfun(@(x)~isempty(x),tsvNames)));
            obj.fCurationWriter = OlyDat.CurationFileWriter(tsvObjs);
        end
        
        function initCurationExpName2ExpIDMap(obj)
            m = containers.Map('KeyType','char','ValueType','double');
            Nexps = obj.fNumExps;
            fdata = obj.fData;
            assert(Nexps>0);
            assert(Nexps==numel(fdata));
            for c = 1:Nexps
                eid = fdata(c).experiment_id;
                expname = fdata(c).experiment_name;
                assert(~m.isKey(eid));
                m(expname) = eid;
            end
            obj.fCurationExpName2ExpIDMap = m;
        end
        
        function initCurationTable(obj)
            
            % construct/open BrowserCurationGUI
            [tf loc] = ismember(obj.fCurationFieldObjsTblOrder,obj.fCurationFields);
            assert(all(tf));
            tblObjs = obj.fCurationFieldObjs(loc);
            tblNames = {tblObjs.TableName}';
            assert(all(cellfun(@(x)~isempty(x),tblNames)));
            obj.fCurationH = BrowserCurationGUI(obj,tblObjs);
            set(obj.fCurationH,'CloseRequestFcn',@(src,evt)set(src,'Visible','off'));
            obj.fCurationGData = guidata(obj.fCurationH);

            % initialize table contents in order of fCurationInfo
            assert(~isempty(obj.fCurationInfo));
            curationInfoStructs = cell2mat(obj.fCurationInfo.values);
            assert(numel(curationInfoStructs)==obj.fNumExps);

            % init table data
            tblData = cell(numel(curationInfoStructs),numel(obj.fCurationGData.tableMetadataCol2Field));
                
            % EID
            eids = {curationInfoStructs.experiment_id}';
            colIdx = obj.fCurationGData.tableMetadataCurationField2Col.experiment_id;
            tblData(:,colIdx) = eids;
            
            % initialize excluded col
            assert(all(obj.fDataTF),'Expect all experiments to be unexcluded at this time.');
            tblData(:,obj.fCurationGData.tableMetadataExcludedColIdx) = {false};
            
            set(obj.fCurationGData.tblCuration,'Data',tblData);
            
            % set up rest of table
            obj.curationRefreshCurationTableFromInfo();
        end
        
        function fname = curationCurationFileHelper(obj,uiFileFcn)
            if ~isempty(obj.fCurationFile)
                hintFile = obj.fCurationFile;
            else
                hintFile = pwd;
            end
            [fname pname] = uiFileFcn(OlyDat.CurationFileWriter.fileExtension,...
                'Select curation file',hintFile);
            if isnumeric(fname)
                % user cancelled
                fname  = [];
            else
                fname = fullfile(pname,fname);
            end
        end
        
        function curationOpenCurationManager(obj)
            set(obj.fCurationH,'Visible','on');
        end
        
        % Refresh curation table information using fCurationInfo (exclude
        % column unaffected). This is done by looking up info for each
        % experiment_id, in table order. Row-order of table is unchanged.
        function curationRefreshCurationTableFromInfo(obj)
            tblData = get(obj.fCurationGData.tblCuration,'Data');
            info = obj.fCurationInfo;
            Nrows = size(tblData,1);
            assert(Nrows==info.size(1)); % maybe unnecessary in future when tbl shows restricted view
            
            fldnames = fieldnames(obj.fCurationGData.tableMetadataCurationField2Col);
            
            for c = 1:Nrows
                row = tblData(c,:);
                eid = row{obj.fCurationGData.tableMetadataCurationField2Col.experiment_id};
                assert(info.isKey(eid));
                s = info(eid);
                for d = 1:numel(fldnames)
                    % custom conversion for table
                    fld = fldnames{d};
                    val = s.(fld);
%                     switch fld
%                         case 'manual_pf'
%                             if isempty(val)
%                                 val = ' '; % for empty option in popupmenu. Is there a workaround for this?
%                             else
%                                 assert(any(strcmp(val,{'U';'P';'F'})));
%                             end
%                         case {'flag_redo' 'flag_review'} % xxx other flag fields?
%                             if isempty(val)
%                                 val = ' ';
%                             else
%                                 assert(any(strcmp(val,{'0';'1'})));
%                             end                            
%                         otherwise
%                             %none
%                     end
                    row{obj.fCurationGData.tableMetadataCurationField2Col.(fld)} = val;
                end
                tblData(c,:) = row;
            end
            
            set(obj.fCurationGData.tblCuration,'Data',tblData);            
        end
                
        % Methods that enable editing of fCurationInfo:
        % * curationLoadCurationFile
        % * curationUpdateCurationData
        % * curationResetCurationInfoToData
        % * curationEditCurationInfoForCurrentExp
        % * curationClearCurationInfoForCurrentExp
        % * curationMarkAllUncuratedAsPass

        % Update fCurationInfo based on infoStruct.
        function curationUpdateCurationData(obj,infoStruct)
            disp('CURRENTLY UNUSED METHOD');
            eid = infoStruct.experiment_id;
            assert(obj.fCurationInfo.isKey(eid));
            existingInfo = obj.fCurationInfo(eid);
            existingInfo = zlclStructMerge(existingInfo,infoStruct);
            obj.fCurationInfo(eid) = existingInfo;
            obj.fCurationMarkedThisSession(eid) = true;
            obj.fCurationUnsavedCurationInfo = true;
        end
        
        % reset the curation info for the experiment ID eids to what is
        % present in fData.
        function curationResetCurationInfoToData(obj,eids)
            
            % ML pre-memory mgr enhancement, this propget is amazingly slow
            % for large fData. So get this once, outside loop.
            allData = obj.fData;
            
            for id = eids(:)'
                info = obj.fCurationInfo(id);
                idx = obj.fEid2Idx(id);
                dat = allData(idx);

                flds = fieldnames(info);
                for c = 1:numel(flds)
                    fld = flds{c};                
                    if isfield(dat,fld)
                        info.(fld) = dat.(fld);
                    else
                        info.(fld) = '';
                    end
                end

                obj.fCurationInfo(id) = info;
            end
        end

        % edit curation fields for all currently selected exps.
        % s: struct with curation fields
        function curationEditCurationInfoForCurrentExp(obj,s)
            assert(isstruct(s));
            assert(~isempty(obj.fSelectedExpIdxs));
            
            assert(~any(isfield(s,obj.fCurationReadonlyFields)));
                        
            eids = obj.fAllSelectedEids;
            for c = 1:numel(eids)
                id = eids(c);
                info = obj.fCurationInfo(id);
                info = zlclStructMerge(info,s);
                obj.fCurationInfo(id) = info;
                obj.fCurationMarkedThisSession(id) = true;
            end
            
            obj.fCurationUnsavedCurationInfo = true;
            obj.curationRefreshCurationTableFromInfo();
        end
        
        % Warning: this method also resets fCurationMarkedThisSession for
        % the current experiment(s) to false. It is as if the experiment was
        % never touched. However, note that if a curation file was loaded
        % at some point, this method does NOT restore the loaded info; the
        % info is initialized based on fData. This may be a little weird but it's
        % too edgy to worry about in the absence of real usage information.
        function curationClearCurationInfoForCurrentExp(obj)
            assert(~isempty(obj.fSelectedExpIdxs));
            
            eids = obj.fAllSelectedEids;
            obj.curationResetCurationInfoToData(eids);
            for id = eids(:)'
                obj.fCurationMarkedThisSession(id) = false; % take note
            end
            
            %obj.fCurationUnsavedCurationInfo = true; Don't touch
            %fCurationUnsavedCurationInfo.
            obj.curationRefreshCurationTableFromInfo();
        end
                
        % Open curation info entry modal dialog, prepopulated with
        % currently selected exps.
        function curationOpenCurationInfoEntry(obj)
            eids = obj.fAllSelectedEids;
            for c = numel(eids):-1:1
                info(c) = obj.fCurationInfo(eids(c));
            end
            cbkFcn = @(s)obj.curationEditCurationInfoForCurrentExp(s);
            feval(obj.fCurationInfoEntryGUIName,info,obj.fUsername,cbkFcn);
        end        
        
        % If there is unsaved curation info, ask user if it is okay to
        % continue (with an operation which will cause them to lose that
        % info).
        function tfOkToContinue = curationCheckForUnsavedCurationInfo(obj)
            tfOkToContinue = true;
            if obj.fCurationUnsavedCurationInfo
                questStr = 'You have unsaved manual curation edits. You will lose those edits if you continue. Would you like to continue anyway?';
                titleStr = 'Unsaved curation data';
                btn = questdlg(questStr,titleStr,'Continue anyway','No','No');
                switch btn
                    case {'' 'No'}
                        tfOkToContinue = false;
                        return;
                    case 'Continue anyway'
                        % none
                end
            end
        end    
        
        % Get curation info marked this session.
        function expName2CurationInfo = curationGetInfoMarkedThisSession(obj)
            allEids = cell2mat(obj.fCurationMarkedThisSession.keys);
            tfMarked = cell2mat(obj.fCurationMarkedThisSession.values);
            markedEids = allEids(tfMarked);
            
            expName2CurationInfo = containers.Map('KeyType','char','ValueType','any');
            for eid = markedEids(:)'
                info = obj.fCurationInfo(eid);
                expname = info.experiment_name;
                assert(~expName2CurationInfo.isKey(expname));
                expName2CurationInfo(expname) = info;
            end
        end
        
% CURRENTLY UNUSED        
%         function expName2CurationInfo = curationGetAllMarkedInfo(obj)
%             expName2CurationInfo = containers.Map('KeyType','char','ValueType','any');
%             allEids = cell2mat(obj.fCurationInfo.keys);
%             for eid = allEids(:)'
%                 info = obj.fCurationInfo(eid);
%                 expname = info.experiment_name;
%                 if obj.curationInfoIsManualMarked(info)
%                     assert(~expName2CurationInfo.isKey(expname));
%                     expName2CurationInfo(expname) = info;
%                 end
%             end
%         end
        
        function curationMarkAllUncuratedAsPass(obj)
            curationInfo = obj.fCurationInfo;
            eids = cell2mat(curationInfo.keys);
            for id = eids(:)'
                info = curationInfo(id);
                if ~obj.curationInfoIsManualMarked(info)
                    info = obj.curationInfoMarkAsManualPass(info);
                    curationInfo(id) = info;
                    obj.fCurationMarkedThisSession(id) = true;
                end
            end
            obj.curationRefreshCurationTableFromInfo();
            obj.fCurationUnsavedCurationInfo = true;
        end
                
        % For the moment this hardcodes against fields.
        function tf = curationInfoIsManualMarked(obj,info) %#ok<MANU>
            tf = ~isempty(info.manual_curator) || ...
                 ~isempty(info.manual_curation_date);             
        end
        
        % For the moment this hardcodes against fields.
        function info = curationInfoMarkAsManualPass(obj,info)
           info.manual_pf = 'P';
           info.manual_curator = obj.fUsername;
           info.manual_curation_date = datestr(now,30);
        end        
        
    end
    
    methods (Static)
        
        function name = mangleGroupName(grpname)
            name = sprintf('Custom_Group_%s',grpname);            
        end
        
    end
    
end

function s = zlclStructMerge(s,news)
newfields = fieldnames(news);
for c = 1:numel(newfields)
    fld = newfields{c};
    s.(fld) = news.(fld);
end
end

function zlclThrowWarningForMissingFields(data,browserStatObjs)
tfBuried = [browserStatObjs.tfBuriedStat]';
fldnames = {browserStatObjs.Name}';
fldnames = fldnames(~tfBuried); % don't bother checking score stats
for c = 1:numel(fldnames)
    if ~isfield(data,fldnames{c})
        warning('OlyDat:Browser:missingField',...
            'Specified Browser stat ''%s'' does not have a corresponding field in the data. Selection of this statistic may cause an error.',...
            fldnames{c});
    end
end
end

function data = zlclAddMissingFields(data,browserStatObjs)
tfBuried = [browserStatObjs.tfBuriedStat]';
fldnames = {browserStatObjs.Name}';
fldnames = fldnames(~tfBuried); % don't bother checking score stats
for c = 1:numel(fldnames)
    if ~isfield(data,fldnames{c})
        warning('OlyDat:Browser:missingField',...
            'Specified experiment detail statistic ''%s'' is not a field in the data. Adding field with values ''<unknown>''.',...
            fldnames{c});
        [data.(fldnames{c})] = deal('<unknown>');
    end
end
end

% function v = zlclInitFO_FILE_SYSTEM_ROOT_DEFAULT
% if ispc
%     v = '\\dm10\flyolympiad\';
% elseif ismac
%     v = '/Volumes/flyolympiad/'; % typical mount location (?)
%     %v = 'smb://dm10/flyolympiad';
% elseif isunix
%     % typical mount location?
%     v = '/groups/sciserv/flyolympiad/';
% else
%     assert(false);
% end
% end

