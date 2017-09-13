classdef DataSelector < handle
%DataSelector Select olympiad data from SAGE.
%   The DataSelector lets you build up a SAGE query in a graphical table.
%   Queries may be saved for later re-use. Once a query has been
%   constructed, running it pulls data from SAGE to the base workspace.
%   Assay-specific dataSet merging, data cleaning, etc is performed as
%   appropriate.
%
%   See also OlyDat.DataSelector/DataSelector, OlyDat.Browser, OlyDat.Analyzer.

    properties (Constant)
        ComparePairs = {...
            'is equal to'     '='; ...
            'is not equal to' '!='; ...
            'is greater than' '>'; ...
            'is less than'    '<';...
            };
    end
        
    properties (Hidden)
        fDataSetFamily;  % scalar SAGE.dataSetFamily
        fDataSetFields;  % col vec of SAGE.dataSetFields
        
        fGuiH; % handle to DataSelectorGUI
        fGData; % handle to guidata struct
        fMTable; % handle to mtable, see createTable
        fJTable; % handle to jtable, see createTable
        
        fCurrentFieldHasValidValues; % scalar logical. If true, the current field (in the "on deck" query) has validValues.
        fCurrentQuerySet; % full filename for currently loaded query
        
        fDataPuller; % OlyDat.DataPuller object
        fDataPullerHasOptions; % logical scalar
    end
    
    properties (Dependent)
        AssayName;
    end
    
    methods
        function v = get.AssayName(obj)
            v = obj.fDataSetFamily.displayName;
        end
    end
    
    methods
        
        function guiWillClose(obj)
            delete(obj);
        end
        
        function delete(obj)
            obj.fMTable = [];
            obj.fJTable = [];
            delete(obj.fGuiH);          
        end
        
        function obj = DataSelector(datasetfamily,datapuller,parentH,pos) %#ok<INUSD>
        %DataSelector Construct a DataSelector and bring up the data selector GUI.
        %   obj = DataSelector(datasetfamily,datapuller) constructs a data
        %   selector. datasetfamily is a scalar SAGE.dataSetFamily object.
        %   datapuller is an OlyDat.DataPuller object.
        
            validateattributes(datasetfamily,{'SAGE.DataSetFamily'},{'scalar'});
            assert(isa(datapuller,'OlyDat.DataPuller'));
                       
            obj.fDataSetFamily = datasetfamily;
            dataSets = datasetfamily.dataSets;
            allFields = [];
            for c = 1:numel(dataSets)
                allFields = [allFields dataSets(c).fields]; %#ok<AGROW>
            end
            allFields = allFields(:);
            allFields = OlyDat.DataSelector.customizeDataSetFields(...
                                    datapuller,datasetfamily,allFields);            
            [~,idx] = unique({allFields.name});
            obj.fDataSetFields = allFields(idx);
            
            tfFS = datapuller.tfPullFSDefault;
            assert(isscalar(tfFS) && (islogical(tfFS)||isnumeric(tfFS)),...
                'DataPuller.tfPullFSDefault property must be a scalar logical.');            
            opts = datapuller.options;
            assert(iscellstr(opts),'DataPuller.options property must be a cellstr of data-pulling options.');
            obj.fDataPullerHasOptions = ~isempty(opts);
            obj.fDataPuller = datapuller;
            
            if nargin < 3
                obj.fGuiH = DataSelectorGUI(obj);
                obj.fGData = guidata(obj.fGuiH);
                set(obj.fGuiH,'CloseRequestFcn',@(src,evt)obj.guiWillClose());
            else
                error('TODO: not supported yet');
            end
                        
            obj.fMTable = createTable(obj.fGData.pnlQuery,{'Field' 'Comparison' 'Value'},cell(0,3),[1 1 0 0]);
            obj.fMTable.setGridColor([.9;.9;1.0]);
            cbFld = javax.swing.JComboBox({obj.fDataSetFields.name}');
            cbCmps = javax.swing.JComboBox(OlyDat.DataSelector.ComparePairs(:,1));
            obj.fJTable = obj.fMTable.getTable;
            obj.fJTable.getColumnModel.getColumn(0).setCellEditor(javax.swing.DefaultCellEditor(cbFld));
            obj.fJTable.getColumnModel.getColumn(1).setCellEditor(javax.swing.DefaultCellEditor(cbCmps));
            
            set(obj.fGData.pumField,'String',{obj.fDataSetFields.name}');
            set(obj.fGData.pumField,'Value',1);
            set(obj.fGData.pumCompare,'String',OlyDat.DataSelector.ComparePairs(:,1));
            obj.newFieldSelection();
            
            obj.setCurrentQuerySet('unnamed query');
        end                   
        
        function qObjs = getQueryObjs(obj)
        %qObjs = getQueryObjs(obj)
        %   Get cell array of SAGE.Query objects from current GUI state. 
        
            dat = cell(obj.fMTable.getData);
            if isempty(dat)
                qObjs = cell(0,1);
                return;
            end
            
            dat(:,2) = cellfun(@(x)obj.cmpTxt2Symbol(x),dat(:,2),'UniformOutput',false);
            Nqry = size(dat,1);
            qObjs = cell(Nqry,1);
            for c = 1:Nqry
                qObjs{c} = obj.queryObj(dat{c,:});
            end
        end
        
        function runQuery(obj)
        %runQuery(obj)
        %   Run the current query and pull data to the base workspace.
        
            % get query objs
            qobjs = obj.getQueryObjs;
            tfEmpty = cellfun(@isempty,qobjs);
            if any(tfEmpty)
                warning('DataSelector:invalidQueries',...
                    'Ignoring one or more invalid query clauses.');
            end
            qobjs = qobjs(~tfEmpty);
            
            % get data var
            dataVarName = get(obj.fGData.etDataVar,'String');
            if ~isvarname(dataVarName)
                error('DataSelector:invalidDataVariableName',...
                    'Data variable must be a valid MATLAB variable name.');
            end
            
            if isa(obj.fDataPuller,'OlyDat.DataPuller')
                if obj.fDataPullerHasOptions
                    val = get(obj.fGData.pumPullerOptions,'Value');
                    str = get(obj.fGData.pumPullerOptions,'String');
                    option = str{val};
                else
                    option = [];
                end
                tfPullFS = get(obj.fGData.cbPullFileSys,'Value');
                if tfPullFS
                    datapath = get(obj.fGData.etDataDir,'String');
                    if ~isempty(datapath) && exist(datapath,'dir')~=7
                        warning('DataSelector:badDataPath',...
                            'Experiment directory ''%s'' not found.',datapath);
                    end
                    data = obj.fDataPuller.pullDataFS(qobjs,datapath,option);
                else
                    data = obj.fDataPuller.pullData(qobjs,option);
                end
            else           
                assert(false,'Unsupported codepath.');
            end
            
            if numel(data) > 0
                msgbox(sprintf('Data for %d experiments pulled.',numel(data)),'Success');
                tfVarExists = evalin('base',sprintf('exist(''%s'',''var'')',dataVarName));
                if tfVarExists
                    warning('Analyzer:overWriteVar','Overwriting variable ''%s''.',dataVarName);
                end
                assignin('base',dataVarName,data);
            else
                msgbox('Data for 0 experiments pulled.','No Data Pulled'); 
            end
            
        end

    end
    
    methods (Hidden)
        
        function setJTableData(obj,dat)
            obj.fMTable.setData(dat);
        end                        
        
        function newFieldSelection(obj)
            idx = get(obj.fGData.pumField,'Value');
            dsf = obj.fDataSetFields(idx);
            if ~isempty(dsf.validValues)
                obj.fCurrentFieldHasValidValues = true;
            else
                obj.fCurrentFieldHasValidValues = false;
            end
            obj.newCompareSelection();
        end
        
        function newCompareSelection(obj)
           idx = get(obj.fGData.pumCompare,'Value');
           str = get(obj.fGData.pumCompare,'String');
           str = str{idx};
           if obj.fCurrentFieldHasValidValues && any(strcmp(str,{'is equal to'}))
               % use enumerated values for field
               set(obj.fGData.etValue,'Visible','off');
               idx = get(obj.fGData.pumField,'Value');
               set(obj.fGData.lbValidValues,'String',obj.fDataSetFields(idx).validValues);
               set(obj.fGData.lbValidValues,'Visible','on');
           else
               set(obj.fGData.lbValidValues,'Visible','off');
               set(obj.fGData.etValue,'Visible','on');
           end
        end
        
        % Get data from specialized entry controls
        % d.field
        % d.cmp
        % d.value (char or cellstr)
        function d = getEnteredData(obj,tfConvertValue)
            if nargin < 2
                tfConvertValue = false;
            end
            idx = get(obj.fGData.pumField,'Value');
            d.field = obj.fDataSetFields(idx).name;
            idx = get(obj.fGData.pumCompare,'Value');
            str = get(obj.fGData.pumCompare,'String');
            d.cmp = str{idx};
            if strcmp(get(obj.fGData.etValue,'Visible'),'on')
                d.value = get(obj.fGData.etValue,'String');
            else
                idx = get(obj.fGData.lbValidValues,'Value');
                str = get(obj.fGData.lbValidValues,'String');
                d.value = str(idx);
            end
            
            if tfConvertValue
                if iscellstr(d.value)
                    val = sprintf('%s|',d.value{:});
                    val = val(1:end-1);
                    d.value = val;
                end
            end            
        end
        
        function replaceQuery(obj)
            d = obj.getEnteredData(true);
            obj.jtStopEditing();
            
            row = obj.fJTable.getSelectedRow;
            if row >= 0
                tm = obj.fMTable.getTableModel;
                tm.setValueAt(d.field,row,0);
                tm.setValueAt(d.cmp,row,1);
                tm.setValueAt(d.value,row,2);
            end
        end
            
        function appendQuery(obj)
            d = obj.getEnteredData(true);
            obj.jtStopEditing();
                        
            newRowData = {d.field d.cmp d.value};
            obj.fMTable.getTableModel.addRow(newRowData);
            
            % Move the selection to Column A of this new row
            obj.fJTable.changeSelection(obj.fJTable.getRowCount-1,0,false,false);
            
            if obj.fJTable.getRowCount==1
                set(obj.fGData.pbDelete,'Enable','on');
                set(obj.fGData.pbDeleteAll,'Enable','on');
            end
        end
        
        function insertQuery(obj)
            d = obj.getEnteredData(true);
            obj.jtStopEditing();
            newRowData = {d.field d.cmp d.value};
            obj.fMTable.getTableModel.insertRow(max(0,obj.fJTable.getSelectedRow),newRowData);
            if obj.fJTable.getRowCount==1
                set(obj.fGData.pbDelete,'Enable','on');
                set(obj.fGData.pbDeleteAll,'Enable','on');
            end
        end
        
        function deleteQuery(obj)
            obj.jtStopEditing();
            jt = obj.fJTable;
            
            % If there are any rows displayed, then delete the currently-selected row
            rowCount = jt.getRowCount;
            if rowCount > 0  % might be==0 during slow processing & user double-click
                currentRow = max(0,jt.getSelectedRow);
                currentCol = max(0,jt.getSelectedColumn);
                obj.fMTable.getTableModel.removeRow(currentRow);
                if currentRow > jt.getRowCount
                    currentRow = currentRow-1;
                end
                jt.changeSelection(currentRow, currentCol, false, false);
            end
            if jt.getRowCount <= 0
                set(obj.fGData.pbDelete,'Enable','off');
                set(obj.fGData.pbDeleteAll,'Enable','off');
            end
        end
                
        function deleteAllQueries(obj)
            obj.jtStopEditing();
            obj.fMTable.setNumRows(0);
            set(obj.fGData.pbDelete,'Enable','off');
            set(obj.fGData.pbDeleteAll,'Enable','off');
        end
                
        %%% jtable
        function jtStopEditing(obj)
            component = obj.fJTable.getEditorComponent;
            if ~isempty(component)
                event = javax.swing.event.ChangeEvent(component);
                obj.fJTable.editingStopped(event);
            end
        end
        %%% end jtable
                    
        function setCurrentQuerySet(obj,fullfname)
            obj.fCurrentQuerySet = fullfname;
            [~,n] = fileparts(fullfname);
            set(obj.fGData.pnlOuter,'Title',sprintf('Query: %s',n));
        end
        
        function fullfname = getCurrentQuerySet(obj)
            fullfname = obj.fCurrentQuerySet;            
        end
        
        function fileWrite(obj)
            dat = cell(obj.fMTable.getData);
            dat(:,2) = cellfun(@(x)obj.cmpTxt2Symbol(x),dat(:,2),'UniformOutput',false);            
            fh = fopen(obj.fCurrentQuerySet,'w');
            if fh==-1
                return;
            end
            for c = 1:size(dat,1)
                fprintf(fh,'%s %s %s\n',dat{c,:});
            end
            fclose(fh);
        end
        
        function fileRead(obj,fname)
            fh = fopen(fname);
            if isequal(fh,-1)
                return;
            end
            dat = textscan(fh,'%s%s%s');
            dat = [dat{:}];
            dat(:,2) = cellfun(@(x)obj.cmpSymbol2Txt(x),dat(:,2),'UniformOutput',false);
            obj.deleteAllQueries();
            obj.fMTable.setData(dat);
            drawnow; % sometimes appears needed to get Delete/DeleteAll buttons to refresh
            if obj.fJTable.getRowCount>0
                set(obj.fGData.pbDelete,'Enable','on');
                set(obj.fGData.pbDeleteAll,'Enable','on');
            end
        end        
      
    end

    methods (Static,Hidden)
        
        function qObj = queryObj(fldStr,cmpStr,valStr)
            if isempty(fldStr)||isempty(cmpStr)||isempty(valStr)
                qObj = [];
                return;
            end
            vals = textscan(valStr,'%s','Delimiter','|');
            vals = vals{1};
            if numel(vals)==1
                qObj = SAGE.Query.Compare(fldStr,cmpStr,vals{1});
            else
                if ~strcmp(cmpStr,'=')
                    error('QueryComponentPanel:multipleVals','Multiple values for non-equality comparison.');
                end
                Nqry = numel(vals);
                qs = cell(Nqry,1);
                for c = 1:Nqry
                    qs{c} = SAGE.Query.Compare(fldStr,cmpStr,vals{c});
                end
                qObj = SAGE.Query.Any(qs{:});             
            end            
        end
        
        function y = cmpTxt2Symbol(x)
            if ~isempty(x)
                assert(ischar(x));
                txts = OlyDat.DataSelector.ComparePairs(:,1);
                tf = strcmp(x,txts);
                assert(nnz(tf)==1,'Unknown comparison type.');
                y = OlyDat.DataSelector.ComparePairs{tf,2};
            else
                y = '';
            end
        end
        
        function y = cmpSymbol2Txt(x)
            assert(ischar(x));
            symbols = OlyDat.DataSelector.ComparePairs(:,2);
            tf = strcmp(x,symbols);
            assert(nnz(tf)==1,'Unknown comparison type.');
            y = OlyDat.DataSelector.ComparePairs{tf,1};
        end        
        
        function dsFieldObjs = customizeDataSetFields(dataPuller,dataSetFam,dsFieldObjs)
            customInfo = dataPuller.getCustomFieldPulldownInfo();
            if isempty(customInfo)
                return;
            end
            assert(isstruct(customInfo) && all(isfield(customInfo,{'Name';'DataSetName'})));
            
            for s = customInfo(:)'
                % remove elements of dsFieldObjs with s.Name
                tf = strcmp(s.Name,{dsFieldObjs.name}');
                dsFieldObjs(tf,:) = [];
                
                % add in custom fieldobj
                newObj = dataSetFam.dataSet(s.DataSetName).field(s.Name);
                dsFieldObjs(end+1,1) = newObj; %#ok<AGROW>
            end
        end
        
    end
end