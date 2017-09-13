classdef CurationFileWriter < handle
%CurationFileWriter Read/Write interface to OlyDat Curation TSV files
   
    properties (Constant)
        fileExtension = '.tsv';
    end
    
    properties (Hidden)
        fFields; %NFieldsx1 array of OlyDat.CurationField
        fExpNameColIdx; % scalar int; idx of 'experiment' field in TSV file
    end
    
    properties (Dependent)
        numFields;
        fieldnames;
    end
    
    methods        
        function v = get.numFields(obj)
            v = numel(obj.fFields);
        end
        function v = get.fieldnames(obj)
            v = {obj.fFields.Name}';
        end
    end
    
    methods
        
        function obj = CurationFileWriter(flds)
            assert(isa(flds,'OlyDat.CurationField'));            
            obj.fFields = flds;
            fnames = {flds.TSVName}';
            tfExpName = strcmp('experiment',fnames);
            assert(nnz(tfExpName)==1,'Fields must include a field ''experiment''.');
            obj.fExpNameColIdx = find(tfExpName);
        end
                
        % curationInfo: map from expName->curation info struct
        % nRows: total number of rows in TSV file (including header)
        function [curationInfo Nrows] = read(obj,fname)
            assert(ischar(fname));
            fname = OlyDat.CurationFileWriter.addFileExtensionIfNeeded(fname);
            if exist(fname,'file')~=2
                error('OlyDat:CurationFileWriter:fileNotFound',...
                    'File ''%s'' not found.',fname);
            end
            fh = fopen(fname,'r');
            if isequal(fh,-1)
                error('OlyDat:CurationFileWriter:cantOpenFile',...
                    'Error opening file ''%s''.',fname);
            end
            
            % verify header/fields
            row = OlyDat.CurationFileWriter.readRow(fh);
            if isempty(row) || ~isequal(row(:),{obj.fFields.TSVName}')
                fclose(fh);
                error('OlyDat:CurationFileWriter:corruptCurationFile',...
                    'Curation file is corrupted; field names missing/incorrect.');
            end
            
            curationInfo = containers.Map('KeyType','char','ValueType','any');
            
            Nrows = 1;
            while 1
                row = OlyDat.CurationFileWriter.readRow(fh);
                if isempty(row)
                    break;
                end
                Nrows = Nrows+1;

                if numel(row)~=obj.numFields
                    warnNoBacktrace('OlyDat:CurationFileWriter:corruptRow',...
                        'Skipping corrupt row number %d in curation file.',Nrows);
                    continue;
                end
                
                expName = row{obj.fExpNameColIdx};
                if isempty(expName)
                    warnNoBacktrace('OlyDat:CurationFileWriter:corruptRow',...
                        'Skipping corrupt row number %d in curation file.',Nrows);
                    continue;
                end
                
                if curationInfo.isKey(expName)
                    warnNoBacktrace('OlyDat:CurationFileWriter:corruptRow',...
                        'Repeated experiment name ''%s'' in curation file at row number %d. Skipping row.',expName,Nrows);
                    continue;
                end
                
                info = struct();
                tfValidRow = true;
                for d = 1:obj.numFields
                    value = strtrim(row{d});
                    
                    % deal with empty vals (also whitespace vals,
                    % strtrimmed to empty)
                    if isempty(value)
                        value = '';
                    end
                    
                    % verify enumerated fields         
                    fld = obj.fFields(d);
                    tfValidValue = fld.validateValue(value);
                    if ~tfValidValue
                        tfValidRow = false;
                        warnNoBacktrace('OlyDat:CurationFileWriter:corruptRow',...
                            'Invalid value for field ''%s'' found in curation file, row %d. Skipping row.',...
                            fld.TSVName,Nrows);
                        break;
                    end
                    
                    info.(fld.Name) = value;
                end
                
                if tfValidRow
                    curationInfo(expName) = info;
                end
                
            end
            
            fclose(fh);            
        end
                
        % curationInfo: map from expName->curInfo
        % rowsWritten.append: number of rows appended (note that appended
        % rows don't necessarily show up at the "bottom" of the file due to
        % sorting by expName)
        % rowsWritten.overwrite: number of rows overwritten
        function [tfSuccess rowsWritten] = append(obj,fname,curationInfo,tfWarnDlgIfDangerous)
            if nargin < 4
                tfWarnDlgIfDangerous = true;
            end
            
            tfSuccess = false;
            rowsWritten = struct();
                
            assert(ischar(fname));
            fname = OlyDat.CurationFileWriter.addFileExtensionIfNeeded(fname);
            assert(exist(fname,'file')==2);
            assert(isa(curationInfo,'containers.Map'));
            assert(isscalar(tfWarnDlgIfDangerous)&&islogical(tfWarnDlgIfDangerous));
            
            % read the curation file
            [existingInfo Nrows] = obj.read(fname);
            
            % warn if there are one or more corrupt rows; these will not
            % get written in the new file.
            if size(existingInfo,1)+1 ~= Nrows
                if tfWarnDlgIfDangerous
                    queststr = sprintf('Curation file ''%s'' has corrupt rows. These rows will be removed from the file if you continue. Continue anyway?',fname);
                    btn = questdlg(queststr,'Corrupt Rows Found','Continue anyway', 'No, cancel', 'No, cancel');
                    switch btn
                        case {'' 'No, cancel'}
                            return;
                        case 'Continue anyway'
                            % none
                    end
                else
                    warnNoBacktrace('OlyDat:CurationFileWriter:corruptRows',...
                        'Curation file ''%s'' has corrupt rows. These rows will be removed from the file.',fname);
                end
            end
            
            if tfWarnDlgIfDangerous
                expNamesInFile = existingInfo.keys;
                newExpNames = curationInfo.keys;
                tf = ismember(newExpNames(:),expNamesInFile(:));
                if any(tf)
                    expNameStr = sprintf('%s, ',newExpNames{tf});
                    expNameStr = expNameStr(1:end-2);
                    warnStr = sprintf('Continuing will overwrite existing data in curation file for experiment(s) %s. Is this okay?',expNameStr);
                    btn = questdlg(warnStr,'Overwriting data','Yes, continue','No, cancel','Yes, continue');
                    switch btn
                        case {'' 'No, cancel'}
                            return;
                        case 'Yes, continue'
                            % none
                    end
                end
            end
            
            % merge in new data
            [rowsWritten.append rowsWritten.overwrite] = ...
                OlyDat.CurationFileWriter.merge(curationInfo,existingInfo,true);

            % backup file and open file for writing
            backupfname = tempname;
            tf = copyfile(fname,backupfname);
            if ~tf
                error('OlyDat:CurationFileWriter','Could not create backup file ''%s''.',backupfname);
            end
            assert(exist(backupfname,'file')==2);
            
            fh = fopen(fname,'w');
            if isequal(fh,-1)
                delete(backupfname);
                error('OlyDat:CurationFileWriter:cantOpenFile',...
                    'Error opening file ''%s'' for writing.',fname);
            end
            
            tfRestoreBackup = false;
            try
                obj.writeCurationFileHeader(fh);
                                
                % for now, always save with expNames sorted
                expNames = sort(existingInfo.keys);
                for c = 1:numel(expNames)
                    info = existingInfo(expNames{c});
                    row = obj.infoStructToRow(info);
                    zlclWriteTabDelimitedRow(fh,row);
                end
            catch ME
                tfRestoreBackup = true;
                errMsg = sprintf('Error while writing to curation file: %s Restoring backup TSV file.\n',ME.message);
            end
            fclose(fh);
            
            if tfRestoreBackup
                warndlg(errMsg,'Save Error');
                [tf msg] = movefile(backupfname,fname); 
                if ~tf
                    % unsuccessful restore of backup
                    warndlg(sprintf('Error restoring backup file: %s. Open backup file ''%s'' directly to access your original data.\n',msg,backupfname),...
                        'Restore Error');
                end
            else
                % Only codepaths that end here were successful appends
                tfSuccess = true;
                delete(backupfname);
            end
        end
        
        % This discards any existing contents if fname already exists.
        % info: map from expName->curation info struct       
        function [tfSuccess rowsWritten] = writeNewFile(obj,fname,info)
            tfSuccess = false; %#ok<NASGU>
            
            assert(ischar(fname));
            fname = OlyDat.CurationFileWriter.addFileExtensionIfNeeded(fname);
            assert(isa(info,'containers.Map'));
            
            fh = fopen(fname,'w');
            if isequal(fh,-1)
                error('OlyDat:CurationFileWriter:cantOpenFile',...
                    'Error opening file ''%s'' for writing.',fname);
            end

            tfErr = false;
            try
                obj.writeCurationFileHeader(fh);
                
                % for now, always save with expNames sorted
                expNames = sort(info.keys);
                for c = 1:numel(expNames)
                    s = info(expNames{c});
                    row = obj.infoStructToRow(s);
                    zlclWriteTabDelimitedRow(fh,row);
                end
            catch ME
                tfErr = true;
            end
            
            fclose(fh);
            if tfErr
                errMsg = sprintf('Error while writing new curation file: %s\n',ME.message);
                warnDlg(errMsg,'Error');
                error('OlyDat:CurationFileWriter','Error writing new TSV file.');
            end
                        
            % At the moment, this method is either successful, or it
            % throws. There is no codepath that returns tfSuccess==false.
            tfSuccess = true; 
            rowsWritten = info.Count;
        end
        
    end
    
    methods (Hidden)
        
        function writeCurationFileHeader(obj,fh)
            header = {obj.fFields.TSVName}';
            zlclWriteTabDelimitedRow(fh,header);
        end
        
        % This discards any existing contents if fname already exists.
        function createNewCurationFile(obj,fname)
            assert(ischar(fname));
            fname = OlyDat.CurationFileWriter.addFileExtensionIfNeeded(fname);
            fh = fopen(fname,'w');
            if fh==-1
                error('OlyDat:CurationFileWriter:cantOpenFile','Error opening file ''%s''.',fname);
            end
            obj.writeCurationFileHeader(fh);
            fclose(fh);
        end
        
        function row = infoStructToRow(obj,info)
            row = cell(1,obj.numFields);
            for c = 1:obj.numFields
                fld = obj.fFields(c);
                if ~isfield(info,fld.Name)
                    warnNoBacktrace('OlyDat:CurationFileWriter:missingField',...
                        'Field ''%s'' missing from curation info structure.',fld.Name);
                    val = ' ';
                else
                    val = info.(fld.Name);
                end
                
                if ~fld.validateValue(val)
                    warnNoBacktrace('OlyDat:CurationFileWriter:invalidField',...
                        'Invalid value for field ''%s''.',fld.Name);
                end
                
                row{c} = val;
            end
        end
        
    end
    
    methods (Static,Hidden)
        
        % Create a new map from an existing map by "changing keys".
        % map: map from oldKey -> val
        % oldKey2newKey: map from oldKey -> newKey. This map must be 1-to-1.
        % newMap: map from newKey -> val
        function newMap = changeKey(map,oldKey2newKey)
            newMap = containers.Map('KeyType',oldKey2newKey.ValueType,'ValueType',map.ValueType);            
            oldKeys = map.keys;
            for c = 1:numel(oldKeys)
                oldk = oldKeys{c};
                val = map(oldk);
                assert(oldKey2newKey.isKey(oldk));
                newk = oldKey2newKey(oldk);
                assert(~newMap.isKey(newk));
                newMap(newk) = val;
            end
        end
        
        % Merge curation infos from src into dest.
        % If tfAppend is false, the EIDs in dest are unchanged; info for an
        % eid is only copied from src to dest if the eid exists in both
        % maps.
        % If tfAppend is true, every EID in src is copied to dest, possibly
        % enlarging dest.
        function [Nappend Noverwrite] = merge(srcInfo,destInfo,tfAppend)
            assert(isa(srcInfo,'containers.Map'));
            assert(isa(destInfo,'containers.Map'));            
            srcEids = srcInfo.keys;
            Nappend = 0;
            Noverwrite = 0;
            for c = 1:numel(srcEids)
                k = srcEids{c};
                isKeyDest = destInfo.isKey(k);
                if tfAppend || isKeyDest
                    destInfo(k) = srcInfo(k);
                    if isKeyDest
                        Noverwrite = Noverwrite+1;
                    else
                        Nappend = Nappend+1;
                    end
                end
            end
        end
                     
        function fname = addFileExtensionIfNeeded(fname)
            [p f e] = fileparts(fname);
            if isempty(e)
                fname = fullfile(p,[f OlyDat.CurationFileWriter.fileExtension]);
            end
        end
        
        function dat = readRow(fh)
            t = fgetl(fh);
            if isequal(t,-1)
                dat = [];
                return;
            end
            
            % Workaround annoying textscan/regexp behavior. If you have
            %   str = sprintf('X\tY\t\tZ');
            %   a = textscan(str,'%s','Delimiter','\t');
            % Then a{1} = {X';'Y';'';'Z'}. However if you have
            %   str = sprintf('X\tY\tZ\t');
            %   a = ...
            % Then a{1} = {'X';'Y';'Z'}.
            %
            % textscan() does not recognize the trailing empty field.
            tfTextscanWorkaround = false;
            if isequal(t(end),sprintf('\t'))
                t(end+1) = '#';
                tfTextscanWorkaround = true;
            end
                                
            dat = textscan(t,'%s','Delimiter','\t');
            dat = dat{1};
            if tfTextscanWorkaround
                assert(isequal(dat{end},'#'));
                dat{end} = '';
            end
        end
        
    end
    
end

function zlclWriteTabDelimitedRow(fh,row)
rowStr = sprintf('%s\t',row{:});
rowStr(end) = sprintf('\n'); % replace final tab with newline
fprintf(fh,rowStr);
end
