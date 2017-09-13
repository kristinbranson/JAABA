% classdef CurationFile < handle
%    
%     properties (Constant)
%         fileExtension = '.xlsx';
%         xlsSheetName = 'Sheet1';
%         headerText = 'OlyDat Curation File. Do not modify the first two rows of this file.';
%         col2Fieldname = {'experiment_id';'experiment_name';'exp_datetime';'manual_pf';'flag_redo';'notes'};
%         fieldTypes = {'numeric';'char';'char';'char';'char';'char'};
%     end
%     
%     methods (Static)       
%         
%         % eid2CurationInfo: map from expID->curation info struct
%         % eid2XlsRow: map from expID->row in xls file
%         % Nrows: total number of rows in xls file (including header)
%         function [eid2CurationInfo eid2XlsRow Nrows] = read(fname)
%             assert(ischar(fname));
%             fname = OlyDat.CurationFile.addFileExtensionIfNeeded(fname);            
%             if exist(fname,'file')~=2
%                 error('OlyDat:CurationFile:fileNotFound',...
%                     'File ''%s'' not found.',fname);
%             end
%             
%             try
%                 [~,~,raw] = xlsread(fname,OlyDat.CurationFile.xlsSheetName);
%             catch ME
%                 throw(ME); % TODO
%             end
%             
%             if ~isequal(OlyDat.CurationFile.headerText,raw{1,1}) || ...
%                ~isequal(OlyDat.CurationFile.col2Fieldname,raw(2,1:numel(OlyDat.CurationFile.col2Fieldname))')     
%                 error('OlyDat:CurationFile:corruptCurationFile',...
%                     'Curation file is corrupted.');
%             end
%             
%             eid2CurationInfo = containers.Map('KeyType','double','ValueType','any');
%             eid2XlsRow = containers.Map('KeyType','double','ValueType','double');
%             
%             Nrows = size(raw,1);
%             for c = 3:Nrows % first two rows are header
%                 eid = raw{c,1}; % hardcoded first col is expID
%                 if isnan(eid) || isempty(eid) || ~isnumeric(eid)
%                     warning('OlyDat:CurationFile:corruptRow',...
%                         'Skipping corrupt row in curation file.');
%                     continue;
%                 end
%                 if eid2CurationInfo.isKey(eid)
%                     warning('OlyDat:CurationFile:corruptRow',...
%                         'Repeated experiment ID in curation file.');
%                     continue;
%                 end
%                 
%                 info = struct();
%                 tfValidRow = true;
%                 for d = 1:numel(OlyDat.CurationFile.col2Fieldname)
%                     value = raw{c,d};
%                     fld = OlyDat.CurationFile.col2Fieldname{d};
%                     
%                     % deal with empty cells (imported as nan)
%                     if isnan(value) 
%                         switch OlyDat.CurationFile.fieldTypes{d}
%                             case 'numeric'
%                                 % none; use nan
%                             case 'char'
%                                 value = '';
%                         end
%                     end
%                     
%                     % custom field checking/conversion
%                     switch fld
%                         case 'manual_pf'
%                             if ~any(strcmp(value,{'U';'P';'F'}))
%                                 tfValidRow = false;
%                             end
%                         case 'flag_redo'
%                             if isnumeric(value)
%                                 value = num2str(value);
%                             end
%                             if ~any(strcmp(value,{'0';'1'}))
%                                 tfValidRow = false;
%                             end
%                         otherwise
%                             % none
%                     end
%                     if ~tfValidRow
%                         warning('OlyDat:CurationFile:corruptRow',...
%                             'Invalid value of ''%s'' found in curation file, row %d. Skipping row.',fld,c);
%                         break;
%                     end
%                     
%                     info.(fld) = value;
%                 end
%                 
%                 if tfValidRow
%                     eid2CurationInfo(eid) = info;
%                     eid2XlsRow(eid) = c;
%                 end
%             end
%             
%         end
%         
%         function append(fname,curationInfo)
%             assert(ischar(fname));
%             fname = OlyDat.CurationFile.addFileExtensionIfNeeded(fname);
%             
%             assert(isa(curationInfo,'containers.Map'));
% 
%             % create a new curation file if necessary
%             if exist(fname,'file')~=2
%                 OlyDat.CurationFile.createNewCurationFile(fname);
%             end
%             
%             % read the curation file
%             [eid2Info eid2XlsRow NXlsRows] = OlyDat.CurationFile.read(fname);
%             
%             % Unfortunately, xlswrite doesn't appear to allow making
%             % noncontiguous edits to a file. Modifications to the curation
%             % file therefore cannot be made in a single operation.
%             %
%             % In the following try-catch block we attempt to write/append
%             % all new curation info to the xls file.
%             
%             % backup existing file
%             backupfname = tempname;
%             tf = copyfile(fname,backupfname);
%             if ~tf
%                 error('OlyDat:CurationFile','Could not create backup file ''%s''.',backupfname);
%             end
%             
%             tfFailedDuringXlsWrite = false;
%             try
%                 
%                 hWB = waitbar(0,'Please wait...','Name','Writing Curation Info');
%                 
%                 % loop over new info
%                 newKeys = curationInfo.keys;
%                 Nfields = numel(OlyDat.CurationFile.col2Fieldname);
%                 oneRowPastEnd = NXlsRows+1;
%                 for c = 1:numel(newKeys)
%                     waitbar(c/numel(newKeys),hWB);
%                     
%                     k = newKeys{c};
%                     
%                     % prepare a row for writing into XLS file.
%                     newInfo = curationInfo(k);
%                     assert(all(isfield(newInfo,OlyDat.CurationFile.col2Fieldname)));
%                     xlsData = cell(1,Nfields);
%                     for d = 1:Nfields
%                         xlsData{d} = newInfo.(OlyDat.CurationFile.col2Fieldname{d});
%                     end
% 
%                     assert(eid2Info.isKey(k)==eid2XlsRow.isKey(k));
%                     if eid2XlsRow.isKey(k)
%                         % key already exists in curation file.
%                         xlsRow = eid2XlsRow(k);
%                     else
%                         % key doesnt exist; append to end of file.
%                         xlsRow = oneRowPastEnd;
%                         oneRowPastEnd = oneRowPastEnd+1;
%                     end
%                     xlsRange = sprintf('A%d',xlsRow);
%                     xlswrite(fname,xlsData,OlyDat.CurationFile.xlsSheetName,xlsRange);
%                 end
%             catch ME
%                 tfFailedDuringXlsWrite = true;
%             end
%             
%             delete(hWB);
%             
%             if tfFailedDuringXlsWrite
%                 try
%                     % may fail, eg if fname is in use by excel
%                     movefile(backupfname,fname); 
%                 catch
%                     delete(backupfname);
%                 end
%                 error('Error while writing new curation file: %s',ME.message);
%             else
%                 delete(backupfname); % backup file should exist
%             end
%         end
%         
%         % Merge curation infos from src into dest. The EIDs in dest are
%         % unchanged; if an eid is in both dest and src, the info from src
%         % is written to dest.
%         function merge(srcInfo,destInfo)
%             assert(isa(srcInfo,'containers.Map'));
%             assert(isa(destInfo,'containers.Map'));            
%             srcEids = srcInfo.keys;
%             for c = 1:numel(srcEids)
%                 k = srcEids{c};
%                 if destInfo.isKey(k)
%                     destInfo(k) = srcInfo(k);
%                 end
%             end
%         end
%         
%     end
%     
%     methods (Static,Hidden)
%         
%         function createNewCurationFile(fname)
%             assert(ischar(fname));
%             fname = OlyDat.CurationFile.addFileExtensionIfNeeded(fname);
%             if exist(fname,'file')==2
%                 error('OlyDat:CurationFile:fileExists',...
%                     'Cannot create curation file ''%s''. File already exists.',fname);
%             end
%             
%             data(2,:) = OlyDat.CurationFile.col2Fieldname';
%             data{1,1} = OlyDat.CurationFile.headerText;
%             try
%                 xlswrite(fname,data,OlyDat.CurationFile.xlsSheetName);
%             catch ME
%                 throw ME;
%             end
%         end
%         
%         function fname = addFileExtensionIfNeeded(fname)
%             [p f e] = fileparts(fname);
%             if isempty(e)
%                 fname = fullfile(p,[f OlyDat.CurationFile.fileExtension]);
%             end
%         end
%         
%     end
%     
%     
% end