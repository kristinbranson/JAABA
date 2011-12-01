function outputLabels(segstarts,segends,labels,behaviors,idx,moviename,matname,trxname,savename,t0,t1,flies)

if ~exist('flies','var'),
  flies = idx;
end

version = 0.1;
LABELIDXOFFSET = -1 + 2; % -1 for 0-indexing, +2 for unknown and bad behaviors. 

fid = fopen(savename,'wb');
if fid < 0,
  error('Could not open file %s for writing',savename);
end

if ~exist('t0','var')
  t0 = 1;
end
if ~exist('t1','var'),
  t1 = max(segends);
end

% format:
% version (double)
% moviename_length (int)
% moviename (char*)
% matname_length (int)
% matname (char*)
% trxname_length (int)
% trxname (char*)
% nflies (int)
% fly ids (int*)
% firstframe (int)
% lastframe (int)
% nbehaviors (int)
% behaviorname1_length (int)
% behaviorname1 (char*)
% ...
% behaviornamen_length (int)
% behaviornamen (char*)
% nbouts (int)
% startframe_1 (int)
% endframe_1 (int)
% behavioridx_1 (int)
% ...
% startframe_nbouts (int)
% endframe_nbouts (int)
% behavioridx_nbouts (int)

% header
fwrite(fid,version,'double');
fwrite(fid,length(moviename),'int');
fwrite(fid,moviename,'char');
fwrite(fid,length(matname),'int');
fwrite(fid,matname,'char');
fwrite(fid,length(trxname),'int');
fwrite(fid,trxname,'char');
fwrite(fid,length(flies),'int');
fwrite(fid,flies,'int');
fwrite(fid,t0-1,'int');
fwrite(fid,t1-1,'int');
fwrite(fid,length(behaviors),'int');

% behaviors
for i = 1:length(behaviors),
  fn = behaviors{i};
  fwrite(fid,length(fn),'int');
  fwrite(fid,fn,'char');
end

% sort by start frame
[segstarts,order] = sort(segstarts);
segends = segends(order);
labels = labels(order);
[keep,labelidx] = ismember(labels,behaviors);

numDiscarded = sum(~keep);
if numDiscarded > 0 % CSC 20110322: warn if annotations are discarded
    warning('Warning: This should never happen: %d annotations use an unsupported labelname and are discarded\n');
end

% number of bouts
fwrite(fid,nnz(keep),'int');

% write each bout
data = [segstarts(keep)-1;segends(keep);labelidx(keep)+LABELIDXOFFSET];
% data = [segstarts(keep)-1;segends(keep)-1;labelidx(keep)+LABELIDXOFFSET]; % CSC 20110321: Fixed Typo?
fwrite(fid,data(:),'int');
  
fclose(fid);
