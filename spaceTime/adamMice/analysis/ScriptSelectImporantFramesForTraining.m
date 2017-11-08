addpath /groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/misc;
addpath /groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/filehandling;
rootdir = '/tier2/hantman';

% datafiles = {'m119nocnoraw.csv'
%   'm119cnoraw.csv'};
% datafiles = {'m119forallen.mat'};
% datafiles = {'m122raw.csv'
%   'm130raw.csv'};
datafiles = {'data/PostProcessed_M122_20141210.csv'
'data/PostProcessed_M130_20141210.csv'};

rawdata = [];
for i = 1:numel(datafiles),
  
  [~,~,ext] = fileparts(datafiles{i});
  switch ext,
    case '.csv'
      rawdatacurr = ReadRawDataFile(datafiles{i});
      rawdata = structappend(rawdata,rawdatacurr);
    case '.mat'
      rawdatacurr = load(datafiles{i});
      rawdata = structappend(rawdata,rawdatacurr.data);
  end
  
end

%% get behavior names and order

% fnsspecial = {'auto_Grab_success','auto_Grab_successtype'};
% fns = fieldnames(rawdata);
% datafns = setdiff(fns(~cellfun(@isempty,regexp(fns,'^auto'))),fnsspecial);
% nstats = numel(datafns);
% 
% % some behaviors don't have counts
% tmp = regexp(datafns,'^auto_([^_]+)_0$','tokens','once');
% behaviors = unique([tmp{:}]);

behaviors = {'Lift','Handopen','Grab','Sup','Atmouth','Chew'};
behaviorfns = {'auto_GS00_Lift_0','auto_GS00_Handopen_0','auto_GS00_Grab_0','auto_GSSS_Sup_0',...
  'auto_GSSS_Atmouth_0','auto_GSSS_Chew_0'};
nbehaviors = numel(behaviors);

%% save frames of each

importantframes = [];
for i = 1:numel(rawdata),
  expdir = strrep(rawdata(i).expfull,'\','/');
  expdir = regexprep(expdir,'^[A-Z]:',rootdir);
  scurr = struct('expdir',expdir,...
    'firstbehavior',rawdata(i).auto_Lift_0,...
    'firstsuccessbehavior',rawdata(i).auto_GS00_Lift_0,...
    'lastsuccessbehavior',nan,...
    'lastbehavior',nan);
  tssuccess = [rawdata(i).auto_GS00_Lift_0,rawdata(i).auto_GS00_Handopen_0,rawdata(i).auto_GS00_Grab_0,rawdata(i).auto_GSSS_Sup_0,rawdata(i).auto_GSSS_Atmouth_0,rawdata(i).auto_GSSS_Chew_0];
  scurr.lastsuccessbehavior = max(tssuccess);
  ts = [rawdata(i).auto_Lift_0,rawdata(i).auto_Lift_1,rawdata(i).auto_GS00_Lift_0,rawdata(i).auto_GSSS_Lift_0,...
    rawdata(i).auto_Handopen_0,rawdata(i).auto_GS00_Handopen_0,rawdata(i).auto_GSSS_Handopen_0,...
    rawdata(i).auto_Grab_0,rawdata(i).auto_GS00_Grab_0,rawdata(i).auto_GSSS_Grab_0,...
    rawdata(i).auto_Sup_0,rawdata(i).auto_GSSS_Sup_0,...
    rawdata(i).auto_Atmouth_0,rawdata(i).auto_GSSS_Atmouth_0,...
    rawdata(i).auto_Chew_0,rawdata(i).auto_GSSS_Chew_0];    
  scurr.lastbehavior = max(ts);
  importantframes = structappend(importantframes,scurr);  
end

save M122M130ImportantFrames20141210.mat importantframes