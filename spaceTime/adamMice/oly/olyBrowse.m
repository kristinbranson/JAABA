function browser = olyBrowse(data)
%startBrowser Start OlyDat Browser for Adam mice exps.

if isempty(data)
  warning('olyBrowse:noData','No data.');
  return;
end

data = ExpPP.convertCustomGroupsToStrs(data);
%data = ExpPP.convertLogicalFieldToStrs(data,'auto_Grab_success','failure','success');

digits = ceil(log10(numel(data)+1)); % number of digits required to number experiments 

for i = 1:numel(data)
  mouseval = str2double(data(i).mouse);  
  if isnan(mouseval)
    mouseval = 0;
  end
  data(i).experiment_id = mouseval*10^digits + i;  
  data(i).experiment_name = data(i).id;
  data(i).line_name = data(i).mouse;
  data(i).mousecondition = [data(i).mouse '#' data(i).condition];
  
  % AL 20150129 semi-hack, add some indicator vars for successrate anls
  data(i).auto_Grab_success1 = double(strcmp(data(i).auto_Grab_successtype,'successful_1grab'));
  data(i).auto_Grab_success2plus = double(strcmp(data(i).auto_Grab_successtype,'successful_2plusgrab'));
end

% Determine stats to include in 'statistic'
btn = questdlg('Select subset of statistics to include:','Select Statistics',...
  'Default (limited)','All','Default (limited)');
if isempty(btn)
  return;
end
switch btn
  case 'Default (limited)'
    fldsStat = {
      'auto_Lift_0'
      'auto_Handopen_0'
      'auto_Grab_0'
      'auto_Sup_0'
      'auto_Atmouth_0'
      'auto_Chew_0'
      'auto_Lift_1'
      'auto_Grab_num'
      'auto_Grab_success'
      'auto_Grab_success1'
      'auto_Grab_success2plus'
      'auto_GS00_Lift_0'
      'auto_GS00_Handopen_0'
      'auto_GS00_Grab_0'
      'auto_del_GS00_Handopen_Lift'
      'auto_del_GS00_Grab_Lift'
      'auto_del_GS00_Grab_Handopen'
      'auto_del_GS01_Grab_GS00_Grab'
      'auto_del_GSSS_Grab_Handopen'
      'auto_del_Lift1_Lift0'
      'auto_del_GSSS_Sup_Grab'
      'auto_del_GSSS_Atmouth_Grab'
      'auto_del_GSSS_Atmouth_Sup'      
      };
  case 'All'
    fldsAll = fieldnames(data);
    tf01 = regexpmatch(fldsAll,'^auto_.+_[01]$');
    tfbl = regexpmatch(fldsAll,'^auto_.+_bl$');
    tfnum = regexpmatch(fldsAll,'^auto_.+_num$');
    tfdel = regexpmatch(fldsAll,'^auto_del');
    tfGS0 = regexpmatch(fldsAll,'_GS0[01]_');
    tfGSSS = regexpmatch(fldsAll,'_GSSS_');
    tfmisc = ismember(fldsAll,{'auto_numGS' 'auto_Grab_success' 'auto_Grab_success1' 'auto_Grab_success2plus'});
    tf = tf01 | tfbl | tfnum | tfdel | tfGS0 | tfGSSS | tfmisc;
    fldsStat = fldsAll(tf);
  otherwise
    assert(false);
end
% Add all corresponding labl stats
fldsStatLabl = regexprep(fldsStat,'^auto_','labl_');
fldsStat = [fldsStat(:);fldsStatLabl(:)];

tf = ismember(fldsStat,fieldnames(data));
if any(~tf)
  warning('olyBrowse:missingFields','Fields not present in data: %s',...
    civilizedStringFromCellArrayOfStrings(fldsStat(~tf)));
end
fldsStat = fldsStat(tf);
% Behavioral quantities
for i = 1:numel(fldsStat)
  bst(i,1) = OlyDat.BrowserStat(fldsStat{i},fldsStat{i}); %#ok<AGROW>
end

% Grouping quantities
uniqueDates = unique({data.date}');
gst(1) =       OlyDat.BrowserStat('','<no group>');
gst(end+1,1) = OlyDat.BrowserStat('date','date',uniqueDates);
allSucc = double([data.auto_Grab_success]);
gst(end+1,1) = OlyDat.BrowserStat('auto_Grab_success','auto_Grab_success',unique(allSucc));
gst(end+1,1) = OlyDat.BrowserStat('auto_Grab_successtype','auto_Grab_successtype',{'successful_1grab' 'successful_2plusgrab' 'unsuccessful'});
allSucc1 = double([data.auto_Grab_success1]);
gst(end+1,1) = OlyDat.BrowserStat('auto_Grab_success1','auto_Grab_success1',unique(allSucc1));
allSucc2plus = double([data.auto_Grab_success2plus]);
gst(end+1,1) = OlyDat.BrowserStat('auto_Grab_success2plus','auto_Grab_success2plus',unique(allSucc2plus));
allmice = {data.mouse}';
gst(end+1,1) = OlyDat.BrowserStat('mouse','mouse',unique(allmice));
[data,gstNew] = lclAddGroupingVars(data,gst(2:end)); % Don't include <no group> group when forming pairwise etc
gst = [gst(1);gstNew];

% Auxiliary quantities
successtypesubsets = ComparisonPlot.SUCCESSTYPE_SUBSETS;
for i = 1:numel(successtypesubsets)
  stype = successtypesubsets{i};
  ast(i,1) = OlyDat.BrowserStat(stype,stype); %#ok<AGROW>
end

% Experiment Detail Quantities
est(1)       = OlyDat.BrowserStat('exp');
est(end+1,1) = OlyDat.BrowserStat('date');
est(end+1,1) = OlyDat.BrowserStat('auto_Grab_successtype');

% Plots
hObj = HandleObj; % shared state for Repeat plots cache
plots = {ComparisonPlot; BoxPlotPlot; MultiStatCompare; RepeatBoxPlot(hObj); RepeatMultiStatCompare(hObj); RepeatSuccess(hObj)};

% DetailHandler
expdetailhandler = DetailHandler;

curFlds = OlyDat.CurationField.empty(0,1);
curFlds(end+1,1) = OlyDat.CurationField('experiment_id',false,'Exp. ID',50,'');
curFlds(end+1,1) = OlyDat.CurationField('line_name',false,'',[],'line');
curFlds(end+1,1) = OlyDat.CurationField('experiment_name',false,'Exp. Name',280,'experiment');
curFlds(end+1,1) = OlyDat.CurationField('manual_pf',true,'Manual P/F',70,[],{'P';'F';'U'});
curFlds(end+1,1) = OlyDat.CurationField('manual_curator',true,'Manual Curator',[],[]);
curFlds(end+1,1) = OlyDat.CurationField('manual_curation_date',true,'Curation Date',125,[]);
curFlds(end+1,1) = OlyDat.CurationField('notes_curation',true,'Curation Notes',240);
curFlds(end+1,1) = OlyDat.CurationField('flag_redo',true,'Flag Redo',65,[],{'';'0';'1'});
curFlds(end+1,1) = OlyDat.CurationField('flag_review',true,'Flag Review',70,[],{'';'0';'1'});

curFldInfo.objs = curFlds;
curFldInfo.tblOrder = {'experiment_id';'experiment_name';'manual_pf';...
    'notes_curation';'flag_redo';'flag_review';'manual_curation_date';...
    'manual_curator'};
curFldInfo.tsvOrder = {'line_name';'experiment_name';'manual_pf';...
    'manual_curator';'manual_curation_date';'notes_curation';'flag_redo';'flag_review'};

curationInfoEntryGUIName = [];

browser = OlyDat.Browser('box',data,plots,bst,gst,ast,est,expdetailhandler,...
    curFldInfo,curationInfoEntryGUIName);

function [data,gst] = lclAddGroupingVars(data,gst)
% Augment gst, vector of OlyDat.BrowserStats for grouping variables
%
% - Add custom groups to gst
% - Form all pairwise grouping variables and add to data and gst

[~,customGrps] = ExpPP.ensureCustomGroupInit(data);
customGrpFlds = cellfun(@(x)sprintf(ExpPP.CUSTOMGROUPNAME_PAT,x),customGrps,'uni',0);
Ncgrp = numel(customGrps);

% Add all custom groups
for i = 1:Ncgrp
  fld = customGrpFlds{i};
  vals = {data.(fld)};
  assert(iscellstr(vals)); % We converted to char vals earlier
  unVals = unique(vals);
  assert(numel(unVals)<=2,'Custom group field takes on more than 2 values');
  gst(end+1,1) = OlyDat.BrowserStat(fld,fld,unVals); %#ok<AGROW>
end

% Check all grouping fields; remove unsuitable fields from gst
for i = numel(gst):-1:1
  gF = gst(i).Name;
  if ~isfield(data,gF)
    warning('olyBrowse:missingField','Data is missing field ''%s''.',gF);
    gst(i,:) = [];
    continue;
  end
  
  [tf,msg] = ExpPP.checkCategoricalVar(data,gF);
  if ~tf
    warning('olyBrowse:invalidGroupingField','Skipping invalid grouping field ''%s'': %s',...
      gF,msg);
    gst(i,:) = [];
    continue;
  end
  
  assert(~isempty(gst(i).ValidValues));
end
  
% Form all pairwise groups
[data,gPairs] = ExpPP.addPairwiseGroups(data,{gst.Name}');
for iPairs = 1:numel(gPairs)
  gF = gPairs{iPairs};
  [tf,~,vPairs] = ExpPP.checkCategoricalVar(data,gF);
  assert(tf);
  gst(end+1,1) = OlyDat.BrowserStat(gF,gF,unique(vPairs)); %#ok<AGROW>
end
