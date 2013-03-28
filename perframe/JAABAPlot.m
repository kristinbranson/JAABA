function varargout = JAABAPlot(varargin)
%JAABAPLOT: Plot the output of JAABA
%
% This program is part of JAABA.
%
% JAABA: The Janelia Automatic Animal Behavior Annotator
% Copyright 2012, Kristin Branson, HHMI Janelia Farm Resarch Campus
% http://jaaba.sourceforge.net/
% bransonk@janelia.hhmi.org
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License (version 3 pasted in LICENSE.txt) for 
% more details.

% Last Modified by GUIDE v2.5 11-Mar-2013 15:51:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @JAABAPlot_OpeningFcn, ...
                   'gui_OutputFcn',  @JAABAPlot_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
   gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% ---
function handles=initialize(handles)

handles.grouplist={};
handles.groupvalue=1;
handles.experimentlist={{}};
handles.experimentvalue={1};
handles.classifierlist={};
handles.classifiervalue=1;
handles.configurations={};
handles.scorefiles={};
handles.analysis='';
handles.behaviorlist={};
handles.behaviornot=0;
handles.behaviorvalue=1;
handles.behaviorlogic=1;
handles.behaviorvalue2=1;
handles.behaviornormalizenot=0;
handles.behaviorvalue3=1;
handles.features={};
handles.featurelist={};
handles.featurevalue=1;
handles.individuals_behavior=[];
handles.individuals_feature=[];
handles.individuallist={'All'};
handles.individualvalue=1;
handles.individualidx='A';
handles.sexdata={};
handles.fps=nan;
handles.classify_forcecompute=false;
handles.behaviorbarchart_style=1;
handles.behaviortimeseries_style=1;
handles.featurehistogram_style=1;
handles.featurehistogram_style2=1;
handles.comparison=0;
handles.logbinsize=0;
handles.nbins=100;
handles.featuretimeseries_style=1;
handles.featuretimeseries_style2=1;
handles.subtractmean=0;
handles.windowradius=10;
handles.boutstats_style=1;
handles.boutstats_style2=1;
handles.omitnan=1;
handles.omitinf=1;
handles.absdprime=1;
handles.centraltendency=1;
handles.dispersion=1;
handles.xoffset=1;
handles.minimumtrajectorylength=1;
handles.convolutionwidth=10;
handles.pvalue=0.01;
handles.interestingfeaturehistograms_cache=[];
handles.interestingfeaturetimeseries_cache=[];
handles.defaultcolors=[1 0 0;  0 0.5 0;  0 0 1;  0 1 1;  1 0 1;  0.749 0.749 0;  0 0 0];
handles.colors=[];


% ---
function handles=load_configuration_file(filename,hObject,eventdata,handles)

handles_saved=load(filename);
handles_saved=handles_saved.handles;
handles.grouplist=handles_saved.grouplist;
handles.groupvalue=handles_saved.groupvalue;
handles.experimentlist=handles_saved.experimentlist;
handles.experimentvalue=handles_saved.experimentvalue;
try
  handles.classifierlist=handles_saved.classifierlist;
  handles.classifiervalue=handles_saved.classifiervalue;
  handles.configurations=handles_saved.configurations;
  handles.analysis=handles_saved.analysis;
  handles.scorefiles=handles_saved.scorefiles;
  handles.behaviorlist=handles_saved.behaviorlist;
  handles.behaviornot=handles_saved.behaviornot;
  handles.behaviorvalue=handles_saved.behaviorvalue;
  handles.behaviorlogic=handles_saved.behaviorlogic;
  handles.behaviorvalue2=handles_saved.behaviorvalue2;
  handles.behaviornormalizenot=handles_saved.behaviornormalizenot;
  handles.behaviorvalue3=handles_saved.behaviorvalue3;
  handles.features=handles_saved.features;
  handles.featurelist=handles_saved.featurelist;
  handles.featurevalue=handles_saved.featurevalue;
  handles.individuals_behavior=handles_saved.individuals_behavior;
  handles.individuals_feature=handles_saved.individuals_feature;
  handles.individuallist=handles_saved.individuallist;
  handles.individualvalue=handles_saved.individualvalue;
  handles.individualidx=handles_saved.individualidx;
  handles.sexdata=handles_saved.sexdata;
  handles.fps=handles_saved.fps;
  handles.classify_forcecompute=handles_saved.classify_forcecompute;
  handles.behaviorbarchart_style=handles_saved.behaviorbarchart_style;
  handles.behaviortimeseries_style=handles_saved.behaviortimeseries_style;
  handles.featurehistogram_style=handles_saved.featurehistogram_style;
  handles.featurehistogram_style2=handles_saved.featurehistogram_style2;
  handles.comparison=handles_saved.comparison;
  handles.logbinsize=handles_saved.logbinsize;
  handles.nbins=handles_saved.nbins;
  handles.featuretimeseries_style=handles_saved.featuretimeseries_style;
  handles.featuretimeseries_style2=handles_saved.featuretimeseries_style2;
  handles.subtractmean=handles_saved.subtractmean;
  handles.windowradius=handles_saved.windowradius;
  handles.boutstats_style=handles_saved.boutstats_style;
  handles.boutstats_style2=handles_saved.boutstats_style2;
  handles.omitnan=handles_saved.omitnan;
  handles.omitinf=handles_saved.omitinf;
  handles.absdprime=handles_saved.absdprime;
  handles.centraltendency=handles_saved.centraltendency;
  handles.dispersion=handles_saved.dispersion;
  handles.xoffset=handles_saved.xoffset;
  handles.minimumtrajectorylength=handles_saved.minimumtrajectorylength;
  handles.convolutionwidth=handles_saved.convolutionwidth;
  handles.pvalue=handles_saved.pvalue;
  handles.interestingfeaturehistograms_cache=handles_saved.interestingfeaturehistograms_cache;
  handles.interestingfeaturetimeseries_cache=handles_saved.interestingfeaturetimeseries_cache;
  handles.defaultcolors=handles_saved.defaultcolors;
  handles.colors=handles_saved.colors;
  handles.table=[];
catch me
  handles=initialize(handles);
  handles.grouplist=handles_saved.grouplist;
  handles.groupvalue=handles_saved.groupvalue;
  handles.experimentlist=handles_saved.experimentlist;
  handles.experimentvalue=handles_saved.experimentvalue;
  handles.colors=zeros(length(handles.grouplist),3);
  handles=update_experiment_data(handles,true,true,true);
  uiwait(warndlg([filename ' is in an old format.  only group and experiment lists are salvageable.  save a new version']));  drawnow;
end


% ---
function ret_val=get_nindividuals_behavior(experiment_path,score_file,classifier_timestamp)

try
  behavior_data=load(fullfile(experiment_path,score_file));
catch
  behavior_data.allScores.t0s=[];
  behavior_data.timestamp=nan;
end
if ((classifier_timestamp ~= behavior_data.timestamp) || (length(behavior_data.allScores.t0s)==0))
  ret_val=-1;
else
  ret_val=length(behavior_data.allScores.t0s);
end


% ---
function ret_val=get_nindividuals_feature(experiment_path,feature_file)

% could do more error checking here by loading them all in

try
  feature_data=load(fullfile(experiment_path,'perframe',[feature_file{1} '.mat']));
catch
  feature_data.data={};
end
if (length(feature_data.data)==0)
  ret_val=-1;
else
  ret_val=length(feature_data.data);
end


% ---
function ret_val=get_fps(filename)

t=load(filename);
if isfield(t.trx(1),'fps'),
  ret_val=t.trx(1).fps;
elseif isfield(t.trx(1),'dt'),
  fps = 1/mean(t.trx(1).dt(1:10));
  if ~isnan(fps)
    ret_val = fps;
  else
    uiwait(warndlg('Trx file does not have recording fps (frames per second). Assuming fps as 30'));
    drawnow;
    ret_val = 30;
  end
else
  uiwait(warndlg('Trx file does not have recording fps (frames per second). Assuming fps as 30'));
  drawnow;
  ret_val = 30;
end


% ---
function handles=update_experiment_data(handles,features,sexdata,individuals)

handlesexperimentlist=[handles.experimentlist{:}];

handlesfeatures=cell(1,length(handlesexperimentlist));
handlessexdata=cell(1,length(handlesexperimentlist));
handlesindividualsbehavior=zeros(length(handlesexperimentlist),length(handles.scorefiles));
handlesindividualsfeature=zeros(1,length(handlesexperimentlist));

parfor ge=1:length(handlesexperimentlist)
  if(features)
    tmp=dir(fullfile(handlesexperimentlist{ge},'perframe','*.mat'));
    [handlesfeatures{ge}{1:length(tmp)}]=deal(tmp.name);
    handlesfeatures{ge}=cellfun(@(x) x(1:(end-4)),handlesfeatures{ge},'uniformoutput',false);
  end

  if(sexdata)
    fullfile(handlesexperimentlist{ge},'perframe','sex.mat');
    if(exist(ans,'file'))
      tmp=load(ans);
      cellfun(@(x) strcmp(x,'M'),tmp.data,'uniformoutput',false);
      handlessexdata(ge)={ans};
    else
      tmp=dir(fullfile(handlesexperimentlist{ge},'perframe','*.mat'));
      tmp=load(fullfile(handlesexperimentlist{ge},'perframe',tmp(1).name));
      cellfun(@(x) nan(1,length(x)),tmp.data,'uniformoutput',false);
      handlessexdata(ge)={ans};
    end
  end

  if(individuals)
    behavior_data=[];
    parfor_tmp=zeros(1,length(handles.scorefiles));
    for s=1:length(handles.scorefiles)
      classifier=load(handles.classifierlist{s});
      parfor_tmp(s)=get_nindividuals_behavior(handlesexperimentlist{ge},handles.scorefiles{s},...
          classifier.classifierTS);
    end
    handlesindividualsbehavior(ge,:)=parfor_tmp;

    handlesindividualsfeature(ge)=get_nindividuals_feature(handlesexperimentlist{ge},handlesfeatures{ge});
  end
end

if(features)
  handles.features={handlesfeatures{:}};
  handles.featurelist=check_for_diff_and_return_intersection(handles.features);
end

if(sexdata)
  handles.sexdata={handlessexdata{:}};
end

if(individuals)
  handles.individuals_behavior=handlesindividualsbehavior;
  handles.individuals_feature=handlesindividualsfeature;
  handles=fillin_individuallist(handles);
end

%classifier=load(handles.classifierlist{1});
%handles.fps=get_fps(fullfile(handlesexperiments{1},classifier.trxfilename));
handles.fps=get_fps(fullfile(handlesexperimentlist{1},'registered_trx.mat'));

handles.interestingfeaturehistograms_cache=[];
handles.interestingfeaturetimeseries_cache=[];


% ---
function update_figure(handles)

if(isempty(handles.grouplist) || length(handles.experimentlist{handles.groupvalue})==0)
  set(handles.ExperimentList,'enable','off');
else
  set(handles.ExperimentList,'enable','on');
end
if(length(handles.grouplist)<2)
  set(handles.ExperimentMove,'enable','off');
else
  set(handles.ExperimentMove,'enable','on');
end
if(isempty(handles.grouplist))
  set(handles.GroupList,'enable','off');
  set(handles.GroupChange,'enable','off');
  set(handles.ExperimentAdd,'enable','off');
  set(handles.ExperimentDelete,'enable','off');
else
  set(handles.GroupList,'enable','on');
  set(handles.GroupChange,'enable','on');
  set(handles.ExperimentAdd,'enable','on');
  set(handles.ExperimentDelete,'enable','on');
end
if(sum(cellfun(@length,handles.experimentlist))==0)
  set(handles.ClassifierAuto,'enable','off');
  set(handles.ClassifierCheck,'enable','off');
else
  set(handles.ClassifierAuto,'enable','on');
  set(handles.ClassifierCheck,'enable','on');
end
if(isempty(handles.classifierlist))
  set(handles.ClassifierList,'enable','off');
  set(handles.ClassifierDelete,'enable','off');
else
  set(handles.ClassifierList,'enable','on');
  set(handles.ClassifierDelete,'enable','on');
end
if((sum(cellfun(@length,handles.experimentlist))==0) || (isempty(handles.classifierlist)))
  set(handles.ClassifierClassify,'enable','off');
else
  set(handles.ClassifierClassify,'enable','on');
end

if(sum(cellfun(@length,handles.experimentlist))==0)
  set(handles.FeatureHistogram,'enable','off');
  set(handles.FeatureTimeSeries,'enable','off');
  set(handles.InterestingFeatureHistograms,'enable','off');
else
  set(handles.FeatureHistogram,'enable','on');
  set(handles.FeatureTimeSeries,'enable','on');
  set(handles.InterestingFeatureHistograms,'enable','on');
end
if((sum(cellfun(@length,handles.experimentlist))==0) || ...
   (isempty(handles.classifierlist)) || ...
   (sum(sum(handles.individuals_behavior==-1))>0) || ...
   (sum(sum(diff(handles.individuals_behavior,[],2)~=0))>0))
  set(handles.BehaviorBarChart,'enable','off');
  set(handles.BehaviorTimeSeries,'enable','off');
  set(handles.BoutStats,'enable','off');
else
  set(handles.BehaviorBarChart,'enable','on');
  set(handles.BehaviorTimeSeries,'enable','on');
  set(handles.BoutStats,'enable','on');
end

set(handles.StyleList,'enable','off');              set(handles.StyleList2,'enable','off');
set(handles.BehaviorNot,'enable','off');            set(handles.BehaviorList,'enable','off');
set(handles.BehaviorLogic,'enable','off');          set(handles.BehaviorList3,'enable','off');
set(handles.BehaviorNormalizeNot,'enable','off');   set(handles.BehaviorList2,'enable','off');
set(handles.FeatureList,'enable','off');            set(handles.IndividualList,'enable','off');
set(handles.ConvolutionWidth,'enable','off');       set(handles.PValue,'enable','off');
set(handles.NBins,'enable','off');                  set(handles.WindowRadius,'enable','off');
set(handles.XOffset,'enable','off');                set(handles.MinimumTrajectoryLength,'enable','off');
set(handles.OmitInf,'enable','off');                set(handles.OmitNaN,'enable','off');
set(handles.AbsDPrime,'enable','off');              set(handles.SubtractMean,'enable','off');
set(handles.LogBinSize,'enable','off');             set(handles.AllFrames,'enable','off');
set(handles.NotDuring,'enable','off');
set(handles.CentralTendency,'enable','off');
set(handles.Dispersion,'enable','off');

bg=get(handles.GroupNew,'backgroundcolor');
fg=get(handles.GroupNew,'foregroundcolor');
set(handles.FeatureHistogram,'backgroundcolor',bg);
set(handles.FeatureHistogram,'foregroundcolor',fg);
set(handles.FeatureTimeSeries,'backgroundcolor',bg);
set(handles.FeatureTimeSeries,'foregroundcolor',fg);
set(handles.BehaviorBarChart,'backgroundcolor',bg);
set(handles.BehaviorBarChart,'foregroundcolor',fg);
set(handles.BehaviorTimeSeries,'backgroundcolor',bg);
set(handles.BehaviorTimeSeries,'foregroundcolor',fg);
set(handles.BoutStats,'backgroundcolor',bg);
set(handles.BoutStats,'foregroundcolor',fg);
set(handles.InterestingFeatureHistograms,'backgroundcolor',bg);
set(handles.InterestingFeatureHistograms,'foregroundcolor',fg);
analysis2=handles.analysis;
switch(handles.analysis)
  case 'feature_histogram'
    set(handles.FeatureHistogram,'backgroundcolor',fg);
    set(handles.FeatureHistogram,'foregroundcolor',bg);
    set(handles.StyleList,'string',handles.featurehistogram_stylelist,...
        'value',handles.featurehistogram_style2);
    set(handles.StyleList2,'string',handles.featurehistogram_stylelist2,...
        'value',handles.featurehistogram_style);
    if(strcmp(get(handles.FeatureHistogram,'enable'),'off'))
      analysis2='';
    else
      set(handles.FeatureList,'enable','on');
      set(handles.PValue,'enable','on');
      set(handles.CentralTendency,'enable','on');
      set(handles.Dispersion,'enable','on');
      set(handles.LogBinSize,'enable','on');
      set(handles.AllFrames,'enable','on');
      set(handles.NotDuring,'enable','on');
      set(handles.NBins,'enable','on');
    end
  case 'feature_timeseries'
    set(handles.FeatureTimeSeries,'backgroundcolor',fg);
    set(handles.FeatureTimeSeries,'foregroundcolor',bg);
    set(handles.StyleList,'string',handles.featuretimeseries_stylelist,...
        'value',handles.featuretimeseries_style);
    if(length(handles.classifierlist)>0)
      set(handles.StyleList2,'string',handles.featuretimeseries_stylelist2,...
          'value',handles.featuretimeseries_style2);
    else
      set(handles.StyleList2,'string',{''},'value',1);
    end
    if(strcmp(get(handles.FeatureTimeSeries,'enable'),'off'))
      analysis2='';
    else
      set(handles.FeatureList,'enable','on');
      set(handles.ConvolutionWidth,'enable','on');
      set(handles.WindowRadius,'enable','on');
      set(handles.XOffset,'enable','on');
      set(handles.CentralTendency,'enable','on');
      set(handles.Dispersion,'enable','on');
      set(handles.SubtractMean,'enable','on');
    end
  case 'behavior_barchart'
    set(handles.BehaviorBarChart,'backgroundcolor',fg);
    set(handles.BehaviorBarChart,'foregroundcolor',bg);
    set(handles.StyleList,'string',handles.behaviorbarchart_stylelist,...
        'value',handles.behaviorbarchart_style);
    set(handles.StyleList2,'string',{''},'value',1);
    if(strcmp(get(handles.BehaviorBarChart,'enable'),'off'))
      analysis2='';
    else
      set(handles.PValue,'enable','on');
      set(handles.CentralTendency,'enable','on');
      set(handles.Dispersion,'enable','on');
    end
  case 'behavior_timeseries'
    set(handles.BehaviorTimeSeries,'backgroundcolor',fg);
    set(handles.BehaviorTimeSeries,'foregroundcolor',bg);
    set(handles.StyleList,'string',handles.behaviortimeseries_stylelist,...
        'value',handles.behaviortimeseries_style);
    set(handles.StyleList2,'string',{''},'value',1);
    if(strcmp(get(handles.BehaviorTimeSeries,'enable'),'off'))
      analysis2='';
    else
      set(handles.ConvolutionWidth,'enable','on');
      set(handles.XOffset,'enable','on');
      set(handles.CentralTendency,'enable','on');
      set(handles.Dispersion,'enable','on');
    end
  case 'bout_stats'
    set(handles.BoutStats,'backgroundcolor',fg);
    set(handles.BoutStats,'foregroundcolor',bg);
    set(handles.StyleList,'string',handles.boutstats_stylelist,...
        'value',handles.boutstats_style);
    set(handles.StyleList2,'string',handles.boutstats_stylelist2,...
        'value',handles.boutstats_style2);
    if(strcmp(get(handles.BoutStats,'enable'),'off'))
      analysis2='';
    else
      set(handles.PValue,'enable','on');
      set(handles.CentralTendency,'enable','on');
      set(handles.Dispersion,'enable','on');
    end
  case 'interesting_feature_histograms'
    set(handles.InterestingFeatureHistograms,'backgroundcolor',fg);
    set(handles.InterestingFeatureHistograms,'foregroundcolor',bg);
    set(handles.StyleList,'string',{''},'value',1);
    set(handles.StyleList2,'string',{''},'value',1);
    if(strcmp(get(handles.InterestingFeatureHistograms,'enable'),'off'))
      analysis2='';
    else
      set(handles.OmitInf,'enable','on');
      set(handles.OmitNaN,'enable','on');
      set(handles.AbsDPrime,'enable','on');
    end
  otherwise
    set(handles.StyleList,'string',{''});
    set(handles.StyleList2,'string',{''});
end

if(~isempty(analysis2))
  set(handles.MinimumTrajectoryLength,'enable','on');
  if((ismember(analysis2,{'feature_histogram','behavior_barchart','behavior_timeseries','bout_stats'}) || ...
     (strcmp(analysis2,'feature_timeseries')&&(handles.featuretimeseries_style2~=1))) && ...
      (length(handles.behaviorlist)>0))
    set(handles.BehaviorNot,'enable','on');
    set(handles.BehaviorList,'enable','on');
    set(handles.BehaviorLogic,'enable','on');
  end
  if((ismember(analysis2,{'behavior_barchart','behavior_timeseries'})) && ...
      (length(handles.behaviorlist)>0))
    set(handles.BehaviorList3,'enable','on');
  end
  if(~isempty(handles.classifierlist) && (handles.behaviorvalue3>1))
    set(handles.BehaviorNormalizeNot,'enable','on');
  end
  if(handles.behaviorlogic>1)
    set(handles.BehaviorList2,'enable','on');
  end
  if(~strcmp(analysis2,'interesting_feature_histograms'))
    set(handles.IndividualList,'enable','on');
    set(handles.StyleList,'enable','on');
  end
  if(length(get(handles.StyleList2,'string'))>1)
    set(handles.StyleList2,'enable','on');
  end
end

if(isempty(handles.grouplist))
  set(handles.GroupList,'String',{''},'Value',1);
else
  tmp=length(handles.grouplist);
  cellstr(strcat(repmat('<html><font color=#"',tmp,1),...
      reshape(dec2hex(round(handles.colors'*255),2)',6,size(handles.colors,1))',...
      repmat('">',tmp,1),handles.grouplist',repmat('</font></html>',tmp,1)));
  set(handles.GroupList,'String',ans,'Value',handles.groupvalue);
end
if(isempty(handles.experimentlist))
  set(handles.ExperimentList,'String',{''},'Value',1);
else
  set(handles.ExperimentList,'String',handles.experimentlist{handles.groupvalue});
  set(handles.ExperimentList,'Value',handles.experimentvalue{handles.groupvalue});
end
if(isempty(handles.classifierlist))
  set(handles.ClassifierList,'String',{''},'Value',1);
else
  set(handles.ClassifierList,'String',handles.classifierlist);
  set(handles.ClassifierList,'Value',handles.classifiervalue);
end
if(isempty(handles.behaviorlist))
  set(handles.BehaviorList,'String',{''},'Value',1);
  set(handles.BehaviorList2,'String',{''},'Value',1);
  set(handles.BehaviorList3,'String',{''},'Value',1);
else
  set(handles.BehaviorList,'String',{handles.behaviorlist{:} 'All behaviors'},'Value',handles.behaviorvalue);
  set(handles.BehaviorList2,'String',handles.behaviorlist,'Value',handles.behaviorvalue2);
  set(handles.BehaviorList3,'String',{'All frames' handles.behaviorlist{:}},'Value',handles.behaviorvalue3);
end
set(handles.BehaviorNot,'Value',handles.behaviornot);
set(handles.BehaviorLogic,'Value',handles.behaviorlogic);
set(handles.BehaviorNormalizeNot,'Value',handles.behaviornormalizenot);
if(isempty(handles.featurelist))
  set(handles.FeatureList,'String',{''},'Value',1);
else
  set(handles.FeatureList,'String',handles.featurelist,'Value',handles.featurevalue);
end
if(isempty(handles.individuallist))
  set(handles.IndividualList,'String',{''},'Value',1);
else
  set(handles.IndividualList,'String',handles.individuallist,'Value',handles.individualvalue);
end
%set(handles.Table,'Data',[]);

menu_classify_forcecompute_set(handles);
button_comparison_set(handles);

set(handles.MinimumTrajectoryLength,'string',handles.minimumtrajectorylength);
set(handles.ConvolutionWidth,'string',handles.convolutionwidth);
set(handles.PValue,'string',handles.pvalue);
set(handles.NBins,'string',handles.nbins);
set(handles.WindowRadius,'string',handles.windowradius);

set(handles.LogBinSize,'value',handles.logbinsize);
set(handles.SubtractMean,'value',handles.subtractmean);
set(handles.AbsDPrime,'value',handles.absdprime);
set(handles.OmitNaN,'value',handles.omitnan);
set(handles.OmitInf,'value',handles.omitinf);
set(handles.XOffset,'value',handles.xoffset);
set(handles.CentralTendency,'value',handles.centraltendency);
set(handles.Dispersion,'value',handles.dispersion);


% --- Executes just before JAABAPlot is made visible.
function JAABAPlot_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

if ~isdeployed,
  SetUpJAABAPath;
end

if(exist('matlabpool')==2 && matlabpool('size')==0)
  
  useparallel = true;
  if isdeployed,
    if ispc || (isunix && ~ismac),
      filename = deployedRelative2Global('ParallelComputingConfiguration_Local_Win4.settings');
      if ~exist(filename,'file'),
        useparallel = false;
      else
        setmcruserdata('ParallelProfile','ParallelComputingConfiguration_Local_Win4.settings');
      end
    end
  end
  if useparallel
    matlabpool('open');
  end

end

handles.featurehistogram_stylelist=...
    {'Central Tendency','Central Tendency & Dispersion','Overlayed per-Exp Means'};
handles.featurehistogram_stylelist2=...
    {'per Frame','Mean per Bout','Median per Bout','Max per Bout','Min per Bout','Std. Dev. per Bout'};
handles.featuretimeseries_stylelist=...
    {'Central Tendency','Central Tendency & Dispersion','Overlayed per-Exp Means'};
handles.featuretimeseries_stylelist2=...
    {'Entire Recording','Onset Triggered','Offset Triggered'};
handles.behaviorbarchart_stylelist=...
    {'per Group','per Experiment, CT & D','per Experiment, box',...
     'per Fly, grouped','per Fly, scatter','per Fly, trajectory length'};
handles.behaviortimeseries_stylelist=...
    {'Central Tendency','Central Tendency & Dispersion','Overlayed per-Exp Means'};
handles.boutstats_stylelist=...
    {'per Experiment, CT & D', 'per Fly, grouped'};
handles.boutstats_stylelist2=...
    {'Bout Length','Inter-bout Length'};

if isdeployed,
  handles.rcfilename = deployedRelative2Global('most_recent_config.mat');
else
  handles.rcfilename = 'most_recent_config.mat';
end

if(exist(handles.rcfilename)==2)
  try
    handles=load_configuration_file(handles.rcfilename,hObject,eventdata,handles);
  catch ME
    errordlg({'Error loading last used configuration. Using default values.',getReport(ME)},'Error loading last used configuration');
    handles=initialize(handles);
  end
else
  handles=initialize(handles);
end
update_figure(handles);


% Choose default command line output for JAABAPlot
handles.output = hObject;

set(hObject,'CloseRequestFcn',@figure_CloseRequestFcn);
set(handles.ExperimentList,'Callback',@ExperimentListCallback);
set(handles.ClassifierList,'Callback',@ClassifierListCallback);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes JAABAPlot wait for user response (see UIRESUME)
% uiwait(handles.figure1);


function ExperimentListCallback(hObject,eventdata)

handles=guidata(hObject);
handles.experimentvalue{handles.groupvalue}=get(handles.ExperimentList,'Value');
handles.interestingfeaturehistograms_cache=[];
handles.interestingfeaturetimeseries_cache=[];
guidata(hObject,handles);


function ClassifierListCallback(hObject,eventdata)

handles=guidata(hObject);
handles.classifiervalue=get(handles.ClassifierList,'Value');
handles.interestingfeaturehistograms_cache=[];
handles.interestingfeaturetimeseries_cache=[];
guidata(hObject,handles);


% ---
function figure_CloseRequestFcn(hObject, eventdata)

handles=guidata(hObject);
try
  save(handles.rcfilename,'handles');
catch ME,
  errordlg({'Error saving last configuration. State will not be saved.',getReport(ME)},'Error saving last configuration');
end
delete(hObject);


% --- Outputs from this function are returned to the command line.
function varargout = JAABAPlot_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in UpdatePlot.
function UpdatePlot_Callback(hObject, eventdata, handles)
% hObject    handle to UpdatePlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)


% --- Executes on selection change in IndividualList.
function IndividualList_Callback(hObject, eventdata, handles)
% hObject    handle to IndividualList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns IndividualList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from IndividualList

handles.individualvalue=get(handles.IndividualList,'Value');
tmp=cellfun(@(x) x(1),get(handles.IndividualList,'String'));
handles.individualidx=tmp(handles.individualvalue);
if handles.individualidx=='G'
  handles.individualidx=handles.individualvalue-find(tmp=='G',1,'first')+1;
end
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function IndividualList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to IndividualList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in FeatureList.
function FeatureList_Callback(hObject, eventdata, handles)
% hObject    handle to FeatureList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns FeatureList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from FeatureList

handles.featurevalue=get(handles.FeatureList,'Value');
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function FeatureList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FeatureList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in BehaviorNot.
function BehaviorNot_Callback(hObject, eventdata, handles)
% hObject    handle to BehaviorNot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of BehaviorNot

handles.behaviornot=get(handles.BehaviorNot,'Value');
guidata(hObject,handles);


% --- Executes on selection change in BehaviorList.
function BehaviorList_Callback(hObject, eventdata, handles)
% hObject    handle to BehaviorList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns BehaviorList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from BehaviorList

handles.behaviorvalue=get(handles.BehaviorList,'Value');
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function BehaviorList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BehaviorList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in BehaviorLogic.
function BehaviorLogic_Callback(hObject, eventdata, handles)
% hObject    handle to BehaviorLogic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns BehaviorLogic contents as cell array
%        contents{get(hObject,'Value')} returns selected item from BehaviorLogic

handles.behaviorlogic=get(handles.BehaviorLogic,'Value');
if(handles.behaviorlogic==1)
  set(handles.BehaviorList2,'enable','off');
else
  set(handles.BehaviorList2,'enable','on');
end
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function BehaviorLogic_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BehaviorLogic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in BehaviorList2.
function BehaviorList2_Callback(hObject, eventdata, handles)
% hObject    handle to BehaviorList2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns BehaviorList2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from BehaviorList2

handles.behaviorvalue2=get(handles.BehaviorList2,'Value');
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function BehaviorList2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BehaviorList2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in BehaviorNormalizeNot.
function BehaviorNormalizeNot_Callback(hObject, eventdata, handles)
% hObject    handle to BehaviorNormalizeNot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of BehaviorNormalizeNot

handles.behaviornormalizenot=get(handles.BehaviorNormalizeNot,'Value');
guidata(hObject,handles);


% --- Executes on selection change in BehaviorList3.
function BehaviorList3_Callback(hObject, eventdata, handles)
% hObject    handle to BehaviorList3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns BehaviorList3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from BehaviorList3

handles.behaviorvalue3=get(handles.BehaviorList3,'Value');
if(handles.behaviorvalue3==1)
  handles.behaviornormalizenot=0;
  set(handles.BehaviorNormalizeNot,'enable','off');
  update_figure(handles);
else
  set(handles.BehaviorNormalizeNot,'enable','on');
end
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function BehaviorList3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BehaviorList3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ExperimentList.
function ExperimentList_Callback(hObject, eventdata, handles)
% hObject    handle to ExperimentList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ExperimentList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ExperimentList



% --- Executes during object creation, after setting all properties.
function ExperimentList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ExperimentList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% ---
function feature_intersection=check_for_diff_and_return_intersection(arg)

if(length(arg)==0)  feature_intersection=[];  return;  end
if(length(arg)==1)  feature_intersection=arg{1};  return;  end

feature_union=unique([arg{:}]);
feature_intersection=arg{1};
for i=2:length(arg)
  feature_intersection=intersect(feature_intersection,arg{i});
end
if(numel(feature_intersection)<numel(feature_union))
  char(setdiff(feature_union,feature_intersection));
  [ans repmat(', ',size(ans,1),1)];
  reshape(ans',1,numel(ans));
  uiwait(errordlg(['feature(s) ' ans(1:end-2) ' are/is not in all experiments and so will be ignored.']));
  drawnow;
end


% ---
function handles=fillin_individuallist(handles)

%if((numel(handles.individuals_behavior)==0) || ...
%    (sum(sum(diff(handles.individuals_behavior,[],2)~=0))>0))  return;  end

cumsum_num_exp_per_group=[0 cumsum(cellfun(@length,handles.experimentlist))];

tmp(1)={'All'};
if(sum(cellfun(@(x) islogical([x{:}]),handles.sexdata)==0)==0)
  tmp(2:3)={'Male' 'Female'};
end
k=1;
for e=1:length(handles.individuals_feature)
  g=find(cumsum_num_exp_per_group<e,1,'last');
  for i=1:handles.individuals_feature(e)
    tmp{end+1}=['Grp ' handles.grouplist{g} ', Exp #' num2str(e-cumsum_num_exp_per_group(g)) ', Indi #' num2str(i)];
    k=k+1;
  end
end
handles.individuallist=tmp;
if(handles.individualvalue>length(handles.individuallist))
  handles.individualvalue=1;
  handles.individualidx='A';
end
set(handles.IndividualList,'String',handles.individuallist,'Value',handles.individualvalue);


% ---
function handles=experiment_add(handles,newexperiments,newgroups,newcolors)

newexperiments=[cellfun(@(x) regexprep(x,'/$',''),newexperiments,'uniformoutput',false)];
tmp=ismember(newexperiments,[handles.experimentlist{:}]);
if(sum(tmp)>0)
  msg{1}='The following experiments have already been added:';
  msg{2}='';
  msg(3:(2+sum(tmp)))=newexperiments(tmp);
  uiwait(errordlg(msg));
  newexperiments(tmp)=[];
  if(~isempty(newgroups))  newgroups(tmp)=[];  end
  if(~isempty(newcolors))  newcolors(tmp)=[];  end
end
if isempty(newexperiments),
  return;
end

set(handles.Status,'string','Thinking...','foregroundcolor','b');
set(handles.figure1,'pointer','watch');
drawnow;

handlesfeatures=cell(1,length(newexperiments));
handlessexdata=cell(1,length(newexperiments));
handlesindividualsbehavior=zeros(length(newexperiments),length(handles.scorefiles));
handlesindividualsfeature=zeros(1,length(newexperiments));
idx_error=false(1,length(newexperiments));
parfor n=1:length(newexperiments)
%for n=1:length(newexperiments)
  tmp=dir(fullfile(newexperiments{n},'perframe','*.mat'));
  if isempty(tmp),
    handlesfeatures{n} = {};
  else
    [handlesfeatures{n}{1:length(tmp)}]=deal(tmp.name);
    handlesfeatures{n}=cellfun(@(x) x(1:(end-4)),handlesfeatures{n},'uniformoutput',false);
  end

  fullfile(newexperiments{n},'perframe','sex.mat');
  if(exist(ans,'file'))
    tmp=load(ans);
    cellfun(@(x) strcmp(x,'M'),tmp.data,'uniformoutput',false);
    handlessexdata(n)={ans};
  else
    tmp=dir(fullfile(newexperiments{n},'perframe','*.mat'));
    if(~isempty(tmp))
      tmp=load(fullfile(newexperiments{n},'perframe',tmp(1).name));
      cellfun(@(x) nan(1,length(x)),tmp.data,'uniformoutput',false);
      handlessexdata(n)={ans};
    else
      idx_error(n)=true;
      continue;
    end
  end

  behavior_data=[];
  parfor_tmp=zeros(1,length(handles.scorefiles));
  for s=1:length(handles.scorefiles)
    classifier=load(handles.classifierlist{s});
    parfor_tmp(s)=get_nindividuals_behavior(newexperiments{n},handles.scorefiles{s},...
        classifier.classifierTS);
  end
  handlesindividualsbehavior(n,:)=parfor_tmp;

  handlesindividualsfeature(n)=get_nindividuals_feature(newexperiments{n},handlesfeatures{n});
end

if(sum(idx_error)>0)
  uiwait(errordlg([{'skipping the following experiments because they have not been tracked:  ' '' ...
      newexperiments{idx_error}}],''));  drawnow;
  newexperiments=newexperiments(~idx_error);
  handlesfeatures=handlesfeatures(~idx_error);
  handlessexdata=handlessexdata(~idx_error);
  if(~isempty(handlesindividualsbehavior))  handlesindividualsbehavior=handlesindividualsbehavior(~idx_error);  end
  handlesindividualsfeature=handlesindividualsfeature(~idx_error);
  if(~isempty(newgroups))  newgroups=newgroups(~idx_error);  end
  if(~isempty(newcolors))  newcolors=newcolors(~idx_error);  end
end
if isempty(newexperiments),
  return;
end

if(isempty(newgroups))
  handles.experimentlist{handles.groupvalue}={handles.experimentlist{handles.groupvalue}{:} newexperiments{:}};
  handles.experimentvalue{handles.groupvalue}=1:length(handles.experimentlist{handles.groupvalue});
else
  [newgroups,i]=sort(newgroups);
  newexperiments=newexperiments(i);
  if(~isempty(newcolors))  newcolors=newcolors(i);  end
  [ng,ia,ic]=unique(newgroups);
  for ngi=1:length(ng)
    k=length(handles.grouplist);
    handles.grouplist{k+1}=ng{ngi};
    if(~isempty(newcolors))
      handles.colors(k+1,1)=hex2dec(newcolors{ia(ngi)}(1:2))/255;
      handles.colors(k+1,2)=hex2dec(newcolors{ia(ngi)}(3:4))/255;
      handles.colors(k+1,3)=hex2dec(newcolors{ia(ngi)}(5:6))/255;
    else
      handles.colors(k+1,:)=[0 0 0];
    end
    find(cellfun(@(x) strcmp(x,ng{ngi}),newgroups));
    handles.experimentlist{k+1}={newexperiments{ans}};
    handles.experimentvalue{k+1}=1:length(handles.experimentlist{k+1});
  end
end

if(isnan(handles.fps))
%if((isnan(handles.fps))&&(length(handles.classifierlist)>0))
%  classifier=load(handles.classifierlist{1});
%  handles.fps=get_fps(fullfile(newexperiments{1},classifier.trxfilename));
  handles.fps=get_fps(fullfile(newexperiments{1},'registered_trx.mat'));
end

handles.features={handles.features{:} handlesfeatures{:}};
handles.sexdata={handles.sexdata{:} handlessexdata{:}};
handles.individuals_behavior=[handles.individuals_behavior; handlesindividualsbehavior];
handles.individuals_feature=[handles.individuals_feature handlesindividualsfeature];

if(isempty(newgroups))
  tmp=length(handles.features);
  idx=[1 : (sum(cellfun(@length,handles.experimentlist(1:handles.groupvalue)))-length(newexperiments)) ...
      ((tmp-length(newexperiments)+1) : tmp) ...
      ((tmp-length(newexperiments)-sum(cellfun(@length,handles.experimentlist((handles.groupvalue+1):end)))+1) : ...
          (tmp-length(newexperiments)))];
  handles.features=handles.features(idx);
  handles.sexdata=handles.sexdata(idx);
  handles.individuals_behavior=handles.individuals_behavior(idx,:);
  handles.individuals_feature=handles.individuals_feature(idx);
end

handles.featurelist=check_for_diff_and_return_intersection(handles.features);
handles=fillin_individuallist(handles);

handles.interestingfeaturehistograms_cache=[];
handles.interestingfeaturetimeseries_cache=[];

update_figure(handles);

set(handles.Status,'string','Ready.','foregroundcolor','g');
set(handles.figure1,'pointer','arrow');
drawnow;


% --- Executes on button press in ExperimentBatch.
function ExperimentBatch_Callback(hObject, eventdata, handles)
% hObject    handle to ExperimentBatch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

persistent directory
if(isempty(directory))  directory=pwd;  end

tmp=directory;
[newexperiments directory]=uigetfile(fullfile(directory,'*.*'),'Select batch file');
if(isnumeric(newexperiments)&&(newexperiments==0))  directory=tmp; return;  end
newexperiments=textread(fullfile(directory,newexperiments),'%s');
sum(cellfun(@exist,newexperiments)~=0);
if(ans==(length(newexperiments)/2))
  newgroups=newexperiments(2:2:end);
  newcolors=[];
  newexperiments=newexperiments(1:2:end);
elseif(ans==(length(newexperiments)/3))
  newgroups=newexperiments(2:3:end);
  newcolors=newexperiments(3:3:end);
  newexperiments=newexperiments(1:3:end);
else
  newgroups=[];
  newcolors=[];
end
if((length(handles.grouplist)==0)&&(isempty(newgroups)))
  uiwait(errordlg('Add a new group before adding ungrouped experiments'));
  return;
end

handles=experiment_add(handles,newexperiments,newgroups,newcolors);

guidata(hObject,handles);


% --- Executes on button press in ExperimentAdd.
function ExperimentAdd_Callback(hObject, eventdata, handles)
% hObject    handle to ExperimentAdd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

persistent directory
if(isempty(directory))  directory=pwd;  end

newexperiments=uipickfiles('prompt','Select experiment directory','filterspec',directory);
if(~iscell(newexperiments) || (length(newexperiments)==0))  return;  end
[directory,~,~]=fileparts(newexperiments{1});

handles=experiment_add(handles,newexperiments,[],[]);

guidata(hObject,handles);


% --- Executes on button press in ExperimentDelete.
function ExperimentDelete_Callback(hObject, eventdata, handles)
% hObject    handle to ExperimentDelete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

if(length(handles.experimentlist)==0)  return;  end

set(handles.Status,'string','Thinking...','foregroundcolor','b');
set(handles.figure1,'pointer','watch');
drawnow;

idx=handles.experimentvalue{handles.groupvalue};
handles.experimentlist{handles.groupvalue}(idx)=[];
handles.experimentvalue{handles.groupvalue}=1:length(handles.experimentlist{handles.groupvalue});

idx=idx+sum(cellfun(@length,handles.experimentlist(1:(handles.groupvalue-1))));
handles.features(idx)=[];
handles.sexdata(idx)=[];
handles.individuals_behavior(idx,:)=[];
handles.individuals_feature(idx)=[];

if(isempty(handles.experimentlist{handles.groupvalue}))
  handles.colors=handles.colors(setdiff(1:size(handles.colors,1),handles.groupvalue),:);
  handles.experimentlist(handles.groupvalue)=[];
  handles.experimentvalue(handles.groupvalue)=[];
  handles.grouplist(handles.groupvalue)=[];
  handles.groupvalue=max(1,min([handles.groupvalue length(handles.grouplist)]));
end

if(isempty(handles.experimentlist))
  handles.grouplist={};
  handles.groupvalue=1;
  handles.experimentlist={{}};
  handles.experimentvalue={1};
end

handles.featurelist=check_for_diff_and_return_intersection(handles.features);
handles=fillin_individuallist(handles);

update_figure(handles);

handles.interestingfeaturehistograms_cache=[];
handles.interestingfeaturetimeseries_cache=[];

guidata(hObject,handles);

set(handles.Status,'string','Ready.','foregroundcolor','g');
set(handles.figure1,'pointer','arrow');
drawnow;


% --- Executes on button press in ExperimentMove.
function ExperimentMove_Callback(hObject, eventdata, handles)
% hObject    handle to ExperimentMove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if(isempty(handles.experimentlist{handles.groupvalue}))  return;  end
[to_group,ok]=listdlg('liststring',...
    {handles.grouplist{setdiff(1:length(handles.grouplist),handles.groupvalue)}},...
    'selectionmode','single');
if(~ok)  return;  end
if(to_group>=handles.groupvalue)  to_group=to_group+1;  end

set(handles.Status,'string','Thinking...','foregroundcolor','b');
set(handles.figure1,'pointer','watch');
drawnow;

from_group=handles.groupvalue;
idx=handles.experimentvalue{from_group};
idxF=idx+sum(cellfun(@length,handles.experimentlist(1:(from_group-1))));
idxT=sum(cellfun(@length,handles.experimentlist(1:to_group)));

handles.experimentlist{to_group}=...
    {handles.experimentlist{to_group}{:} handles.experimentlist{from_group}{idx}};
handles.experimentvalue{to_group}=1:length(handles.experimentlist{to_group});
handles.experimentlist{from_group}(idx)=[];
handles.experimentvalue{from_group}=1:length(handles.experimentlist{from_group});
%if(isempty(handles.experimentlist{from_group}))
%  set(handles.ExperimentList,'String',{''},'Value',1);
%else
%  set(handles.ExperimentList,'String',handles.experimentlist{from_group});
%  set(handles.ExperimentList,'Value',handles.experimentvalue{from_group});
%end
if(idxF(1)<=idxT)  idxT=idxT-length(idxF);  end
tmp=setdiff(1:length(handles.features),idxF);
tmp=[tmp(1:idxT) idxF tmp((idxT+1):end)];
%handles.behaviors=handles.behaviors(tmp);
handles.features=handles.features(tmp);
handles.sexdata=handles.sexdata(tmp);
handles.individuals_behavior=handles.individuals_behavior(tmp,:);
handles.individuals_feature=handles.individuals_feature(tmp);

handles=fillin_individuallist(handles);

update_figure(handles);

handles.interestingfeaturehistograms_cache=[];
handles.interestingfeaturetimeseries_cache=[];

guidata(hObject,handles);

set(handles.Status,'string','Ready.','foregroundcolor','g');
set(handles.figure1,'pointer','arrow');
drawnow;


% --- Executes on selection change in GroupList.
function GroupList_Callback(hObject, eventdata, handles)
% hObject    handle to GroupList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns GroupList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from GroupList

handles.groupvalue=get(handles.GroupList,'Value');
set(handles.ExperimentList,'String',handles.experimentlist{handles.groupvalue},...
    'Value',handles.experimentvalue{handles.groupvalue});
update_figure(handles);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function GroupList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GroupList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in GroupNew.
function GroupNew_Callback(hObject, eventdata, handles)
% hObject    handle to GroupNew (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

inputdlg({'Name:'},'Create new experiment group');
if((isempty(ans))||strcmp(ans,''))  return;  end
handles.colors(end+1,:)=...
    uisetcolor(handles.defaultcolors(1+mod(length(handles.grouplist),size(handles.defaultcolors,1)),:));
handles.grouplist{end+1}=char(ans);
handles.groupvalue=length(handles.grouplist);
handles.experimentlist{handles.groupvalue}={};
handles.experimentvalue{handles.groupvalue}=[];
update_figure(handles);
guidata(hObject,handles);


% --- Executes on button press in GroupChange.
function GroupChange_Callback(hObject, eventdata, handles)
% hObject    handle to GroupChange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

inputdlg({'Name:'},'New name for this group',1,handles.grouplist(handles.groupvalue));
if((isempty(ans))||strcmp(ans,''))  return;  end
handles.colors(handles.groupvalue,:)=uisetcolor(handles.colors(handles.groupvalue,:));
handles.grouplist{handles.groupvalue}=char(ans);
update_figure(handles);
guidata(hObject,handles);


% ---
function handles=classifier_add(handles,newclassifiers)

handlesconfigurations=cell(1,length(newclassifiers));
handlesbehaviorlist=cell(1,length(newclassifiers));
handlesscorefiles=cell(1,length(newclassifiers));
handlesindividualsbehavior=zeros(sum(cellfun(@length,handles.experimentlist)),length(newclassifiers));
for c=1:length(newclassifiers)
  classifier=load(newclassifiers{c});
  if(~isfield(classifier,'postprocessparams'))
    uiwait(errordlg(['not a valid classifier file.  skipping ' newclassifiers{c}],''));  drawnow;
    newclassifiers{c}='';
    continue;
  end
  params=[];  params.behaviors.names='';
  try
    handlesconfigurations{c}=classifier.configfilename;
    [~,~,ext] = fileparts(classifier.configfilename);
    if strcmpi(ext,'.xml');
      params=ReadXMLConfigParams(handlesconfigurations{c});
    else
      params = load(handlesconfigurations{c});
    end
  catch
    directory = fileparts(newclassifiers{c});
    fullfile(directory,classifier.configfilename);
    if(exist(ans,'file'))
      [pa,na,ex]=fileparts(ans);
    else
      pa=directory;  na='';  ex='';
    end
    [~,tmp,~]=fileparts(newclassifiers{c});
    [configfile tmp]=uigetfile(pa,['Select configuration file for ' tmp],[na ex]);
    if(isnumeric(configfile) && isnumeric(tmp) && (configfile==0) && (tmp==0))
      uiwait(errordlg(['skipping ' newclassifiers{c}],''));
      newclassifiers{c}='';
      continue;
    else
      try
        handlesconfigurations{c}=fullfile(tmp,configfile);
        if strcmpi(ext,'.xml');
          params=ReadXMLConfigParams(handlesconfigurations{c});
        else
          params = load(handlesconfigurations{c});
        end
      catch
        uiwait(errordlg(['problem loading config file.  skipping ' newclassifiers{c}],''));
        newclassifiers{c}='';
        continue;
      end
    end
  end
  if(~isfield(params.behaviors,'names') || ~isfield(params.file,'scorefilename'))
    uiwait(errordlg(['not a valid config file.  skipping ' newclassifiers{c}],''));  drawnow;
    newclassifiers{c}='';
    continue;
  end

  if iscell(params.behaviors.names),
    handlesbehaviorlist{c} = params.behaviors.names{1};
  else
    handlesbehaviorlist{c}=params.behaviors.names;
  end
  handlesscorefiles{c}=params.file.scorefilename;

  ee=0;  behavior_data=[];  parfor_tmp=zeros(sum(cellfun(@length,handles.experimentlist)),1);
  for g=1:length(handles.grouplist)
    for e=1:length(handles.experimentlist{g})
      parfor_tmp(ee+e)=get_nindividuals_behavior(handles.experimentlist{g}{e},handlesscorefiles{c},...
          classifier.classifierTS);
    end
    ee=ee+length(handles.experimentlist{g});
  end
  handlesindividualsbehavior(:,c)=parfor_tmp;
end

handlesexperimentlist=[handles.experimentlist{:}];
if((isnan(handles.fps))&&(length(handlesexperimentlist)>0))
  handles.fps=get_fps(fullfile(handlesexperimentlist{1},classifier.trxfilename));
end

idx=find(~cellfun(@isempty,newclassifiers));
handles.classifierlist={handles.classifierlist{:} newclassifiers{idx}};
handles.configurations={handles.configurations{:} handlesconfigurations{idx}};
handles.behaviorlist={handles.behaviorlist{:} handlesbehaviorlist{idx}};
handles.scorefiles={handles.scorefiles{:} handlesscorefiles{idx}};
handles.individuals_behavior=[handles.individuals_behavior handlesindividualsbehavior(:,idx)];

handles.classifiervalue=1:length(handles.classifierlist);

handles=fillin_individuallist(handles);

handles.interestingfeaturehistograms_cache=[];
handles.interestingfeaturetimeseries_cache=[];


% --- Executes on button press in ClassifierBatch.
function ClassifierBatch_Callback(hObject, eventdata, handles)
% hObject    handle to ClassifierBatch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

persistent directory
if(isempty(directory))  directory=pwd;  end

tmp=directory;
[newclassifiers directory]=uigetfile(fullfile(directory,'*.*'),'Select classifier files');
if(isnumeric(newclassifiers)&&(newclassifiers==0))  directory=tmp; return;  end
newclassifiers=textread(fullfile(directory,newclassifiers),'%s');

set(handles.Status,'string','Thinking...','foregroundcolor','b');
set(handles.figure1,'pointer','watch');
drawnow;

handles=classifier_add(handles,newclassifiers);
update_figure(handles);

guidata(hObject,handles);

set(handles.Status,'string','Ready.','foregroundcolor','g');
set(handles.figure1,'pointer','arrow');
drawnow;


% --- Executes on button press in ClassifierAdd.
function ClassifierAdd_Callback(hObject, eventdata, handles)
% hObject    handle to ClassifierAdd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

persistent directory
if(isempty(directory))  directory=pwd;  end

tmp=directory;
[newclassifiers directory]=uigetfile(directory,'Select classifier files','multiselect','on');
if(isnumeric(newclassifiers)&&(newclassifiers==0))  directory=tmp; return;  end
if(~iscell(newclassifiers))  newclassifiers={newclassifiers};  end
newclassifiers=cellfun(@(x) fullfile(directory,x),newclassifiers,'uniformoutput',false);

set(handles.Status,'string','Thinking...','foregroundcolor','b');
set(handles.figure1,'pointer','watch');
drawnow;

handles=classifier_add(handles,newclassifiers);
update_figure(handles);

guidata(hObject,handles);

set(handles.Status,'string','Ready.','foregroundcolor','g');
set(handles.figure1,'pointer','arrow');
drawnow;


% --- Executes on button press in ClassifierDelete.
function ClassifierDelete_Callback(hObject, eventdata, handles)
% hObject    handle to ClassifierDelete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.Status,'string','Thinking...','foregroundcolor','b');
set(handles.figure1,'pointer','watch');
drawnow;

idx=handles.classifiervalue;
handles.classifierlist(idx)=[];
handles.classifiervalue=1:max(1,length(handles.classifierlist));
handles.configurations(idx)=[];
handles.behaviorlist(idx)=[];
handles.behaviorvalue=1;
handles.behaviorvalue2=1;
handles.behaviorvalue3=1;
handles.scorefiles(idx)=[];
handles.individuals_behavior(:,idx)=[];

%if(length(handles.classifierlist)==0)
%  handles.fps=nan;
%end

handles=fillin_individuallist(handles);
update_figure(handles);

handles.interestingfeaturehistograms_cache=[];
handles.interestingfeaturetimeseries_cache=[];

guidata(hObject,handles);

set(handles.Status,'string','Ready.','foregroundcolor','g');
set(handles.figure1,'pointer','arrow');
drawnow;


% --- Executes on button press in ClassifierAuto.
function ClassifierAuto_Callback(hObject, eventdata, handles)
% hObject    handle to ClassifierAuto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.Status,'string','Thinking...','foregroundcolor','b');
set(handles.figure1,'pointer','watch');
drawnow;

handlesexperimentlist=[handles.experimentlist{:}];

classifiers_found=cell(1,length(handlesexperimentlist));
classifiers_notfound=cell(1,length(handlesexperimentlist));
parfor ge=1:length(handlesexperimentlist)
  tmp=dir(fullfile(handlesexperimentlist{ge},'*.mat'));
  possiblescorefiles=setdiff({tmp.name},handles.scorefiles);
  classifiers_found{ge}={};
  classifiers_notfound{ge}={};
  for p=1:length(possiblescorefiles)
    scoresfilename = fullfile(handlesexperimentlist{ge},possiblescorefiles{p});
    if ~exist(scoresfilename,'file'),
      continue;
    end
    tmp=load(fullfile(handlesexperimentlist{ge},possiblescorefiles{p}));
    if(~isfield(tmp,'classifierfilename'))
      continue;
    end
    if(exist(tmp.classifierfilename)==2)
      classifiers_found{ge}{end+1}=tmp.classifierfilename;
    else
      classifiers_notfound{ge}{end+1}=tmp.classifierfilename;
    end
  end
end
classifiers_found=unique([classifiers_found{:}]);
classifiers_notfound=unique([classifiers_notfound{:}]);

classifiers_found=setdiff(classifiers_found,handles.classifierlist);
handles=classifier_add(handles,classifiers_found);

if(length(classifiers_notfound)>2)
  uiwait(errordlg({'Could not find these classifiers:' '' [classifiers_notfound{:}]}));
end

update_figure(handles);

guidata(hObject,handles);

set(handles.Status,'string','Ready.','foregroundcolor','g');
set(handles.figure1,'pointer','arrow');
drawnow;


% --- Executes on button press in ClassifierCheck.
function ClassifierCheck_Callback(hObject, eventdata, handles)
% hObject    handle to ClassifierCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.Status,'string','Thinking...','foregroundcolor','b');
set(handles.figure1,'pointer','watch');
drawnow;

handlesexperimentlist=[handles.experimentlist{:}];

table=cell(length(handlesexperimentlist),length(handles.scorefiles));

parfor ge=1:length(handlesexperimentlist)
%for ge=1:length(handlesexperimentlist)
  behavior_data=[];
  parfor_tmp=cell(1,length(handles.scorefiles));
  for c=1:length(handles.classifierlist)
    classifier=load(handles.classifierlist{c});
    try
      behavior_data=load(fullfile(handlesexperimentlist{ge},handles.scorefiles{c}));
    catch ME,
%      warning(getReport(ME));
      parfor_tmp{c}='missing';
      continue;
    end
    if (classifier.classifierTS ~= behavior_data.timestamp)
      parfor_tmp{c}=[datestr(classifier.classifierTS) ' ~= ' datestr(behavior_data.timestamp)];
    else
      parfor_tmp{c}='up-to-date';
    end
  end
  table(ge,:)=parfor_tmp;
  tmp=dir(fullfile(handlesexperimentlist{ge},'*.mat'));
  extrascorefiles=setdiff({tmp.name},handles.scorefiles);
  table2{ge}={};
  for s=1:length(extrascorefiles)
    tmp=load(fullfile(handlesexperimentlist{ge},extrascorefiles{s}));
    if(~isfield(tmp,'classifierfilename'))
      continue;
    end
    table2{ge}{end+1}=extrascorefiles{s};
  end
end

[~,tmp,~]=cellfun(@(x) fileparts(x),handlesexperimentlist,'uniformoutput',false);
table=[tmp' table];
[~,tmp,~]=cellfun(@(x) fileparts(x),handles.classifierlist,'uniformoutput',false);
table=[{'' tmp{:}}; table];

tmp=unique([table2{:}]);
if(length(tmp)>0)
  table{1,end+1}='   ';
  table(1,(end+1):(end+length(tmp)))=tmp;
end
for i=1:length(table2)
  for j=1:length(table2{i})
    idx=find(cellfun(@(x) strcmp(x,table2{i}{j}),table(1,:)));
    table{i+1,idx}='extra';
  end
end

if(length(handles.experimentlist)>1)
  %tmp=zeros(1,length(handles.experimentlist)+length([handles.experimentlist{:}]));
  tmp=zeros(1,1+length([handles.experimentlist{:}]));
  cumsum(cellfun(@length,handles.experimentlist));
  tmp(2+ans(1:end-1))=1;
  cumsum(tmp);
  %tmp=ans+(1:(-1+length(handles.experimentlist)+length([handles.experimentlist{:}])));
  tmp=ans+(1:(1+length([handles.experimentlist{:}])));
  tmp2(tmp,:)=[table(:,:)];
  table=tmp2;
end

%figure('menubar','none','toolbar','none','numbertitle','off',...
%    'name','classifier check');
%uitable('units','normalized','position',[0 0 1 1],...
%    'Data',table,'ColumnName',{''},'RowName',[]);
%6*max(cellfun(@length,table),[],1);
%set(t,'ColumnWidth',mat2cell(ans,1,ones(1,length(ans))));
%set(handles.Table,'ColumnWidth','auto');

hf=figure('menubar','none','toolbar','none','numbertitle','off',...
    'name','classifier check');
ht=uitable('data',table,'columnwidth',num2cell(8*max(cellfun(@length,table))),...
    'rowname',[],'columnname',{''},'rowstriping','off');
extT=get(ht,'extent');
posF=get(hf,'position');
%set(hf,'position',[posF(1) posF(2) 16*(posF(4)<extT(4))+extT(3) posF(4)]);
set(hf,'position',[posF(1) posF(2) min(posF(3),16*(posF(4)<extT(4))+extT(3)) min(posF(4),extT(4))]);
%set(ht,'position',[0 0 16*(posF(4)<extT(4))+extT(3) posF(4)]);
set(ht,'units','normalized','position',[0 0 1 1]);

%handles.table_data=table;
%handles.table='classifier';

set(handles.Status,'string','Ready.','foregroundcolor','g');
set(handles.figure1,'pointer','arrow');
drawnow;


% --- Executes on button press in ClassifierClassify.
function ClassifierClassify_Callback(hObject, eventdata, handles)
% hObject    handle to ClassifierClassify (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.Status,'string','Thinking...','foregroundcolor','b');
set(handles.figure1,'pointer','watch');
drawnow;

%h=waitbar(0,'This will likely take awhile...',...
%    'CreateCancelBtn','fid=fopen(fullfile(tempdir,''cancel.txt''),''w''); fclose(fid);');
h=waitbar(0,'This will likely take awhile...');
fid=fopen(fullfile(tempdir,'progressbar.txt'),'w');
fwrite(fid,length(handles.behaviorlist)*sum(cellfun(@length,handles.experimentvalue)),'uint32');
fclose(fid);
t = timer('TimerFcn',{@progress_bar,h}, 'Period', 3, 'ExecutionMode', 'fixedRate');
start(t);

handlesexperimentlist=[handles.experimentlist{:}];

parfor ge=1:length(handlesexperimentlist)
%for ge=1:length(handlesexperimentlist)
  JAABADetect(handlesexperimentlist{ge},...
      'classifierfiles',handles.classifierlist(handles.classifiervalue),...
      'configfiles',handles.configurations(handles.classifiervalue),...
      'forcecompute',handles.classify_forcecompute);
  fid=fopen(fullfile(tempdir,'progressbar.txt'),'a');  fwrite(fid,1,'uint32');  fclose(fid);
end

stop(t);
delete(t);
delete(h);
delete(fullfile(tempdir,'progressbar.txt'));

%if(exist(fullfile(tempdir,'cancel.txt')))
%  delete(fullfile(tempdir,'cancel.txt'));
%  set(handles.Status,'string','Ready.','foregroundcolor','g');
%  drawnow;
%  set(handles.figure1,'pointer','arrow');
%  return;
%end

handles=update_experiment_data(handles,false,false,true);
update_figure(handles);

handles.interestingfeaturehistograms_cache=[];

guidata(hObject,handles);

set(handles.Status,'string','Ready.','foregroundcolor','g');
set(handles.figure1,'pointer','arrow');
drawnow;


% --- Executes on selection change in ClassifierList.
function ClassifierList_Callback(hObject, eventdata, handles)
% hObject    handle to ClassifierList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ClassifierList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ClassifierList


% --- Executes during object creation, after setting all properties.
function ClassifierList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ClassifierList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in FeatureHistogram.
function FeatureHistogram_Callback(hObject, eventdata, handles)
% hObject    handle to FeatureHistogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

handles.analysis='feature_histogram';
update_figure(handles);
guidata(hObject,handles);


% --- Executes on button press in FeatureTimeSeries.
function FeatureTimeSeries_Callback(hObject, eventdata, handles)
% hObject    handle to FeatureTimeSeries (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

handles.analysis='feature_timeseries';
update_figure(handles);
guidata(hObject,handles);


% --- Executes on button press in BehaviorBarChart.
function BehaviorBarChart_Callback(hObject, eventdata, handles)
% hObject    handle to BehaviorBarChart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.analysis='behavior_barchart';
update_figure(handles);
guidata(hObject,handles);


% --- Executes on button press in BehaviorTimeSeries.
function BehaviorTimeSeries_Callback(hObject, eventdata, handles)
% hObject    handle to BehaviorTimeSeries (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.analysis='behavior_timeseries';
update_figure(handles);
guidata(hObject,handles);

%
% --- Executes on button press in BoutStats.
function BoutStats_Callback(hObject, eventdata, handles)
% hObject    handle to BoutStats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.analysis='bout_stats';
update_figure(handles);
guidata(hObject,handles);


% --- Executes on button press in InterestingFeatureHistograms.
function InterestingFeatureHistograms_Callback(hObject, eventdata, handles)
% hObject    handle to InterestingFeatureHistograms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

handles.analysis='interesting_feature_histograms';
update_figure(handles);
guidata(hObject,handles);


% --- Executes on button press in Plot.
function Plot_Callback(hObject, eventdata, handles)
% hObject    handle to Plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

switch(handles.analysis)
  case 'feature_histogram'
    handles=feature_histogram_plot(handles);
  case 'feature_timeseries'
    handles=feature_timeseries_plot(handles);
  case 'behavior_barchart'
    handles=behavior_barchart_plot(handles);
  case 'behavior_timeseries'
    handles=behavior_timeseries_plot(handles);
  case 'bout_stats'
    handles=bout_stats_plot(handles);
  case 'interesting_feature_histograms'
    handles=interesting_feature_histograms_plot(handles);
end

guidata(hObject,handles);


% ---
function ret_val=get_label(feature_name,feature_units)

ret_val=[strrep(char(feature_name),'_','-') ' ('];
if(~isempty(feature_units.num))
  if(length(feature_units.num)>1)
    ret_val=[ret_val '['];
  end
  tmp=char(feature_units.num);
  tmp=[tmp repmat('-',size(tmp,1),1)];
  tmp=reshape(tmp',1,numel(tmp));
  tmp=tmp(1:end-1);
  ret_val=[ret_val tmp];
  if(length(feature_units.num)>1)
    ret_val=[ret_val ']'];
  end
  if(~isempty(feature_units.den))
    ret_val=[ret_val ' / ' ];
  end
end
if(~isempty(feature_units.den))
  if(length(feature_units.den)>1)
    ret_val=[ret_val '['];
  end
  tmp=char(feature_units.den);
  tmp=[tmp repmat('-',size(tmp,1),1)];
  tmp=reshape(tmp',1,numel(tmp));
  tmp=tmp(1:end-1);
  ret_val=[ret_val tmp];
  if(length(feature_units.den)>1)
    ret_val=[ret_val ']'];
  end
  if(isempty(feature_units.num))
    ret_val=[ret_val ' ^ -1 ' ];
  end
end
ret_val=[ret_val ')'];


% ---
function print_csv_help(fid,type,tstr,xstr,ystr)

tstr=strrep(tstr,'%','%%');
xstr=strrep(xstr,'%','%%');
ystr=strrep(ystr,'%','%%');

fprintf(fid,['%% type=' type '\n%% title=' tstr '\n%% xlabel=' xstr '\n%% ylabel=' ystr '\n\n']);


% ---
function print_csv_data(fid,data)

idx=1;
limit=2^14-1;
while idx<length(data)
  fprintf(fid,'%g, ',data(idx:min(idx+limit-1,end)));
  fprintf(fid,'\n');
  idx=idx+limit;
end


% ---
function h=plot_it(ha,xdata,ydata,style,centraltendency,dispersion,color,linewidth,fid,experimentlist)

fprintf(fid,'%% xdata\n');
print_csv_data(fid,xdata);
%fprintf(fid,'%g, ',xdata);  fprintf(fid,'\n');
if(style~=3)
  switch(centraltendency)
    case 1
      data_ct=nanmean(ydata,1);
      str_ct='mean';
    case 2
      data_ct=nanmedian(ydata,1);
      str_ct='median';
    case 3
      data_ct=mode(ydata,1);
      str_ct='mode';
  end
end
switch(style)
  case 1
    h=plot(ha,xdata,data_ct,'color',color,'linewidth',linewidth);
    fprintf(fid,['%% ydata, ' str_ct '\n']);
    print_csv_data(fid,data_ct);
%    fprintf(fid,'%g, ',data_ct);  fprintf(fid,'\n');
  case 2
    switch(dispersion)
      case 1
        tmp=nanstd(ydata,[],1);
        data_dp=data_ct+tmp;
        data_dn=data_ct-tmp;
        str_dp=[str_ct ' + std dev'];
        str_dn=[str_ct ' - std dev'];
      case 2
        tmp=nanstd(ydata,[],1)./sqrt(sum(~isnan(ydata),1));
        data_dp=data_ct+tmp;
        data_dn=data_ct-tmp;
        str_dp=[str_ct ' + std err'];
        str_dn=[str_ct ' - std err'];
      case 3
        data_dp=prctile(ydata,95);
        data_dn=prctile(ydata,5);
        str_dp='95%%';
        str_dn='5%%';
      case 4
        data_dp=prctile(ydata,75);
        data_dn=prctile(ydata,25);
        str_dp='75%%';
        str_dn='25%%';
    end
    h=plot(ha,xdata,data_ct,'color',color,'linewidth',3*linewidth);
    idx=isnan(data_dp) | isnan(data_dn);
    xdata=xdata(~idx);  data_dp=data_dp(~idx);  data_dn=data_dn(~idx);
    color2=(get(h,'color')+[4 4 4])/5;
    k=1;  m=0;  step=10000;
    while(k<=length(xdata))
      idx=k:min(k+step,length(xdata));
      patch([xdata(idx) fliplr(xdata(idx))],[data_dp(idx) fliplr(data_dn(idx))],color2,'edgecolor','none','parent',ha);
      k=k+step+1;  m=m+1;
    end
    get(ha,'children');  set(ha,'children',circshift(ans,-m));  % send to back
    %get(gca,'children');  set(gca,'children',ans([m+1 1:m (m+2):end]));  % send to back
    fprintf(fid,['%% ydata, ' str_dp '\n']);
    print_csv_data(fid,data_dp);
%    fprintf(fid,'%g, ',data_dp);  fprintf(fid,'\n');
    fprintf(fid,['%% ydata, ' str_dn '\n']);
    print_csv_data(fid,data_dn);
%    fprintf(fid,'%g, ',data_dn);  fprintf(fid,'\n');
    fprintf(fid,['%% ydata, ' str_ct '\n']);
    print_csv_data(fid,data_ct);
%    fprintf(fid,'%g, ',data_ct);  fprintf(fid,'\n');
  case 3
    h=plot(ha,xdata,ydata','color',color,'linewidth',linewidth);
    h=h(1);
    for e=1:size(ydata,1)
      fprintf(fid,['%% ydata, experiment ' experimentlist{e} '\n']);
      print_csv_data(fid,ydata(e,:));
%      fprintf(fid,'%g, ',ydata(e,:));  fprintf(fid,'\n');
    end
end
fprintf(fid,'\n');


% --- 
function [during not_during]=calculate_feature_histogram(behavior_data,behavior_logic,behavior_data2,...
    feature_data,sexdata,individual,perwhat,behaviornot)

%if(isempty(behavior_data.allScores.scores))
if(isempty(feature_data.data))
  during=[];
  not_during=[];
  return;
end

if(iscell(feature_data.data{1}))
  vals=unique([feature_data.data{:}]);
  if(length(vals)>2)  error('uhoh');  end
  feature_data.data=cellfun(@(x) strcmp(x,vals{1}),feature_data.data,'uniformoutput',false);
end

if(isempty(behavior_data))
  idx=1:length(feature_data.data);
  if(~isnan(individual))  idx=individual;  end
  [cellfun(@(x,y) x(true & y(1:length(x))), feature_data.data(idx), sexdata(idx), 'uniformoutput',false)];
  during=[ans{:}];
  not_during=[];
  return;
end

during={};  not_during={};
for i=1:length(behavior_data.allScores.t0s)  % individual
  if((~isnan(individual)) && (i~=individual))  continue;  end

  tmp1=zeros(1,length(feature_data.data{i}));
  tmp1(behavior_data.allScores.t0s{i}-behavior_data.allScores.tStart(i)+1)=1;
  tmp1(behavior_data.allScores.t1s{i}-behavior_data.allScores.tStart(i)+1)=-1;
  tmp1=logical(cumsum(tmp1(1:length(feature_data.data{i}))));
  
  if(behavior_logic>1)
    tmp2=zeros(1,length(feature_data.data{i}));
    tmp2(behavior_data2.allScores.t0s{i}-behavior_data2.allScores.tStart(i)+1)=1;
    tmp2(behavior_data2.allScores.t1s{i}-behavior_data2.allScores.tStart(i)+1)=-1;
    tmp2=logical(cumsum(tmp2(1:length(feature_data.data{i}))));
  end

  if(behaviornot)  tmp1=~tmp1;  end

  switch(behavior_logic)
    case(1)
      partition_idx=tmp1;
    case(2)
      partition_idx=tmp1 & tmp2;
    case(3)
      partition_idx=tmp1 & ~tmp2;
    case(4)
      partition_idx=tmp1 | tmp2;
    case(5)
      partition_idx=tmp1 | ~tmp2;
  end
  if length(partition_idx)==(length(feature_data.data{i})+1)  % n-1 features
    partition_idx=partition_idx(1:end-1);
  end

  if(perwhat==1)  % per frame
    during{i}=feature_data.data{i}(partition_idx & sexdata{i}(1:length(partition_idx)));
    not_during{i}=feature_data.data{i}((~partition_idx) & sexdata{i}(1:length(partition_idx)));
  else  % per bout
    partition_idx=[0 partition_idx 0];
    start=find(~partition_idx(1:(end-1)) &  partition_idx(2:end));
    stop =find( partition_idx(1:(end-1)) & ~partition_idx(2:end));
    during{i}=[];  not_during{i}=[];
    if(length(start)>0)
      for j=1:length(start)
        if(sum(sexdata{i}(start(j):(stop(j)-1))) < ((stop(j)-start(j))/2))  continue;  end
        switch(perwhat)
          case 2
            during{i}(j)=mean(feature_data.data{i}(start(j):(stop(j)-1)));
          case 3
            during{i}(j)=median(feature_data.data{i}(start(j):(stop(j)-1)));
          case 4
            during{i}(j)=max(feature_data.data{i}(start(j):(stop(j)-1)));
          case 5
            during{i}(j)=min(feature_data.data{i}(start(j):(stop(j)-1)));
          case 6
            during{i}(j)=std(feature_data.data{i}(start(j):(stop(j)-1)));
        end
      end
      if(length(start)>1)
        for j=1:(length(start)-1)
          if(sum(sexdata{i}(stop(j):(start(j+1)-1))) < ((start(j+1)-stop(j))/2))  continue;  end
          switch(perwhat)
            case 2
              not_during{i}(j)=mean(feature_data.data{i}(stop(j):start(j+1)));
            case 3
              not_during{i}(j)=median(feature_data.data{i}(stop(j):start(j+1)));
            case 4
              not_during{i}(j)=max(feature_data.data{i}(stop(j):start(j+1)));
            case 5
              not_during{i}(j)=min(feature_data.data{i}(stop(j):start(j+1)));
            case 6
              not_during{i}(j)=std(feature_data.data{i}(stop(j):start(j+1)));
          end
        end
      end
    end
  end
end

during=[during{:}];
not_during=[not_during{:}];


% ---
function tmp=calculate_statistics2(table_data,behaviorlist,grouplist,comparison,fid,crit)

nrows=9+3*(comparison>0);
if(((size(table_data{1},1)==2)&&(comparison==0)) || ((size(table_data{1},1)==1)&&(comparison>0)))
  nrows=7+3*(comparison>0);
end

switch(comparison)
  case 0, ctxt='';
  case 1, ctxt=' allframes';
  case 2, ctxt=' notduring';
end

tmp={'behavior'};
for b=1:length(table_data)

  grouplistcopy=grouplist;
  %g=1;
  %while g<=length(table_data{b})
  %  if(isempty(table_data{b}{g}) || (sum(~isnan(table_data{b}{g}))==0))
  %    table_data{b}(g)=[];
  %    grouplistcopy(g)=[];
  %  else
  %    g=g+1;
  %  end
  %end

  tmp{nrows*b+1,1}=behaviorlist{b};
  for g=1:size(table_data{b},1)
    for c=1:(1+(comparison>0))
      idx=2:4;  if(comparison>0)  idx=c+(1:2:5);  end
      p=nan;
      if(sum(~isnan(table_data{b}{g,c}))>0)
        [~,p,~,~]=kstest(table_data{b}{g,c});
      end
      tmp{nrows*b+idx(1),g}=color_red(length(table_data{b}{g,c}),0);
      tmp{nrows*b+idx(2),g}=color_red(p,p<crit);
      tmp{nrows*b+idx(3),g}=color_red(nanstd(table_data{b}{g,c}),0);
      ctxt2='';
      if(comparison>0)
        ctxt2=ctxt;  if(c==1)  ctxt2=' during';  end
      end
      if(b==1)
        tmp(idx,g)=repmat({[grouplistcopy{g} ctxt2]},3,1);
        if(g==1)
          idx=2:4;  if(comparison>0)  idx=[2:3; 4:5; 6:7]';  end
          tmp(idx(:,1),size(table_data{b},1)+1)=repmat({'(sample size)'},size(idx,1),1);
          tmp(idx(:,2),size(table_data{b},1)+1)=repmat({'(K-S normal)'},size(idx,1),1);
          tmp(idx(:,3),size(table_data{b},1)+1)=repmat({'(std. dev.)'},size(idx,1),1);
        end
      end
    end
  end

  if(((size(table_data{1},1)==2)&&(comparison==0)) || ((size(table_data{1},1)==1)&&(comparison>0)))
  
    if(b==1)  tmp{5+3*(comparison>0),1}='t-test';  end
    if(comparison>0)
      [~,p]=ttest2(table_data{b}{1,1},table_data{b}{1,2});
    else
      [~,p]=ttest2(table_data{b}{1,1},table_data{b}{2,1});
    end
    tmp{nrows*b+5+3*(comparison>0),1}=color_red(p,p<crit);

    if(b==1)  tmp{6+3*(comparison>0),1}='Wilcoxen';  end
    if(comparison>0)
      p=ranksum(table_data{b}{1,1},table_data{b}{1,2});
    else
      p=ranksum(table_data{b}{1,1},table_data{b}{2,1});
    end
    tmp{nrows*b+6+3*(comparison>0),1}=color_red(p,p<crit);

  else

    f1=cellfun(@(x,y) y*ones(1,length(x)),table_data{b},...
         num2cell(repmat((1:size(table_data{b},1))',1,size(table_data{b},2))),...
        'uniformoutput',false);
    f2=cellfun(@(x,y) y*ones(1,length(x)),table_data{b},...
         num2cell(repmat((1:size(table_data{b},2)) ,size(table_data{b},1),1)),...
        'uniformoutput',false);
    if(comparison>0)
      [p,table,stats]=anovan([table_data{b}{:}],{[f1{:}] [f2{:}]},'model','interaction','display','off');
    else
      [p,table,stats]=anova1([table_data{b}{:}],[f1{:}],'off');
    end
    c=multcompare(stats,'alpha',crit,'display','off');
    for(pp=1:length(p))
      tmp{nrows*b+5+3*(comparison>0),pp}=color_red(p(pp),p(pp)<crit);
    end
    if(b==1)
      if(comparison>0)
        tmp(8,1:4)={'group','comparison','interaction','(ANOVA)'};
      else
        tmp(5,1)={'(ANOVA)'};
      end
    end
    for g2=1:size(c,1)
      if(b==1)
        tmp((6:8)+3*(comparison>0),g2)=repmat({strcat(grouplistcopy{c(g2,1)},'-',grouplistcopy{c(g2,2)})},3,1);
      end
      foo=(sign(c(g2,3))==sign(c(g2,5)));
      tmp{nrows*b+6+3*(comparison>0),g2}=color_red(c(g2,3),foo);
      tmp{nrows*b+7+3*(comparison>0),g2}=color_red(c(g2,4),foo);
      tmp{nrows*b+8+3*(comparison>0),g2}=color_red(c(g2,5),foo);
    end
    pooh=size(c,1);
    if(comparison>0)
      [p,table,stats]=anovan([table_data{b}{:}],{[f2{:}] [f1{:}]},'model','interaction','display','off');
      c=multcompare(stats,'alpha',crit,'display','off');
      if(b==1)
        tmp((6:8)+3*(comparison>0),pooh+1)=repmat({strcat('during-',ctxt)},3,1);
      end
      foo=(sign(c(1,3))==sign(c(1,5)));
      tmp{nrows*b+6+3*(comparison>0),pooh+1}=color_red(c(1,3),foo);
      tmp{nrows*b+7+3*(comparison>0),pooh+1}=color_red(c(1,4),foo);
      tmp{nrows*b+8+3*(comparison>0),pooh+1}=color_red(c(1,5),foo);
    end
    if(b==1)
      tmp((6:8)+3*(comparison>0),pooh+1+(comparison>0))={'(5%, Tukey post-hoc)'; '(mean)'; '(95%)'};
    end
  end
end
tmp{end+1,1}='';

fprintf(fid,'%% statistics\n');
for j=1:(size(tmp,1)/nrows)
  for i=1:nrows
    if(j==1)  fprintf(fid,'%% ');  end
    tmp(nrows*(j-1)+i,:);
    ans(cellfun(@isempty,ans))='';
    cellfun(@(x) regexprep(x,'<[^<]*>','*'),ans,'uniformoutput',false);
    fprintf(fid,'%s, ',ans{:});
    fprintf(fid,'\n');
  end
end


% ---
function [behavior_data,behavior_data2,behavior_data3,feature_data]=...
    cull_short_trajectories(handles,behavior_data,behavior_data2,behavior_data3,feature_data)

if(~isempty(behavior_data))
  find((behavior_data.allScores.tEnd-behavior_data.allScores.tStart)./handles.fps < handles.minimumtrajectorylength);
  behavior_data.allScores.scores(ans)=[];
  behavior_data.allScores.t0s(ans)=[];
  behavior_data.allScores.t1s(ans)=[];
end
if(~isempty(behavior_data2))
  find((behavior_data2.allScores.tEnd-behavior_data2.allScores.tStart)./handles.fps < handles.minimumtrajectorylength);
  behavior_data2.allScores.scores(ans)=[];
  behavior_data2.allScores.t0s(ans)=[];
  behavior_data2.allScores.t1s(ans)=[];
end
if(~isempty(behavior_data3))
  find((behavior_data3.allScores.tEnd-behavior_data3.allScores.tStart)./handles.fps < handles.minimumtrajectorylength);
  behavior_data3.allScores.scores(ans)=[];
  behavior_data3.allScores.t0s(ans)=[];
  behavior_data3.allScores.t1s(ans)=[];
end
if(~isempty(feature_data))
  find((cellfun(@length,feature_data.data)./handles.fps) < handles.minimumtrajectorylength);
  feature_data.data(ans)=[];
end


% ---
function handles=feature_histogram_plot(handles)

set(handles.Status,'string','Thinking...','foregroundcolor','b');
set(handles.figure1,'pointer','watch');
drawnow;

handlesexperimentlist=[handles.experimentlist{:}];

cumsum_num_indi_per_exp=[0 cumsum(handles.individuals_feature)];
cumsum_num_exp_per_group=[0 cumsum(cellfun(@length,handles.experimentlist))];
mat2cell(cumsum_num_exp_per_group(1:end-1),1,ones(1,length(cumsum_num_exp_per_group)-1));
cellfun(@(x,y) x+y,handles.experimentvalue,ans,'uniformoutput',false);
selected_exp=[ans{:}];
cumsum_num_selexp_per_group=[0 cumsum(cellfun(@length,handles.experimentvalue))];

%ggee=1:length(handlesexperimentlist);
ggee=selected_exp;
individual=handles.individualidx;
if(isnumeric(individual))
  ggee=find(cumsum_num_indi_per_exp<individual,1,'last');
  individual=individual-cumsum_num_indi_per_exp(ggee);
end

fid=fopen('most_recent_figure.csv','w');

handles.type='feature histogram';

hf=figure('toolbar','figure');

bb=handles.behaviorvalue;
if(bb==(length(handles.behaviorlist)+1))  bb=1:(bb-1);  end
if(strcmp(get(handles.BehaviorList,'enable'),'off'))  bb=0;  end

behavior_logic=handles.behaviorlogic;
score_file2=[];
if((length(bb)>1) || (bb>0))  score_file2=handles.scorefiles{handles.behaviorvalue2};  end
feature_value=handles.featurevalue;
feature_list=handles.featurelist;
comparison=handles.comparison;
if((length(bb)==1) && (bb==0))  comparison=0;  end
nbins=handles.nbins;
style=handles.featurehistogram_style2;
centraltendency=handles.centraltendency;
dispersion=handles.dispersion;
behaviornot=handles.behaviornot;

h=[];
table_data={};
for b=bb

  tstr='all frames';
  if(b>0)
    tstr='';  if(behaviornot)  tstr='NOT ';  end
    tstr=[tstr char(strrep(handles.behaviorlist(b),'_','-'))];
  end
  switch(handles.behaviorlogic)
    case 2
      tstr=[tstr ' AND '];
    case 3
      tstr=[tstr ' AND NOT '];
    case 4
      tstr=[tstr ' OR '];
  end
  if(handles.behaviorlogic>1)
    tstr=[tstr char(strrep(handles.behaviorlist(handles.behaviorvalue2),'_','-'))];
  end
  units=load(fullfile(handlesexperimentlist{ggee(1)},'perframe',...
      [feature_list{feature_value} '.mat']),'units');
  xstr=get_label(feature_list(feature_value),units.units);
  ystr='normalized';

  print_csv_help(fid,handles.type,tstr,xstr,ystr);

  if(length(bb)>1)
    ceil(sqrt(length(bb)));
    ha=subplot(ceil(length(bb)/ans),ans,b,'parent',hf);
  else
    ha=subplot(1,1,1,'parent',hf);
  end
  hold(ha,'on');

  score_file=[];
  if(b>0)  score_file=handles.scorefiles{b};  end

  num_indi=0;
  during_data=cell(1,length(ggee));
  not_during_data=cell(1,length(ggee));
  %for gei=1:numel(ggee)
  parfor gei=1:numel(ggee)
    ge = ggee(gei);

    if(b>0)
      behavior_data=load(fullfile(handlesexperimentlist{ge},score_file));
      if(behavior_logic>1)
        behavior_data2=load(fullfile(handlesexperimentlist{ge},score_file2));
      else
        behavior_data2=[];
      end
    else
      behavior_data=[];
      behavior_data2=[];
    end
    feature_data=load(fullfile(handlesexperimentlist{ge},'perframe',...
          [feature_list{feature_value} '.mat']));

    [behavior_data,behavior_data2,~,feature_data]=...
        cull_short_trajectories(handles,behavior_data,behavior_data2,[],feature_data);
    num_indi=num_indi+length(feature_data.data);

    tmp2=handles.sexdata{ge};
    for i=1:length(tmp2)
      switch(individual)
        case('M')
        case('F')
          tmp2{i}=~tmp2{i};
        otherwise
          tmp2{i}=ones(1,length(tmp2{i}));
      end
    end
    tmploop=nan;  if isnumeric(individual)  tmploop=individual;  end

    [during_data{gei} not_during_data{gei}]=calculate_feature_histogram(...
        behavior_data,behavior_logic,behavior_data2,feature_data,tmp2,tmploop,...
        handles.featurehistogram_style,handles.behaviornot);

    if(comparison==1)
      not_during_data{gei}=[during_data{gei} not_during_data{gei}];
    end
  end

  if(num_indi==0)
    delete(hf);
    set(handles.Status,'string','Ready.','foregroundcolor','g');
    set(handles.figure1,'pointer','arrow');
    uiwait(errordlg('no valid data.  check minimum trajectory length.'));  drawnow;
    return;
  end

  max(cellfun(@(x) size(x,2),during_data));
  cellfun(@(x) [x nan(size(x,1),ans-size(x,2))],during_data,'uniformoutput',false);
  during_data=cat(1,ans{:});
  max(cellfun(@(x) size(x,2),not_during_data));
  cellfun(@(x) [x nan(size(x,1),ans-size(x,2))],not_during_data,'uniformoutput',false);
  not_during_data=cat(1,ans{:});

  low=[];  high=[];  nearzero=[];
  if(~isempty(during_data))
    low=min(min(during_data));
    high=max(max(during_data));
    unique(reshape(abs(during_data),1,prod(size(during_data))));
    nearzero=ans(1);  if(ans(1)==0)  nearzero=ans(2);  end
  end
  if((comparison>0) && ~isempty(not_during_data))
    low=min([low min(not_during_data)]);
    high=max([high max(not_during_data)]);
    unique(reshape(abs(not_during_data),1,prod(size(not_during_data))));
    tmp=ans(1);  if(ans(1)==0)  tmp=ans(2);  end
    nearzero=min(tmp,nearzero);
  end

  if(~isempty(low) && ~isempty(high) && ~isempty(nearzero))
    if(handles.logbinsize)
      if((low>=0) && (high>0))
        bins=logspace(log10(max(low,nearzero)),log10(high),nbins);
      elseif((low<0) && (high<=0))
        bins=fliplr(-logspace(log10(max(abs(high),nearzero)),log10(abs(low)),nbins));
      elseif((low<0) && (high>0))
        bins=[fliplr(-logspace(log10(nearzero),log10(abs(low)),nbins)) ...
            logspace(log10(nearzero),log10(abs(high)),nbins)];
      end
      plot(ha,bins,zeros(1,length(bins)),'k.');
    else
      bins=linspace(low,high,nbins);
    end
  end

  for g=1:length(handles.grouplist)
    color=handles.colors(g,:);

    fprintf(fid,['%% group ' handles.grouplist{g} '\n']);

    if ischar(individual)
      idx=(cumsum_num_selexp_per_group(g)+1):(cumsum_num_selexp_per_group(g+1));
    else
      find(cumsum_num_exp_per_group<ggee,1,'last');
      if(ans~=g)  continue;  end
      idx=1;
    end

    if(comparison>0)
      fprintf(fid,['%% not during\n']);
      if(~isempty(not_during_data))
        hist_not_during=hist(not_during_data(idx,:)',bins);
        if(size(hist_not_during,1)==1)  hist_not_during=hist_not_during';  end
        hist_not_during.*repmat(([0 diff(bins)]+[diff(bins) 0])'/2,1,size(hist_not_during,2));
        hist_not_during=hist_not_during./repmat(sum(ans,1),size(hist_not_during,1),1);
        plot_it(ha,bins,hist_not_during',style,centraltendency,dispersion,color,1,...
          fid,handlesexperimentlist(idx));
      end
    end
    fprintf(fid,['%% during\n']);
    if(~isempty(during_data))
      hist_during=hist(during_data(idx,:)',bins);
      if(size(hist_during,1)==1)  hist_during=hist_during';  end
      hist_during.*repmat(([0 diff(bins)]+[diff(bins) 0])'/2,1,size(hist_during,2));
      hist_during=hist_during./repmat(sum(ans,1),size(hist_during,1),1);
      linewidth=1;  if(comparison>0)  linewidth=2;  end
      h(g)=plot_it(ha,bins,hist_during',style,centraltendency,dispersion,color,linewidth,...
          fid,handlesexperimentlist(idx));
    end
  end

  table_data{end+1}={};
  fprintf(fid,'\n%% raw data\n');
  for g=1:length(handles.grouplist)
    if ischar(individual)
      idx=(cumsum_num_selexp_per_group(g)+1):(cumsum_num_selexp_per_group(g+1));
    else
      find(cumsum_num_exp_per_group<ggee,1,'last');
      if(ans~=g)  continue;  end
      idx=1;
    end
    fprintf(fid,['%% group ' handles.grouplist{g} '\n']);
    for e=1:length(idx)
      fprintf(fid,'%% experiment %s\n',handlesexperimentlist{selected_exp(idx(e))});
      if(comparison>0)
        fprintf(fid,'%% not during\n');
        if(~isempty(not_during_data))
          fprint_csv(fid,not_during_data(idx(e),:));
        end
      end
      fprintf(fid,'%% during\n');
      if(~isempty(during_data))
        fprint_csv(fid,during_data(idx(e),:));
      end
      fprintf(fid,'\n');
    end

    if(~isempty(during_data))  table_data{end}{g,1}=nanmean(during_data(idx,:),2)';  end
    if(comparison>0)
      if(~isempty(not_during_data))  table_data{end}{g,2}=nanmean(not_during_data(idx,:),2)';  end
    end
  end

  title(ha,tstr,'interpreter','none');
  xlabel(ha,xstr,'interpreter','none');
  ylabel(ha,ystr,'interpreter','none');
  axis(ha,'tight');  zoom(ha,'reset');
end

idx=find(h>0);
if ischar(individual)
  legend(ha,h(idx),[cellfun(@(x) [x ' ' handles.individuallist{handles.individualvalue}],...
      handles.grouplist,'uniformoutput',false)],'interpreter','none');
else
  legend(ha,h(idx),handles.individuallist(handles.individualvalue),'interpreter','none');
end

uicontrol(hf,'style','pushbutton','string','Params','position',[5 5 60 20],...
    'callback',@figure_params_callback);
if(ischar(individual) && ((length(handles.grouplist)>1) || (comparison>0)) && ...
    (~isempty(during_data)) && (~isempty(not_during_data)))
  uicontrol(hf,'style','pushbutton','string','Stats','position',[70 5 50 20],...
      'callback',@figure_stats_callback);
  if((length(bb)==1) && (bb==0))
    {'all frames'};
  else
    handles.behaviorlist(bb);
  end
  handles.statistics=calculate_statistics2(table_data,ans,handles.grouplist,...
      comparison,fid,handles.pvalue);
end

fclose(fid);

guidata(hf,handles);

set(handles.Status,'string','Ready.','foregroundcolor','g');
set(handles.figure1,'pointer','arrow');
drawnow;


% ---
function progress_bar(~,~,h)

fid=fopen(fullfile(tempdir,'progressbar.txt'),'r');
den=fread(fid,1,'uint32');
fseek(fid,0,1);
num=ftell(fid)/4;
waitbar(num/den,h);
fclose(fid);
drawnow;


% ---
function handles=interesting_feature_histograms_plot(handles)

set(handles.Status,'string','Thinking...','foregroundcolor','b');
set(handles.figure1,'pointer','watch');
drawnow;

if(isempty(handles.interestingfeaturehistograms_cache))
  table_data={};

  cumsum_num_exp_per_group=[0 cumsum(cellfun(@length,handles.experimentlist))];
  mat2cell(cumsum_num_exp_per_group(1:end-1),1,ones(1,length(cumsum_num_exp_per_group)-1));
  cellfun(@(x,y) x+y,handles.experimentvalue,ans,'uniformoutput',false);
  selected_exp=[ans{:}];
  cumsum_num_selexp_per_group=[0 cumsum(cellfun(@length,handles.experimentvalue))];

  handlesexperimentlist=[handles.experimentlist{:}];
  handlesexperimentlist=handlesexperimentlist(selected_exp);

  nexperiments=length(handlesexperimentlist);
  nbehaviors=length(handles.behaviorlist);
  nfeatures=length(handles.featurelist);

  h=waitbar(0,'This will likely take awhile...',...
      'CreateCancelBtn','fid=fopen(fullfile(tempdir,''cancel.txt''),''w''); fclose(fid);');
  fid=fopen(fullfile(tempdir,'progressbar.txt'),'w');
  fwrite(fid,nfeatures*sum(cellfun(@length,handles.experimentvalue)),'uint32');
  fclose(fid);
  t = timer('TimerFcn',{@progress_bar,h}, 'Period', 3, 'ExecutionMode', 'fixedRate');
  start(t);

  num_indi=0;
  table_data=zeros(nexperiments,max(1,nbehaviors),nfeatures,9);
  parfor ge=1:nexperiments
  %for ge=1:nexperiments
    behavior_data={};
    for b=1:nbehaviors
      behavior_data{b}=load(fullfile(handlesexperimentlist{ge},handles.scorefiles{b}));
      [behavior_data{b},~,~,~]=cull_short_trajectories(handles,behavior_data{b},[],[],[]);
      num_indi=num_indi+length(behavior_data{b}.allScores.scores);
    end

    bad{ge}={};
    parfor_tmp=zeros(nbehaviors,nfeatures,9);
    for f=1:nfeatures
      if(exist(fullfile(tempdir,'cancel.txt')))  break;  end
      feature_data=load(fullfile(handlesexperimentlist{ge},'perframe',...
          [handles.featurelist{f} '.mat']));
      [~,~,~,feature_data]=cull_short_trajectories(handles,[],[],[],feature_data);

      sexdata={};
      for s=1:length(feature_data.data)
        sexdata{s}=ones(1,length(feature_data.data{s}));
      end
      for b=1:nbehaviors
        if(exist(fullfile(tempdir,'cancel.txt')))  break;  end

        [during not_during]=calculate_feature_histogram(behavior_data{b},1,[],...
            feature_data,sexdata,nan,handles.featurehistogram_style,0);
        parfor_tmp(b,f,:)=[mean(during) mean(not_during) mean([during not_during]) ...
            std(during) std(not_during) std([during not_during]) ...
            length(during) length(not_during) length([during not_during])];
      end
      if(nbehaviors==0)
        if(exist(fullfile(tempdir,'cancel.txt')))  break;  end

        [during not_during]=calculate_feature_histogram([],1,[],...
            feature_data,sexdata,nan,handles.featurehistogram_style,0);
        parfor_tmp(1,f,:)=[mean(during) mean(not_during) mean([during not_during]) ...
            std(during) std(not_during) std([during not_during]) ...
            length(during) length(not_during) length([during not_during])];
      end
      fid=fopen(fullfile(tempdir,'progressbar.txt'),'a');  fwrite(fid,1,'uint32');  fclose(fid);
    end
    table_data(ge,:,:,:)=parfor_tmp;
  end

  stop(t);
  delete(t);
  delete(h);
  delete(fullfile(tempdir,'progressbar.txt'));

  if(exist(fullfile(tempdir,'cancel.txt')))
    delete(fullfile(tempdir,'cancel.txt'));
    set(handles.Status,'string','Ready.','foregroundcolor','g');
    set(handles.figure1,'pointer','arrow');
    drawnow;
    return;
  end

  if(num_indi==0)
    handles.interestingfeaturehistograms_cache=nan;
  else
    tmp2=[];
    for g=1:length(handles.grouplist)
      gg=cumsum_num_selexp_per_group(g)+(1:length(handles.experimentlist{g}));
      if(nbehaviors>0)
        tmp2=[tmp2; ...
            repmat(g,nbehaviors*nfeatures,1) ...
            reshape(squeeze(sum(table_data(gg,:,:,7))),nbehaviors*nfeatures,1) ...
            repmat(-2,nbehaviors*nfeatures,1) ...
            reshape(squeeze(sum(table_data(gg,:,:,8))),nbehaviors*nfeatures,1) ...
            reshape(repmat(1:nbehaviors,nfeatures,1),nbehaviors*nfeatures,1) ...
            reshape(repmat(1:nfeatures,nbehaviors,1)',nbehaviors*nfeatures,1) ...
            reshape(squeeze((mean(table_data(gg,:,:,1))-mean(table_data(gg,:,:,2)))./ ...
              sqrt((mean(table_data(gg,:,:,4).^2)+mean(table_data(gg,:,:,5).^2))/2))', ...
              nbehaviors*nfeatures,1) ...
            nan(nbehaviors*nfeatures,1)];
        tmp2=[tmp2; ...
            repmat(g,nbehaviors*nfeatures,1) ...
            reshape(squeeze(sum(table_data(gg,:,:,7))),nbehaviors*nfeatures,1) ...
            repmat(-1,nbehaviors*nfeatures,1) ...
            reshape(squeeze(sum(table_data(gg,:,:,9))),nbehaviors*nfeatures,1) ...
            reshape(repmat(1:nbehaviors,nfeatures,1),nbehaviors*nfeatures,1) ...
            reshape(repmat(1:nfeatures,nbehaviors,1)',nbehaviors*nfeatures,1) ...
            reshape(squeeze((mean(table_data(gg,:,:,1))-mean(table_data(gg,:,:,3)))./ ...
              sqrt((mean(table_data(gg,:,:,4).^2)+mean(table_data(gg,:,:,6).^2))/2))', ...
              nbehaviors*nfeatures,1) ...
            nan(nbehaviors*nfeatures,1)];
      end

      if(g==length(handles.grouplist))  break;  end
      for g2=(g+1):length(handles.grouplist)
        gg2=cumsum_num_selexp_per_group(g2)+(1:length(handles.experimentlist{g2}));
        if(nbehaviors>0)
          tmp2=[tmp2; ...
              repmat(g,nbehaviors*nfeatures,1) ...
              reshape(squeeze(sum(table_data(gg,:,:,7))),nbehaviors*nfeatures,1) ...
              repmat(g2,nbehaviors*nfeatures,1) ...
              reshape(squeeze(sum(table_data(gg2,:,:,7))),nbehaviors*nfeatures,1) ...
              reshape(repmat(1:nbehaviors,nfeatures,1),nbehaviors*nfeatures,1) ...
              reshape(repmat(1:nfeatures,nbehaviors,1)',nbehaviors*nfeatures,1) ...
              reshape(squeeze((mean(table_data(gg,:,:,1))-mean(table_data(gg2,:,:,1)))./ ...
                sqrt((mean(table_data(gg,:,:,4).^2)+mean(table_data(gg2,:,:,4).^2))/2))', ...
                nbehaviors*nfeatures,1) ...
              reshape(squeeze((mean(table_data(gg,:,:,3))-mean(table_data(gg2,:,:,3)))./ ...
                sqrt((mean(table_data(gg,:,:,6).^2)+mean(table_data(gg2,:,:,6).^2))/2))', ...
                nbehaviors*nfeatures,1)];
          tmp2=[tmp2; ...
              repmat(g,nbehaviors*nfeatures,1) ...
              reshape(squeeze(sum(table_data(gg,:,:,8))),nbehaviors*nfeatures,1) ...
              repmat(g2,nbehaviors*nfeatures,1) ...
              reshape(squeeze(sum(table_data(gg2,:,:,8))),nbehaviors*nfeatures,1) ...
              reshape(repmat(-(1:nbehaviors),nfeatures,1),nbehaviors*nfeatures,1) ...
              reshape(repmat(1:nfeatures,nbehaviors,1)',nbehaviors*nfeatures,1) ...
              reshape(squeeze((mean(table_data(gg,:,:,2))-mean(table_data(gg2,:,:,2)))./ ...
                sqrt((mean(table_data(gg,:,:,5).^2)+mean(table_data(gg2,:,:,5).^2))/2))', ...
                nbehaviors*nfeatures,1) ...
              reshape(squeeze((mean(table_data(gg,:,:,3))-mean(table_data(gg2,:,:,3)))./ ...
                sqrt((mean(table_data(gg,:,:,6).^2)+mean(table_data(gg2,:,:,6).^2))/2))', ...
                nbehaviors*nfeatures,1)];
        end
        tmp2=[tmp2; ...
            repmat(g,nfeatures,1) ...
            squeeze(sum(table_data(gg,1,:,9))) ...
            repmat(g2,nfeatures,1) ...
            squeeze(sum(table_data(gg2,1,:,9))) ...
            zeros(nfeatures,1) ...
            (1:nfeatures)' ...
            squeeze((mean(table_data(gg,1,:,3))-mean(table_data(gg2,1,:,3)))./ ...
              sqrt((mean(table_data(gg,1,:,6).^2)+mean(table_data(gg2,1,:,6).^2))/2)), ...
            nan(nfeatures,1)];
      end
    end
    handles.interestingfeaturehistograms_cache=tmp2;
  end
else
  tmp2=handles.interestingfeaturehistograms_cache;
end

if(isnan(handles.interestingfeaturehistograms_cache))
  set(handles.Status,'string','Ready.','foregroundcolor','g');
  set(handles.figure1,'pointer','arrow');
  uiwait(errordlg('no valid data.  check minimum trajectory length.'));  drawnow;
  return;
end

if(handles.omitnan)
  idx=find(~isnan(tmp2(:,7)));
  tmp2=tmp2(idx,:);
end
if(handles.omitinf)
  idx=find(~isinf(tmp2(:,7)));
  tmp2=tmp2(idx,:);
end
if(handles.absdprime)
  tmp2(:,7)=abs(tmp2(:,7));
  tmp2(:,8)=abs(tmp2(:,8));
end
tmp2=sortrows(tmp2,-7);

tmp=cell(size(tmp2,1),8);
tmp(:,1)=handles.grouplist(tmp2(:,1));
tmp(:,2)=cellstr(num2str(tmp2(:,2),'%-d'));
idx=(tmp2(:,3)>0);    tmp(idx,3)=handles.grouplist(tmp2(idx,3));
idx=(tmp2(:,3)==-1);  tmp(idx,3)=cellstr('all frames');
idx=(tmp2(:,3)==-2);  tmp(idx,3)=cellstr('not during');
tmp(:,4)=cellstr(num2str(tmp2(:,4),'%-d'));
idx=(tmp2(:,5)>0);   tmp(idx,5)=handles.behaviorlist(tmp2(idx,5));
idx=(tmp2(:,5)==0);  tmp(idx,5)=cellstr('all frames');
idx=(tmp2(:,5)<0);   tmp(idx,5)=cellstr([repmat('not ',sum(idx),1) char(handles.behaviorlist(-tmp2(idx,5)))]);
tmp(:,6)=handles.featurelist(tmp2(:,6));
tmp(:,7)=cellstr(num2str(tmp2(:,7)));
idx=find(~isnan(tmp2(:,8)));
tmp(idx,8)=cellstr(num2str(tmp2(idx,8)));

handles2.figure1=handles.figure1;
handles2.figure2=figure('menubar','none','toolbar','none','numbertitle','off',...
    'name','interesting feature histograms');
handles2.Table=uitable('Data',tmp,...
    'ColumnName',{'Group' 'n' 'Group2' 'n2' 'Behavior' 'Feature' 'd''' 'd''-AF'},...
    'ColumnWidth',num2cell(7*max(cellfun(@length,tmp))),...
    'RowName',[],'RowStriping','on','BackgroundColor',[1 1 1; 0.95 0.95 0.95],...
    'CellSelectionCallback',@CellSelectionCallback);
extT=get(handles2.Table,'extent');
posF=get(handles2.figure2,'position');
%set(hf,'position',[posF(1) posF(2) 16*(posF(4)<extT(4))+extT(3) posF(4)]);
set(handles2.figure2,'position',[posF(1) posF(2) min(posF(3),16*(posF(4)<extT(4))+extT(3)) min(posF(4),extT(4))]);
%set(ht,'position',[0 0 16*(posF(4)<extT(4))+extT(3) posF(4)]);
set(handles2.Table,'units','normalized','position',[0 0 1 1]);
handles2.table_data=tmp2;
handles2.table='histogram';
guidata(handles2.figure2,handles2);

set(handles.Status,'string','Ready.','foregroundcolor','g');
set(handles.figure1,'pointer','arrow');
drawnow;


% --- Executes on button press in InterestingFeatureTimeSeries.
function InterestingFeatureTimeSeries_Callback(hObject, eventdata, handles)
% hObject    handle to InterestingFeatureTimeSeries (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- 
function data=calculate_entiretimeseries(behavior_data,feature_data,sexdata,individual,xoffset);

if(isempty(feature_data.data))
  data=[];
  return;
end

if(iscell(feature_data.data{1}))
  vals=unique([feature_data.data{:}]);
  if(length(vals)>2)  error('uhoh');  end
  feature_data.data=cellfun(@(x) strcmp(x,vals{1}),feature_data.data,'uniformoutput',false);
end

if(isempty(behavior_data))
  behavior_data.allScores.tStart=1;
  behavior_data.allScores.tEnd=max(cellfun(@length,feature_data.data))+1;
end

k=1;
if(xoffset==1)
  data=nan(length(feature_data.data),max(behavior_data.allScores.tEnd));
else
  data=nan(length(feature_data.data),max(behavior_data.allScores.tEnd)-min(behavior_data.allScores.tStart));
end
for i=1:length(feature_data.data)  % individual
  if((~isnan(individual)) && (i~=individual))  continue;  end
  if(sum(sexdata{i}) < length(sexdata{i})/2)  continue;  end
  switch(xoffset)
    case(1)
      behavior_data.allScores.tStart(i);
    case(2)
      1;
    case(3)
      behavior_data.allScores.tStart(i)-min(behavior_data.allScores.tStart)+1;
  end
  data(k,ans:(ans+length(feature_data.data{i})-1))=feature_data.data{i};
  k=k+1;
end


% --- 
function triggered_data=calculate_triggeredtimeseries(behavior_data,behavior_logic,behavior_data2, ...
    feature_data,sexdata,individual,timing,windowradius,subtractmean,behaviornot)

timing=2+xor(timing-2,behaviornot);

if(iscell(feature_data.data{1}))
  vals=unique([feature_data.data{:}]);
  if(length(vals)>2)  error('uhoh');  end
  feature_data.data=cellfun(@(x) strcmp(x,vals{1}),feature_data.data,'uniformoutput',false);
end

k=1;
triggered_data=[];
for i=1:length(behavior_data.allScores.t0s)  % individual
  if((~isnan(individual)) && (i~=individual))  continue;  end
  feature_data_padded=[nan(1,windowradius) feature_data.data{i} nan(1,windowradius)];
  sexdata_padded=[nan(1,windowradius) sexdata{i} nan(1,windowradius)];
  if(behavior_logic>1)
    idx2=[behavior_data2.allScores.t0s{i}'-behavior_data2.allScores.tStart(i) ...
       behavior_data2.allScores.t1s{i}'-behavior_data2.allScores.tStart(i)];
  end
  for j=1:length(behavior_data.allScores.t0s{i})  % bout
    switch(timing)
      case(2)
        idx=behavior_data.allScores.t0s{i}(j)-behavior_data.allScores.tStart(i);
      case(3)
        idx=behavior_data.allScores.t1s{i}(j)-behavior_data.allScores.tStart(i);
    end
    if((behavior_logic==2) && (diff(sum(idx2<idx))==0))   continue;  end
    if((behavior_logic==3) && (diff(sum(idx2<idx))==-1))  continue;  end
    if(behavior_logic>3) error('or & or not are not implemented yet');  end
    if(sum(sexdata_padded(idx+(0:(2*windowradius))+(3-timing))) < windowradius)  continue;  end
    triggered_data(k,:)=feature_data_padded(idx+(0:(2*windowradius))+(3-timing));
    k=k+1;
  end
end
if(k==1)
  triggered_data(k,:)=nan(1,2*windowradius+1);
end
if(subtractmean)
  triggered_data=triggered_data-...
      repmat(nanmean(triggered_data(:,1:windowradius),2),1,size(triggered_data,2));
end


% ---
function handles=feature_timeseries_plot(handles)

set(handles.Status,'string','Thinking...','foregroundcolor','b');
set(handles.figure1,'pointer','watch');
drawnow;

handlesexperimentlist=[handles.experimentlist{:}];

cumsum_num_indi_per_exp=[0 cumsum(handles.individuals_feature)];
cumsum_num_exp_per_group=[0 cumsum(cellfun(@length,handles.experimentlist))];
mat2cell(cumsum_num_exp_per_group(1:end-1),1,ones(1,length(cumsum_num_exp_per_group)-1));
cellfun(@(x,y) x+y,handles.experimentvalue,ans,'uniformoutput',false);
selected_exp=[ans{:}];
cumsum_num_selexp_per_group=[0 cumsum(cellfun(@length,handles.experimentvalue))];

%ggee=1:length(handlesexperimentlist);
ggee=selected_exp;
individual=handles.individualidx;
if isnumeric(individual)
  ggee=find(cumsum_num_indi_per_exp<individual,1,'last');
  individual=individual-cumsum_num_indi_per_exp(ggee);
end

fid=fopen('most_recent_figure.csv','w');

handles.type='feature time series';

hf=figure('toolbar','figure');

bb=handles.behaviorvalue;
if(handles.featuretimeseries_style2==1)  bb=1;  end
if(bb==(length(handles.behaviorlist)+1))  bb=1:(bb-1);  end
if(strcmp(get(handles.BehaviorList,'enable'),'off'))  bb=0;  end

behavior_logic=handles.behaviorlogic;
score_file2=[];
if((length(bb)>1) || (bb>0))  score_file2=handles.scorefiles{handles.behaviorvalue2};  end
feature_value=handles.featurevalue;
feature_list=handles.featurelist;
sexdata=handles.sexdata;
timing=handles.featuretimeseries_style2;
xoffset=handles.xoffset;
if((length(bb)==1) && (bb==0))
  timing=1;
  xoffset=2;
end
style=handles.featuretimeseries_style;
centraltendency=handles.centraltendency;
dispersion=handles.dispersion;
convolutionwidth=round(handles.convolutionwidth*handles.fps);
subtractmean=handles.subtractmean;
windowradius=handles.windowradius;
behaviornot=handles.behaviornot;

h=[];
for b=bb

  if(length(bb)>1)
    ceil(sqrt(length(bb)));
    ha=subplot(ceil(length(bb)/ans),ans,b,'parent',hf);
  else
    ha=subplot(1,1,1,'parent',hf);
  end
  hold(ha,'on');

  score_file=[];
  if(b>0)  score_file=handles.scorefiles{b};  end

  num_indi=0;
  raw_data=cell(1,length(handlesexperimentlist));
  data=cell(1,length(handlesexperimentlist));
  parfor gei=1:length(ggee)
  %for gei=1:length(ggee)
    ge = ggee(gei);

    if(b>0)
      behavior_data=load(fullfile(handlesexperimentlist{ge},score_file));
      if(behavior_logic>1)
        behavior_data2=load(fullfile(handlesexperimentlist{ge},score_file2));
      else
        behavior_data2=[];
      end
    else
      behavior_data=[];
      behavior_data2=[];
    end
    feature_data=load(fullfile(handlesexperimentlist{ge},'perframe',...
        [feature_list{feature_value} '.mat']));

    [behavior_data,behavior_data2,~,feature_data]=...
        cull_short_trajectories(handles,behavior_data,behavior_data2,[],feature_data);
    num_indi=num_indi+length(feature_data.data);

    tmp2=sexdata{ge};
    for i=1:length(tmp2)
      switch(individual)
        case('M')
        case('F')
          tmp2{i}=~tmp2{i};
        otherwise
          tmp2{i}=ones(1,length(tmp2{i}));
      end
    end
    tmp=nan;  if isnumeric(individual)  tmp=individual;  end

    if(timing==1)
      calculate_entiretimeseries(behavior_data,feature_data,tmp2,tmp,xoffset);
      raw_data{gei}=nanmean(ans,1);
      if(~isempty(raw_data{gei}))
        conv(raw_data{gei},ones(1,convolutionwidth),'same');
        data{gei}=ans./conv(ones(1,length(ans)),ones(1,convolutionwidth),'same');
      else
        data{gei}=raw_data{gei};
      end
    else
      calculate_triggeredtimeseries(behavior_data,behavior_logic,behavior_data2,...
          feature_data,tmp2,tmp,timing,windowradius,subtractmean,behaviornot);
      raw_data{gei}=nanmean(ans,1);
      data{gei}=raw_data{gei};
    end
  end

  if(num_indi==0)
    delete(hf);
    set(handles.Status,'string','Ready.','foregroundcolor','g');
    set(handles.figure1,'pointer','arrow');
    uiwait(errordlg('no valid data.  check minimum trajectory length.'));  drawnow;
    return;
  end

  if(timing==1)
    max(cellfun(@(x) size(x,2),data));
    cellfun(@(x) [x nan(size(x,1),ans-size(x,2))],data,'uniformoutput',false);
    ydata=cat(1,ans{:});
    xdata=1:size(ydata,2);
  else
    ydata=cat(1,data{:});
    xdata=-windowradius:windowradius;
  end

  tstr='';
  if(timing>1)
    tstr='';  if(behaviornot)  tstr='NOT ';  end
    tstr=[tstr char(strrep(handles.behaviorlist(b),'_','-'))];
    switch(handles.behaviorlogic)
      case 2
        tstr=[tstr ' AND '];
      case 3
        tstr=[tstr ' AND NOT '];
      case 4
        tstr=[tstr ' OR '];
      case 5
        tstr=[tstr ' OR NOT '];
    end
    if(behavior_logic>1)
      tstr=[tstr char(strrep(handles.behaviorlist(handles.behaviorvalue2),'_','-'))];
    end
  end
  time_base=xdata./handles.fps;
  xstr='time (sec)';
  if(time_base(end)>300)
    time_base=time_base./60;
    xstr='time (min)';
  end
  if(time_base(end)>300)
    time_base=time_base./60;
    xstr='time (hr)';
  end
  if(time_base(end)>120)
    time_base=time_base./24;
    xstr='time (d)';
  end
  units=load(fullfile(handlesexperimentlist{ggee(1)},'perframe',...
      [feature_list{feature_value} '.mat']),'units');
  ystr=get_label(feature_list(feature_value),units.units);

  print_csv_help(fid,handles.type,tstr,xstr,ystr);

  for g=1:length(handles.grouplist)
    color=handles.colors(g,:);

    if ischar(individual)
      idx=(cumsum_num_selexp_per_group(g)+1):(cumsum_num_selexp_per_group(g+1));
    else
      find(cumsum_num_exp_per_group<ggee,1,'last');
      if(ans~=g)  continue;  end
      idx=1;
    end

    fprintf(fid,['%% group ' handles.grouplist{g} '\n']);
    h(g)=plot_it(ha,time_base,ydata(idx,:),style,centraltendency,dispersion,color,1,fid,handlesexperimentlist(idx));
  end

  fprintf(fid,'\n%% raw data\n');
  for g=1:length(handles.grouplist)
    if ischar(individual)
      idx=(cumsum_num_selexp_per_group(g)+1):(cumsum_num_selexp_per_group(g+1));
    else
      find(cumsum_num_exp_per_group<ggee,1,'last');
      if(ans~=g)  continue;  end
      idx=1;
    end
    fprintf(fid,['%% group ' handles.grouplist{g} '\n']);
    for e=1:length(idx)
      fprintf(fid,'%% experiment %s\n',handlesexperimentlist{selected_exp(idx(e))});
      print_csv_data(fid,raw_data{idx(e)});
      fprintf(fid,'\n');
    end
  end

  xlabel(ha,xstr,'interpreter','none');
  ylabel(ha,ystr,'interpreter','none');
  title(ha,tstr,'interpreter','none');
  axis(ha,'tight');;  zoom(ha,'reset');

end
idx=find(h>0);
if ischar(individual)
  legend(ha,h(idx),[cellfun(@(x) [x ' ' handles.individuallist{handles.individualvalue}],...
      handles.grouplist,'uniformoutput',false)],'interpreter','none');
else
  legend(ha,h(idx),handles.individuallist(handles.individualvalue),'interpreter','none');
end

uicontrol(hf,'style','pushbutton','string','Params','position',[5 5 60 20],...
    'callback',@figure_params_callback);

fclose(fid);

guidata(hf,handles);

set(handles.Status,'string','Ready.','foregroundcolor','g');
set(handles.figure1,'pointer','arrow');
drawnow;


% ---
function table_data=calculate_interesting_timeseries(experiment_value,experiment_list,...
    behavior_list,feature_list,windowradius)

table_data=cell(length(behavior_list),length(feature_list),2);
parfor b=1:length(behavior_list)
%for b=1:length(behavior_list)
  parfor_tmp=cell(length(feature_list),2);
  for f=1:length(feature_list)
    for t=2:3  % timing
      parfor_tmp{f,t-1}=[];
      for e=1:length(experiment_value)
        behavior_data=load(fullfile(experiment_list{experiment_value(e)},[behavior_list{b} '.mat']));
        feature_data=load(fullfile(experiment_list{experiment_value(e)},'perframe',...
          [feature_list{f} '.mat']));
        sexdata={};
        for i=1:length(feature_data.data)
          sexdata{i}=ones(1,length(feature_data.data{i}));
        end

        tmp=calculate_triggeredtimeseries(behavior_data,1,[],feature_data,sexdata,nan,t,windowradius,1);
        parfor_tmp{f,t-1}=[parfor_tmp{f,t-1}; tmp];
      end
      parfor_tmp{f,t-1}=[...
          sqrt(nanmean(parfor_tmp{f,t-1}(:,1:windowradius).^2,2))...
          sqrt(nanmean(parfor_tmp{f,t-1}(:,(windowradius+1):end).^2,2))];
    end
  end
  table_data(b,:,:)=parfor_tmp;
  disp([num2str(b) ' of ' num2str(length(behavior_list))]);
end


% --- Executes on button press in InterestingTimeSeries.
function InterestingTimeSeries_Callback(hObject, eventdata, handles)
% hObject    handle to InterestingTimeSeries (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

set(handles.Status,'string','Thinking...','foregroundcolor','b');
set(handles.figure1,'pointer','watch');
drawnow;

experiment_value=get(handles.ExperimentList,'Value');
experiment_list=get(handles.ExperimentList,'String');
experiment_value2=get(handles.ExperimentList2,'Value');
experiment_list2=get(handles.ExperimentList2,'String');
behavior_list=get(handles.BehaviorList,'String');
feature_list=get(handles.FeatureList,'String');
windowradius=handles.windowradius;

if(isempty(handles.interestingfeaturetimeseries_cache))
  if(length(experiment_value)>0)
    table_data=calculate_interesting_timeseries(experiment_value,experiment_list,...
        behavior_list,feature_list,windowradius);
  end
  if(length(experiment_value2)>0)
    table_data2=calculate_interesting_timeseries(experiment_value2,experiment_list2,...
        behavior_list,feature_list,windowradius);
  end

  if((length(experiment_value)>0) && (length(experiment_value2)>0))
    tmp=(cellfun(@(x) nanmean(x(:,2)),table_data2) - cellfun(@(x) nanmean(x(:,2)),table_data)) ./ ...
        sqrt(cellfun(@(x) nanstd(x(:,2)),table_data2).^2 + cellfun(@(x) nanstd(x(:,2)),table_data).^2);
  elseif(length(experiment_value)>0)
    tmp=cellfun(@(x) (nanmean(x(:,2))-nanmean(x(:,1)))/sqrt(nanstd(x(:,2)).^2+nanstd(x(:,1)).^2),table_data);
  elseif(length(experiment_value2)>0)
    tmp=cellfun(@(x) (nanmean(x(:,2))-nanmean(x(:,1)))/sqrt(nanstd(x(:,2)).^2+nanstd(x(:,1)).^2),table_data2);
  end
  [tmp2(:,1) tmp2(:,2) tmp2(:,3)]=ind2sub(size(tmp),1:prod(size(tmp)));
  tmp2(:,4)=tmp(1:prod(size(tmp)));
  tmp2=sortrows(tmp2,-4);

  handles.interestingfeaturetimeseries_cache=tmp2;
else
  tmp2=handles.interestingfeaturetimeseries_cache;
end

tmp=cell(size(tmp2));
tmp(:,1)=behavior_list(tmp2(:,1));
tmp(:,2)=feature_list(tmp2(:,2));
tmp(:,3)=repmat({'Onset'},size(tmp2,1),1);
  find(tmp2(:,3)==2);
  tmp(ans,3)=repmat({'Offset'},size(ans),1);
tmp(:,4)=num2cell(tmp2(:,4));

set(handles.Table,'Data',tmp);
set(handles.Table,'ColumnName',{'Behavior' 'Feature' 'Timing' 'dRMS'});
set(handles.Table,'ColumnWidth',{150 150 50 50 50});
set(handles.Table,'RowStriping','on','BackgroundColor',[1 1 1; 0.95 0.95 0.95]);

handles.table_data=tmp2;
handles.table='timeseries';
guidata(hObject,handles);

set(handles.Status,'string','Ready.','foregroundcolor','g');
set(handles.figure1,'pointer','arrow');
drawnow;


% ---
function ret_val=print_stats(data,stat)

switch(stat)
  case(0)
    ret_val=num2str(data,3);
  case(1)
    ret_val=num2str(mean(data),3);
  case(2)
    ret_val=[num2str(mean(data),3) ' (' num2str(std(data),3) ')'];
  case(3)
    ret_val=[num2str(mean(data),3) ' (' num2str(std(data)./sqrt(length(data)),3) ')'];
  case(4)
    ret_val=num2str(median(data),3);
  case(5)
    ret_val=[num2str(median(data),3) ' (' num2str(prctile(data,5),3) '-' num2str(prctile(data,95),3) ')'];
  case(6)
    ret_val=[num2str(median(data),3) ' (' num2str(prctile(data,25),3) '-' num2str(prctile(data,75),3) ')'];
end


% ---
function [ct,dp,dn]=calculate_ct_d(data,centraltendency,dispersion)

switch(centraltendency)
  case 1
    ct=nanmean(data);
  case 2
    ct=nanmedian(data);
  case 3
    ct=mode(data);
end
switch(dispersion)
  case 1
    tmp=std(data);
    dp=ct+tmp;
    dn=ct-tmp;
  case 2
    tmp=std(data)./sqrt(length(data));
    dp=ct+tmp;
    dn=ct-tmp;
  case 3
    dp=prctile(data,95);
    dn=prctile(data,5);
  case 4
    dp=prctile(data,75);
    dn=prctile(data,25);
end


% ---
function [h,h2]=errorbarplot(ha,x,b,dn,dp,color)

h = bar(ha,x,b,'grouped');
set(h,'facecolor',color);
fxdata = get(h, 'XData');
fydata = get(h ,'YData');
if(~iscell(fxdata)) xdata{1,1} = fxdata; else xdata = fxdata; end
if(~iscell(fydata)) ydata{1,1} = fydata; else ydata = fydata; end

% Determine number of bars
sizz = size(b);
nb = sizz(1)*sizz(2);
xb = [];
yb = [];
for i = 1:length(xdata),
    xb = [xb xdata{i,1}];
    yb = [yb ydata{i,1}];
end;

%% To find the center of each bar - need to look at the output vectors xb, yb
%% find where yb is non-zero - for each bar there is a pair of non-zero yb values.
%% The center of these values is the middle of the bar
%
%nz = find(yb);
%for i = 1:nb,
%    center(i) = (xb(nz(i*2))-xb(nz((i*2)-1)))/2 + xb(nz((i*2)-1));
%end;

% To place the error bars - use the following:

%errdata=[.1 .2; .3 .4; .5 .6];

hold(ha,'on');
h2=errorbar(ha,xb, b, dn, dp);
set(h2(1),'linewidth',1);            % This changes the thickness of the errorbars
set(h2(1),'color','k');              % This changes the color of the errorbars
set(h2(1),'linestyle','none');       % This removes the connecting


% ---
function ret_val=color_red(arg,flag)

ret_val=num2str(arg);
if(flag)
  ret_val=['<html><font color="red">', ret_val, '</font></html>'];
end


% ---
function tmp=calculate_statistics(table_data,behaviorlist,grouplist,fid,crit)

nrows=13;  if(length(table_data{1})==2)  nrows=7;  end

tmp={'behavior'};
for b=1:length(table_data)

  grouplistcopy=grouplist;
  g=1;
  while g<=length(table_data{b})
    if(isempty(table_data{b}{g}) || (sum(~isnan(table_data{b}{g}))==0))
      table_data{b}(g)=[];
      grouplistcopy(g)=[];
    else
      g=g+1;
    end
  end

  tmp{nrows*b+1,1}=behaviorlist{b};
  for g=1:length(table_data{b})
    if(b==1)  tmp(2:4,g)=repmat({grouplistcopy{g}},3,1);  end
    p=nan;
    if(sum(~isnan(table_data{b}{g}))>0)
      [~,p,~,~]=kstest(table_data{b}{g});
    end
    tmp{nrows*b+2,g}=color_red(length(table_data{b}{g}),0);
    tmp{nrows*b+3,g}=color_red(p,p<crit);
    tmp{nrows*b+4,g}=color_red(nanstd(table_data{b}{g}),0);
  end
  if(b==1)
    tmp{2,length(table_data{b})+1}='(sample size)';
    tmp{3,length(table_data{b})+1}='(K-S normal)';
    tmp{4,length(table_data{b})+1}='(std. dev.)';
  end

  if(length(table_data{b})==2)
  
    if(b==1)  tmp{5,1}='t-test';  end
    [~,p]=ttest2(table_data{b}{1},table_data{b}{2});
    tmp{nrows*b+5,1}=color_red(p,p<crit);

    if(b==1)  tmp{6,1}='Wilcoxen';  end
    p=ranksum(table_data{b}{1},table_data{b}{2});
    tmp{nrows*b+6,1}=color_red(p,p<crit);

  else

    foo=cellfun(@(x,y) y*ones(1,length(x)),table_data{b},num2cell(1:length(table_data{b})),...
        'uniformoutput',false);
    [p,table,stats]=anova1([table_data{b}{:}],[foo{:}],'off');
    c=multcompare(stats,'alpha',crit,'display','off');
    tmp{nrows*b+5,1}=color_red(p,p<crit);
    if(b==1)  tmp{5,1}='ANOVA';  end
    for g2=1:size(c,1)
      if(b==1)
        tmp(6:8,g2)=repmat({strcat(grouplistcopy{c(g2,1)},'-',grouplistcopy{c(g2,2)})},3,1);
      end
      foo=(sign(c(g2,3))==sign(c(g2,5)));
      tmp{nrows*b+6,g2}=color_red(c(g2,3),foo);
      tmp{nrows*b+7,g2}=color_red(c(g2,4),foo);
      tmp{nrows*b+8,g2}=color_red(c(g2,5),foo);
    end
    if(b==1)
      tmp(6:8,size(c,1)+1)={'(5%, Tukey post-hoc)'; '(mean)'; '(95%)'};
    end

    foo=cellfun(@(x,y) y*ones(1,length(x)),table_data{b},num2cell(1:length(table_data{b})),...
        'uniformoutput',false);
    [p,table,stats]=kruskalwallis([table_data{b}{:}],[foo{:}],'off');
    c=multcompare(stats,'alpha',crit,'display','off');
    tmp{nrows*b+9,1}=color_red(p,p<crit);
    if(b==1)  tmp{9,1}='Kruskal-Wallis';  end
    for g2=1:size(c,1)
      if(b==1)
        tmp(10:12,g2)=repmat({strcat(grouplistcopy{c(g2,1)},'-',grouplistcopy{c(g2,2)})},3,1);
      end
      foo=(sign(c(g2,3))==sign(c(g2,5)));
      tmp{nrows*b+10,g2}=color_red(c(g2,3),foo);
      tmp{nrows*b+11,g2}=color_red(c(g2,4),foo);
      tmp{nrows*b+12,g2}=color_red(c(g2,5),foo);
    end
    if(b==1)
      tmp(10:12,size(c,1)+1)={'(5%, Tukey post-hoc)'; '(mean)'; '(95%)'};
    end
  end
end
tmp{end+1,1}='';

fprintf(fid,'%% statistics\n');
for j=1:(size(tmp,1)/nrows)
  for i=1:nrows
    if(j==1)  fprintf(fid,'%% ');  end
    tmp(nrows*(j-1)+i,:);
    ans(cellfun(@isempty,ans))='';
    cellfun(@(x) regexprep(x,'<[^<]*>','*'),ans,'uniformoutput',false);
    fprintf(fid,'%s, ',ans{:});
    fprintf(fid,'\n');
  end
end


% ---
function handles=behavior_barchart_plot(handles)

set(handles.Status,'string','Thinking...','foregroundcolor','b');
set(handles.figure1,'pointer','watch');
drawnow;

handlesexperimentlist=[handles.experimentlist{:}];

cumsum_num_indi_per_exp=[0 cumsum(handles.individuals_feature)];
cumsum_num_exp_per_group=[0 cumsum(cellfun(@length,handles.experimentlist))];
mat2cell(cumsum_num_exp_per_group(1:end-1),1,ones(1,length(cumsum_num_exp_per_group)-1));
cellfun(@(x,y) x+y,handles.experimentvalue,ans,'uniformoutput',false);
selected_exp=[ans{:}];
cumsum_num_selexp_per_group=[0 cumsum(cellfun(@length,handles.experimentvalue))];

%ggee=1:length(handlesexperimentlist);
ggee=selected_exp;
individual=handles.individualidx;
if isnumeric(individual)
  ggee=find(cumsum_num_indi_per_exp<individual,1,'last');
  individual=individual-cumsum_num_indi_per_exp(ggee);
end

fid=fopen('most_recent_figure.csv','w');

handles.type='behavior bar chart';

hf=figure('toolbar','figure');

bb=handles.behaviorvalue;
if(bb==(length(handles.behaviorlist)+1))  bb=1:(bb-1);  end

behavior_logic=handles.behaviorlogic;
score_file2=handles.scorefiles{handles.behaviorvalue2};
score_file3=[];
if(handles.behaviorvalue3>1)
  score_file3=handles.scorefiles{handles.behaviorvalue3-1};
end
sexdata=handles.sexdata;
%perwhat=handles.behaviorbarchart_style;
behaviornot=handles.behaviornot;

h={};
table_data={};
for b=bb

  tstr='';
  ystr='';  if(behaviornot)  ystr='NOT ';  end
  ystr=[ystr char(strrep(handles.behaviorlist(b),'_','-'))];
  switch(handles.behaviorlogic)
    case 2
      ystr=[ystr ' AND '];
    case 3
      ystr=[ystr ' AND NOT '];
    case 4
      ystr=[ystr ' OR '];
  end
  if(handles.behaviorlogic>1)
    ystr=[ystr char(strrep(handles.behaviorlist(handles.behaviorvalue2),'_','-'))];
  end
  ystr=[ystr ' (%)'];
  xstr='group';

  print_csv_help(fid,handles.type,tstr,xstr,ystr);

  if(length(bb)>1)
    ceil(sqrt(length(bb)));
    ha=subplot(ceil(length(bb)/ans),ans,b,'parent',hf);
  else
    ha=subplot(1,1,1,'parent',hf);
  end
  hold(ha,'on');

  score_file=handles.scorefiles{b};

  num_indi=0;
  collated_data=cell(1,length(ggee));
  parfor gei=1:numel(ggee)
  %for gei=1:numel(ggee)
    ge = ggee(gei);

    %if(ischar(individual)&&(~ismember(ge,selected_exp)))  continue;  end

    behavior_data=load(fullfile(handlesexperimentlist{ge},score_file));
    behavior_data2=[];
    if(behavior_logic>1)
      behavior_data2=load(fullfile(handlesexperimentlist{ge},score_file2));
    end
    behavior_data3=[];
    if(handles.behaviorvalue3>1)
      behavior_data3=load(fullfile(handlesexperimentlist{ge},score_file3));
    end

    [behavior_data,behavior_data2,behavior_data3,~]=...
        cull_short_trajectories(handles,behavior_data,behavior_data2,behavior_data3,[]);
    num_indi=num_indi+length(behavior_data.allScores.scores);

    traj_len=behavior_data.allScores.tEnd-behavior_data.allScores.tStart;

    frames_labelled=nan(1,length(behavior_data.allScores.t0s));
    frames_total=nan(1,length(behavior_data.allScores.t0s));
    sex=nan(1,length(behavior_data.allScores.t0s));

    for i=1:length(behavior_data.allScores.t0s)  % individual
      tmp1=zeros(1,behavior_data.allScores.tEnd(i)-behavior_data.allScores.tStart(i)+1);
      tmp1(behavior_data.allScores.t0s{i}-behavior_data.allScores.tStart(i)+1)=1;
      tmp1(behavior_data.allScores.t1s{i}-behavior_data.allScores.tStart(i)+1)=-1;
      tmp1=logical(cumsum(tmp1));

      tmp2=[];
      if(behavior_logic>1)
        tmp2=zeros(1,behavior_data2.allScores.tEnd(i)-behavior_data2.allScores.tStart(i)+1);
        tmp2(behavior_data2.allScores.t0s{i}-behavior_data.allScores.tStart(i)+1)=1;
        tmp2(behavior_data2.allScores.t1s{i}-behavior_data.allScores.tStart(i)+1)=-1;
        tmp2=logical(cumsum(tmp2));
      end

      tmp3=[];
      if(handles.behaviorvalue3>1)
        tmp3=zeros(1,behavior_data3.allScores.tEnd(i)-behavior_data3.allScores.tStart(i)+1);
        tmp3(behavior_data3.allScores.t0s{i}-behavior_data.allScores.tStart(i)+1)=1;
        tmp3(behavior_data3.allScores.t1s{i}-behavior_data.allScores.tStart(i)+1)=-1;
        tmp3=logical(cumsum(tmp3));
        if(handles.behaviornormalizenot)
          tmp3=~tmp3;
        end
      end

      if(behaviornot)  tmp1=~tmp1;  end

      partition_idx=[];
      switch(behavior_logic)
        case(1)
          partition_idx=tmp1;
        case(2)
          partition_idx=tmp1 & tmp2;
        case(3)
          partition_idx=tmp1 & ~tmp2;
        case(4)
          partition_idx=tmp1 | tmp2;
        case(5)
          partition_idx=tmp1 | ~tmp2;
      end

      % KB: for some reason partition_idx was of size > trajectory length
      if numel(partition_idx) > numel(sexdata{ge}{i}),
        warning('More frames of behaviors detected than frames in trajectory by %d',numel(partition_idx)-numel(sexdata{ge}{i}));
        partition_idx = partition_idx(1:numel(sexdata{ge}{i}));
      end

      
      sex(i)=sum(sexdata{ge}{i}(1:length(partition_idx))) > (length(partition_idx)/2);
      frames_labelled(i)=sum(partition_idx);
      if(handles.behaviorvalue3==1)
        frames_total(i)=length(partition_idx);
      else
        frames_total(i)=sum(tmp3);
      end
    end

    collated_data{gei}={frames_labelled frames_total sex traj_len};
  end

  %idx=cellfun(@isempty,collated_data);
  %collated_data=collated_data(~idx);

  if(num_indi==0)
    delete(hf);
    set(handles.Status,'string','Ready.','foregroundcolor','g');
    set(handles.figure1,'pointer','arrow');
    uiwait(errordlg('no valid data.  check minimum trajectory length.'));  drawnow;
    return;
  end

  switch(individual)
    case 'A'
      frames_labelled=cellfun(@(x) x{1},collated_data,'uniformoutput',false);
      frames_total=cellfun(@(x) x{2},collated_data,'uniformoutput',false);
      traj_len=cellfun(@(x) x{4},collated_data,'uniformoutput',false);
    case {'M'}
      frames_labelled=cellfun(@(x) x{1}(x{3}==1),collated_data,'uniformoutput',false);
      frames_total=cellfun(@(x) x{2}(x{3}==1),collated_data,'uniformoutput',false);
      traj_len=cellfun(@(x) x{4}(x{3}==1),collated_data,'uniformoutput',false);
    case {'F'}
      frames_labelled=cellfun(@(x) x{1}(x{3}==0),collated_data,'uniformoutput',false);
      frames_total=cellfun(@(x) x{2}(x{3}==0),collated_data,'uniformoutput',false);
      traj_len=cellfun(@(x) x{4}(x{3}==0),collated_data,'uniformoutput',false);
    otherwise
      frames_labelled=cellfun(@(x) x{1}(individual),collated_data,'uniformoutput',false);
      frames_total=cellfun(@(x) x{2}(individual),collated_data,'uniformoutput',false);
      traj_len=cellfun(@(x) x{4}(individual),collated_data,'uniformoutput',false);
  end

  exp_separators=[];  maxy=0;  k=[];  m=0;  table_data{end+1}=[];
  for g=1:length(handles.grouplist)
    color=handles.colors(g,:);

    if ischar(individual)
      idx=(cumsum_num_selexp_per_group(g)+1):(cumsum_num_selexp_per_group(g+1));
    else
      find(cumsum_num_exp_per_group<ggee,1,'last');
      if(ans~=g)  continue;  end
      idx=1;
    end

    xticklabels{g}=handles.grouplist{g};

    switch(handles.behaviorbarchart_style)
      case 1  % per group
        table_data{end}(g)=100*sum([frames_labelled{idx}])./sum([frames_total{idx}]);
        h{g}=bar(ha,g,table_data{end}(g));
        set(h{g},'facecolor',color);
      case 2  % per experiment, error bars
        table_data{end}{g}=100*cellfun(@sum,frames_labelled(idx))./cellfun(@sum,frames_total(idx));
        [ct(g),dp(g),dn(g)]=...
            calculate_ct_d(table_data{end}{g},handles.centraltendency,handles.dispersion);
        h{g}=errorbarplot(ha,g,ct(g),ct(g)-dn(g),dp(g)-ct(g),color);
      case 3  % per experiment, box
        table_data{end}{g}=100*cellfun(@sum,frames_labelled(idx))./cellfun(@sum,frames_total(idx));
        h{g}=boxplot(ha,table_data{end}{g},'positions',g,'widths',0.5,'colors',color);
      case 4  % per fly, grouped
        fprintf(fid,['%% data, %s\n'],xticklabels{g});
        cumsum(cellfun(@length,frames_labelled(idx)))';
        exp_separators=[exp_separators; ans+sum(k)];
        table_data{end}{g}=100.*[frames_labelled{idx}]./[frames_total{idx}];
        maxy=max([maxy table_data{end}{g}]);
        h{g}=bar(ha,(1:length(table_data{end}{g}))+sum(k),table_data{end}{g},...
            'barwidth',1,'edgecolor','none');
        set(h{g},'facecolor',color);
        k(end+1)=length(table_data{end}{g});
        fprintf(fid,'%g, ',[table_data{end}{g}]);
        fprintf(fid,'\n');
      case 5  % per fly, stern-style
        fprintf(fid,['%% data, %s\n'],xticklabels{g});
        table_data{end}{g}=cell(1,length(frames_labelled(idx)));
        for e=idx
          table_data{end}{g}{e}=100.*frames_labelled{e}./frames_total{e};
          [ct,dp,dn]=calculate_ct_d(table_data{end}{g}{e},...
              handles.centraltendency,handles.dispersion);
          h{g}=plot(ha,m,ct,'o','color',color);
          plot(ha,[m m],[dp dn],'-','color',color);
          plot(ha,m+(1:length(table_data{end}{g}{e})),table_data{end}{g}{e},'.','color',color);
          m=m+16+length(table_data{end}{g}{e});
        end
        [ct,dp,dn]=calculate_ct_d([table_data{end}{g}{:}],...
            handles.centraltendency,handles.dispersion);
        plot(ha,m,ct,'o','color',color,'markersize',9);
        plot(ha,[m m],[dp dn],'-','color',color,'linewidth',3);
        m=m+24;
        k(end+1)=24+16*length(table_data{end}{g})+length([table_data{end}{g}{:}]);
        fprintf(fid,'%g, ',[table_data{end}{g}{:}]);
        fprintf(fid,'\n');
      case 6  % per fly, trajectory length
        fprintf(fid,['%% data, %s\n'],xticklabels{g});
        cumsum(cellfun(@length,traj_len(idx)))';
        exp_separators=[exp_separators; ans+sum(k)];
        table_data{end}{g}=[traj_len{idx}];
        maxy=max([maxy table_data{end}{g}]);
        h{g}=bar(ha,(1:length(table_data{end}{g}))+sum(k),table_data{end}{g},...
            'barwidth',1,'edgecolor','none');
        set(h{g},'facecolor',color);
        k(end+1)=length(table_data{end}{g});
        fprintf(fid,'%g, ',[table_data{end}{g}]);
        fprintf(fid,'\n');
    end
  end

  switch(handles.behaviorbarchart_style)
    case 1  % per group
      fprintf(fid,['%% xdata\n']);  fprintf(fid,'%s, ',xticklabels{:});  fprintf(fid,'\n');
      fprintf(fid,['%% ydata, per group\n']);  fprintf(fid,'%g, ',table_data{end});  fprintf(fid,'\n');
    case 2  % per experiment, error bars
      fprintf(fid,['%% xdata\n']);  fprintf(fid,'%s, ',xticklabels{:});  fprintf(fid,'\n');
      fprintf(fid,['%% ydata, CT+D\n']);  fprintf(fid,'%g, ',dp);  fprintf(fid,'\n');
      fprintf(fid,['%% ydata, CT-D\n']);  fprintf(fid,'%g, ',dn);  fprintf(fid,'\n');
      fprintf(fid,['%% ydata, CT\n']);  fprintf(fid,'%g, ',ct);  fprintf(fid,'\n');
    case {4,6}  % per fly, grouped
      l=exp_separators(1:2:(end-1));
      r=exp_separators(2:2:end);
      hh=patch(0.5+[l r r l l]',repmat([0 0 maxy*1.05 maxy*1.05 0]',1,floor(length(exp_separators)/2)),...
          [0.95 0.95 0.95],'parent',ha);
      set(hh,'edgecolor','none');
      set(ha,'children',flipud(get(ha,'children')));
      k=round(cumsum(k)-k/2);
    case 5  % per fly, stern-style
      k=round(cumsum(k)-k/2);
  end

  fprintf(fid,'\n%% raw data\n');
  for g=1:length(handles.grouplist)
    if ischar(individual)
      idx=(cumsum_num_selexp_per_group(g)+1):(cumsum_num_selexp_per_group(g+1));
    else
      find(cumsum_num_exp_per_group<ggee,1,'last');
      if(ans~=g)  continue;  end
      idx=1;
    end
    fprintf(fid,['%% group ' handles.grouplist{g} '\n']);
    for e=1:length(idx)
      fprintf(fid,'%% experiment %s\n',handlesexperimentlist{selected_exp(idx(e))});
      fprintf(fid,'%% frames labelled\n');
      fprintf(fid,'%d, ',frames_labelled{idx(e)});
      fprintf(fid,'\n');
      fprintf(fid,'%% frames total\n');
      fprintf(fid,'%d, ',frames_total{idx(e)});
      fprintf(fid,'\n');
      fprintf(fid,'\n');
    end
  end

  if(isempty(k))  k=1:length(frames_labelled);  end
  ylabel(ha,ystr,'interpreter','none');
  set(ha,'xtick',k,'xticklabel',xticklabels);
  axis(ha,'tight');  vt=axis;
  axisalmosttight;  vat=axis;
  if(handles.behaviorbarchart_style==4)
    axis(ha,[vat(1) vat(2) 0 vt(4)]);
  else
    axis(ha,[vat(1) vat(2) 0 vat(4)]);
  end
  %fprintf(fid,'\n');

end
%idx=find(~cellfun(@isempty,h));
%if ischar(individual)
%  legend(cellfun(@(x) x(1),h(idx)),handles.grouplist);
%else
%  legend(cellfun(@(x) x(1),h(idx)),handles.individuallist(handles.individualvalue));
%end

%if((ismember(handles.behaviorbarchart_style,[2 3])) && (individual<4))
uicontrol(hf,'style','pushbutton','string','Params','position',[5 5 60 20],...
    'callback',@figure_params_callback);
if((ismember(handles.behaviorbarchart_style,[2 3 4])) && ischar(individual) && (length(handles.grouplist)>1))
  uicontrol(hf,'style','pushbutton','string','Stats','position',[70 5 50 20],...
      'callback',@figure_stats_callback);
  handles.statistics=calculate_statistics(table_data,handles.behaviorlist(bb),handles.grouplist,...
      fid,handles.pvalue);
%  set(handles.Table,'Data',tmp);
%  set(handles.Table,'ColumnWidth','auto');
%  set(handles.Table,'ColumnName',{''});
%else
%  set(handles.Table,'Data',[]);
%  set(handles.Table,'ColumnName',{});
end

fclose(fid);

guidata(hf,handles);

set(handles.Status,'string','Ready.','foregroundcolor','g');
set(handles.figure1,'pointer','arrow');
drawnow;


% ---
function handles=behavior_timeseries_plot(handles)

set(handles.Status,'string','Thinking...','foregroundcolor','b');
set(handles.figure1,'pointer','watch');
drawnow;

handlesexperimentlist=[handles.experimentlist{:}];

cumsum_num_indi_per_exp=[0 cumsum(handles.individuals_feature)];
cumsum_num_exp_per_group=[0 cumsum(cellfun(@length,handles.experimentlist))];
mat2cell(cumsum_num_exp_per_group(1:end-1),1,ones(1,length(cumsum_num_exp_per_group)-1));
cellfun(@(x,y) x+y,handles.experimentvalue,ans,'uniformoutput',false);
selected_exp=[ans{:}];
cumsum_num_selexp_per_group=[0 cumsum(cellfun(@length,handles.experimentvalue))];

%ggee=1:length(handlesexperimentlist);
ggee=selected_exp;
individual=handles.individualidx;
if isnumeric(individual)
  ggee=find(cumsum_num_indi_per_exp<individual,1,'last');
  individual=individual-cumsum_num_indi_per_exp(ggee);
end

fid=fopen('most_recent_figure.csv','w');

handles.type='behavior time series';

hf=figure('toolbar','figure');

bb=handles.behaviorvalue;
if(bb==(length(handles.behaviorlist)+1))  bb=1:(bb-1);  end

behavior_logic=handles.behaviorlogic;
score_file2=handles.scorefiles{handles.behaviorvalue2};
score_file3=[];
if(handles.behaviorvalue3>1)
  score_file3=handles.scorefiles{handles.behaviorvalue3-1};
end
sexdata=handles.sexdata;
convolutionwidth=round(handles.convolutionwidth*handles.fps);
style=handles.behaviortimeseries_style;
centraltendency=handles.centraltendency;
dispersion=handles.dispersion;
xoffset=handles.xoffset;
behaviornot=handles.behaviornot;

h=[];
for b=bb

  if(length(bb)>1)
    ceil(sqrt(length(bb)));
    ha=subplot(ceil(length(bb)/ans),ans,b,'parent',hf);
  else
    ha=subplot(1,1,1,'parent',hf);
  end
  hold(ha,'on');

  score_file=handles.scorefiles{b};

  num_indi=0;
  raw_data=cell(length(ggee),1);
  behavior_cumulative=cell(length(ggee),1);
  parfor gei=1:numel(ggee),
  %for gei=1:numel(ggee),
    ge = ggee(gei);
    
    %if(ischar(individual)&&(~ismember(ge,selected_exp)))  continue;  end

    behavior_data=load(fullfile(handlesexperimentlist{ge},score_file));
    behavior_data2=[];
    if(behavior_logic>1)
      behavior_data2=load(fullfile(handlesexperimentlist{ge},score_file2));
    end
    behavior_data3=[];
    if(handles.behaviorvalue3>1)
      behavior_data3=load(fullfile(handlesexperimentlist{ge},score_file3));
    end

    [behavior_data,behavior_data2,behavior_data3,~]=...
        cull_short_trajectories(handles,behavior_data,behavior_data2,behavior_data3,[]);
    num_indi=num_indi+length(behavior_data.allScores.scores);

    if(xoffset==1)
      parfor_tmp=zeros(2,max(behavior_data.allScores.tEnd));
    else
      parfor_tmp=zeros(2,max(behavior_data.allScores.tEnd)-min(behavior_data.allScores.tStart));
    end

    for i=1:length(behavior_data.allScores.t0s)   % individual
      if(isnumeric(individual)&&(individual~=i))  continue;  end

      tmp1=zeros(1,behavior_data.allScores.tEnd(i));
      tmp1(behavior_data.allScores.t0s{i})=1;
      tmp1(behavior_data.allScores.t1s{i})=-1;
      tmp1=logical(cumsum(tmp1));

      tmp2=[];
      if(behavior_logic>1)
        tmp2=zeros(1,behavior_data2.allScores.tEnd(i));
        tmp2(behavior_data2.allScores.t0s{i})=1;
        tmp2(behavior_data2.allScores.t1s{i})=-1;
        tmp2=logical(cumsum(tmp2));
      end

      tmp3=[];
      if(handles.behaviorvalue3>1)
        tmp3=zeros(1,behavior_data3.allScores.tEnd(i));
        tmp3(behavior_data3.allScores.t0s{i})=1;
        tmp3(behavior_data3.allScores.t1s{i})=-1;
        tmp3=logical(cumsum(tmp3));
        if(handles.behaviornormalizenot)
          tmp3=~tmp3;
        end
      end

      if(behaviornot)  tmp1=~tmp1;  end

      partition_idx=zeros(2,length(tmp1));
      switch(behavior_logic)
        case(1)
          partition_idx(1,:)=tmp1;
        case(2)
          partition_idx(1,:)=tmp1 & tmp2;
        case(3)
          partition_idx(1,:)=tmp1 & ~tmp2;
        case(4)
          partition_idx(1,:)=tmp1 | tmp2;
        case(5)
          partition_idx(1,:)=tmp1 | ~tmp2;
      end
      if(handles.behaviorvalue3==1)
        partition_idx(2,behavior_data.allScores.tStart(i):behavior_data.allScores.tEnd(i))=1;
      else
        partition_idx(2,:)=tmp3;
      end

      [ones(1,behavior_data.allScores.tEnd(i))];
      switch(individual)
        case('M')
          [ones(1,behavior_data.allScores.tStart(i)-1) sexdata{ge}{i}];
        case('F')
          [ones(1,behavior_data.allScores.tStart(i)-1) ~sexdata{ge}{i}];
      end
      partition_idx(1,:) = partition_idx(1,:) & ans;
      partition_idx(2,:) = partition_idx(2,:) & ans;
      
      idx1 = find(partition_idx(1,:));
      idx2 = find(partition_idx(2,:));
      switch(xoffset)
        case(1)
        case(2)
          idx1 = idx1 - behavior_data.allScores.tStart(i) + 1;
          idx2 = idx2 - behavior_data.allScores.tStart(i) + 1;
        case(3)
          idx1 = idx1 - min(behavior_data.allScores.tStart) + 1;
          idx2 = idx2 - min(behavior_data.allScores.tStart) + 1;
      end
      parfor_tmp(1,idx1)=parfor_tmp(1,idx1)+1;
      parfor_tmp(2,idx2)=parfor_tmp(2,idx2)+1;
    end
    raw_data{gei}=parfor_tmp(1,:)./parfor_tmp(2,:);
    find(parfor_tmp(2,:)==0);  parfor_tmp(2,ans)=nan;
    behavior_cumulative{gei}=conv(parfor_tmp(1,:),ones(1,convolutionwidth),'same') ./ ...
       conv(parfor_tmp(2,:),ones(1,convolutionwidth),'same');
  end

  if(num_indi==0)
    delete(hf);
    set(handles.Status,'string','Ready.','foregroundcolor','g');
    set(handles.figure1,'pointer','arrow');
    uiwait(errordlg('no valid data.  check minimum trajectory length.'));  drawnow;
    return;
  end

  max(cellfun(@(x) size(x,2),behavior_cumulative));
  cellfun(@(x) [x nan(size(x,1),ans-size(x,2))],behavior_cumulative,'uniformoutput',false);
  behavior_cumulative=cat(1,ans{:});

  tstr='';
  time_base=(1:size(behavior_cumulative,2))./handles.fps;
  xstr='time (sec)';
  if(time_base(end)>300)
    time_base=time_base./60;
    xstr='time (min)';
  end
  if(time_base(end)>300)
    time_base=time_base./60;
    xstr='time (hr)';
  end
  if(time_base(end)>120)
    time_base=time_base./24;
    xstr='time (d)';
  end
  ystr='';  if(behaviornot)  ystr='NOT ';  end
  ystr=[ystr char(strrep(handles.behaviorlist(b),'_','-'))];
  switch(handles.behaviorlogic)
    case 2
      ystr=[ystr ' AND '];
    case 3
      ystr=[ystr ' AND NOT '];
    case 4
      ystr=[ystr ' OR '];
  end
  if(handles.behaviorlogic>1)
    ystr=[ystr char(strrep(handles.behaviorlist(handles.behaviorvalue2),'_','-'))];
  end
  ystr=[ystr ' (%)'];

  print_csv_help(fid,handles.type,tstr,xstr,ystr);

  for g=1:length(handles.grouplist)
    color=handles.colors(g,:);

    fprintf(fid,['%% group ' handles.grouplist{g} '\n']);

    if ischar(individual)
      idx=(cumsum_num_selexp_per_group(g)+1):(cumsum_num_selexp_per_group(g+1));
    else
      find(cumsum_num_exp_per_group<ggee,1,'last');
      if(ans~=g)  continue;  end
      idx=1;
    end

    h(g)=plot_it(ha,time_base,100.*behavior_cumulative(idx,:),...
        style,centraltendency,dispersion,color,1,fid,handlesexperimentlist(idx));
  end

  fprintf(fid,'\n%% raw data\n');
  for g=1:length(handles.grouplist)
    if ischar(individual)
      idx=(cumsum_num_selexp_per_group(g)+1):(cumsum_num_selexp_per_group(g+1));
    else
      find(cumsum_num_exp_per_group<ggee,1,'last');
      if(ans~=g)  continue;  end
      idx=1;
    end
    fprintf(fid,['%% group ' handles.grouplist{g} '\n']);
    for e=1:length(idx)
      fprintf(fid,'%% experiment %s\n',handlesexperimentlist{selected_exp(idx(e))});
      print_csv_data(fid,raw_data{idx(e)});
%      fprintf(fid,'%g, ',raw_data{idx(e)});
%      fprintf(fid,'\n');
      fprintf(fid,'\n');
    end
  end

  xlabel(ha,xstr,'interpreter','none');
  ylabel(ha,ystr,'interpreter','none');
  axis(ha,'tight');  zoom(ha,'reset');

end
idx=find(h>0);
if ischar(individual)
  legend(ha,h(idx),[cellfun(@(x) [x ' ' handles.individuallist{handles.individualvalue}],...
      handles.grouplist,'uniformoutput',false)],'interpreter','none');
else
  legend(ha,h(idx),handles.individuallist(handles.individualvalue),'interpreter','none');
end

uicontrol(hf,'style','pushbutton','string','Params','position',[5 5 60 20],...
    'callback',@figure_params_callback);

fclose(fid);

guidata(hf,handles);

set(handles.Status,'string','Ready.','foregroundcolor','g');
set(handles.figure1,'pointer','arrow');
drawnow;


% --- 
function [bout_lengths sex inter_bout_lengths inter_sex]=...
    calculate_boutstats(behavior_data,behavior_logic,behavior_data2,sexdata,behaviornot)

bout_lengths=cell(1,length(behavior_data.allScores.t0s));
inter_bout_lengths=cell(1,length(behavior_data.allScores.t0s));
sex=cell(1,length(behavior_data.allScores.t0s));
inter_sex=cell(1,length(behavior_data.allScores.t0s));
for i=1:length(behavior_data.allScores.t0s)  % individual
  tmp1=zeros(1,behavior_data.allScores.tEnd(i)-behavior_data.allScores.tStart(i)+1);
  tmp1(behavior_data.allScores.t0s{i}-behavior_data.allScores.tStart(i)+1)=1;
  tmp1(behavior_data.allScores.t1s{i}-behavior_data.allScores.tStart(i)+1)=-1;
  tmp1=logical(cumsum(tmp1));

  tmp2=[];
  if(behavior_logic>1)
    tmp2=zeros(1,behavior_data2.allScores.tEnd(i)-behavior_data2.allScores.tStart(i)+1);
    tmp2(behavior_data2.allScores.t0s{i}-behavior_data.allScores.tStart(i)+1)=1;
    tmp2(behavior_data2.allScores.t1s{i}-behavior_data.allScores.tStart(i)+1)=-1;
    tmp2=logical(cumsum(tmp2));
  end

  if(behaviornot)  tmp1=~tmp1;  end

  partition_idx=[];
  switch(behavior_logic)
    case(1)
      partition_idx=tmp1;
    case(2)
      partition_idx=tmp1 & tmp2;
    case(3)
      partition_idx=tmp1 & ~tmp2;
    case(4)
      partition_idx=tmp1 | tmp2;
    case(5)
      partition_idx=tmp1 | ~tmp2;
  end

  % inter-bout NOT <behavior> is not quite the same as bout <behavior>

  partition_idx=[0 partition_idx 0];
  start=1+find(~partition_idx(1:(end-1)) &  partition_idx(2:end))-1;
  stop =  find( partition_idx(1:(end-1)) & ~partition_idx(2:end))-1;
  if(length(start)>0)
    bout_lengths{i}=stop-start+1;
    sex{i}=zeros(1,length(bout_lengths{i}));
    for j=1:length(bout_lengths{i})
      sex{i}(j)=sum(sexdata{i}(start(j):stop(j))) > (bout_lengths{i}(j)/2);
    end
    if(length(start)>1)
      inter_bout_lengths{i}=start(2:end)-stop(1:(end-1));
      inter_sex{i}=zeros(1,length(inter_bout_lengths{i}));
      for j=1:length(inter_bout_lengths{i})
        inter_sex{i}(j)=sum(sexdata{i}(stop(j):start(j+1))) > (inter_bout_lengths{i}(j)/2);
      end
    end
  end
end


% ---
function handles=bout_stats_plot(handles)

set(handles.Status,'string','Thinking...','foregroundcolor','b');
set(handles.figure1,'pointer','watch');
drawnow;

handlesexperimentlist=[handles.experimentlist{:}];

cumsum_num_indi_per_exp=[0 cumsum(handles.individuals_feature)];
cumsum_num_exp_per_group=[0 cumsum(cellfun(@length,handles.experimentlist))];
mat2cell(cumsum_num_exp_per_group(1:end-1),1,ones(1,length(cumsum_num_exp_per_group)-1));
cellfun(@(x,y) x+y,handles.experimentvalue,ans,'uniformoutput',false);
selected_exp=[ans{:}];
cumsum_num_selexp_per_group=[0 cumsum(cellfun(@length,handles.experimentvalue))];

%ggee=1:length(handlesexperimentlist);
ggee=selected_exp;
individual=handles.individualidx;
if isnumeric(individual)
  ggee=find(cumsum_num_indi_per_exp<individual,1,'last');
  individual=individual-cumsum_num_indi_per_exp(ggee);
end

fid=fopen('most_recent_figure.csv','w');

handles.type='bout stats';

hf=figure('toolbar','figure');

bb=handles.behaviorvalue;
if(bb==(length(handles.behaviorlist)+1))  bb=1:(bb-1);  end

behavior_logic=handles.behaviorlogic;
score_file2=handles.scorefiles{handles.behaviorvalue2};
sexdata=handles.sexdata;
behaviornot=handles.behaviornot;

h=[];
table_data={};
for b=bb

  if(length(bb)>1)
    ceil(sqrt(length(bb)));
    ha=subplot(ceil(length(bb)/ans),ans,b,'parent',hf);
  else
    ha=subplot(1,1,1,'parent',hf);
  end
  hold(ha,'on');

  score_file=handles.scorefiles{b};

  num_indi=0;
  collated_data=cell(1,length(ggee));
  parfor gei=1:length(ggee)
  %for gei=1:length(ggee)
    ge = ggee(gei);

    behavior_data=load(fullfile(handlesexperimentlist{ge},score_file));
    behavior_data2=[];
    if(behavior_logic>1)
      behavior_data2=load(fullfile(handlesexperimentlist{ge},score_file2));
    end

    [behavior_data,behavior_data2,~,~]=cull_short_trajectories(handles,behavior_data,behavior_data2,[],[]);
    num_indi=num_indi+length(behavior_data.allScores.scores);

    [bout_lengths sex inter_bout_lengths inter_sex]=...
        calculate_boutstats(behavior_data,behavior_logic,behavior_data2,sexdata{ge},behaviornot);
    collated_data{gei}={bout_lengths sex inter_bout_lengths inter_sex};
  end

  if(num_indi==0)
    delete(hf);
    set(handles.Status,'string','Ready.','foregroundcolor','g');
    set(handles.figure1,'pointer','arrow');
    uiwait(errordlg('no valid data.  check minimum trajectory length.'));  drawnow;
    return;
  end

  tstr='';  if(behaviornot)  tstr='NOT ';  end
  tstr=[tstr char(strrep(handles.behaviorlist(b),'_','-'))];
  switch(handles.behaviorlogic)
    case 2
      tstr=[tstr ' AND '];
    case 3
      tstr=[tstr ' AND NOT '];
    case 4
      tstr=[tstr ' OR '];
    case 5
      tstr=[tstr ' OR NOT'];
  end
  if(handles.behaviorlogic>1)
    tstr=[tstr char(strrep(handles.behaviorlist(handles.behaviorvalue2),'_','-'))];
  end
  ystr='bout length (sec)';
  if(handles.boutstats_style2==2)  ystr=['inter-' ystr];  end
  xstr='group';

  print_csv_help(fid,handles.type,tstr,xstr,ystr);

  idx=cellfun(@isempty,collated_data);
  collated_data=collated_data(~idx);

  idx=handles.boutstats_style2*2-1;
  switch(individual)
    case 'A'
      length_data=cellfun(@(x) x{idx},collated_data,'uniformoutput',false);
    case {'M'}
      for ge=1:length(collated_data)
        length_data{ge}=cellfun(@(x,y) x(y==1),...
            collated_data{ge}{idx},collated_data{ge}{idx+1},'uniformoutput',false);
      end
    case {'F'}
      for ge=1:length(collated_data)
        length_data{ge}=cellfun(@(x,y) x(y==0),...
            collated_data{ge}{idx},collated_data{ge}{idx+1},'uniformoutput',false);
      end
    otherwise
      length_data=cellfun(@(x) x{idx}(individual),collated_data,'uniformoutput',false);
  end

  exp_separators=[];  maxy=0;  k=[];  m=0;  table_data{end+1}=[];
  for g=1:length(handles.grouplist)
    color=handles.colors(g,:);

    if ischar(individual)
      idx=(cumsum_num_selexp_per_group(g)+1):(cumsum_num_selexp_per_group(g+1));
    else
      find(cumsum_num_exp_per_group<ggee,1,'last');
      if(ans~=g)  continue;  end
      idx=1;
    end

    xticklabels{g}=handles.grouplist{g};

    switch(handles.boutstats_style)
      case 1  % per experiment, error bars
        table_data{end}{g}=cellfun(@(x) nanmean([x{:}]./handles.fps),length_data(idx));
        [ct(g),dp(g),dn(g)]=...
            calculate_ct_d(table_data{end}{g},handles.centraltendency,handles.dispersion);
        h(g)=errorbarplot(ha,g,ct(g),ct(g)-dn(g),dp(g)-ct(g),color);
      case 2  % per fly, grouped
        cumsum(cellfun(@length,length_data(idx)))';
        exp_separators=[exp_separators; ans+sum(k)];
        table_data{end}{g}=cellfun(@nanmean,[length_data{idx}])./handles.fps;
        maxy=max([maxy table_data{end}{g}]);
        h(g)=bar(ha,(1:length(table_data{end}{g}))+sum(k),table_data{end}{g},...
            'barwidth',1,'edgecolor','none');
        set(h(g),'facecolor',color);
        k(end+1)=length(table_data{end}{g});
        fprintf(fid,['%% data, %s\n'],xticklabels{g});
        fprintf(fid,'%g, ',[table_data{end}{g}]);
        fprintf(fid,'\n');
    end
  end

  switch(handles.boutstats_style)
    case 1  % per experiment, error bars
      fprintf(fid,['%% xdata\n']);  fprintf(fid,'%s, ',xticklabels{:});  fprintf(fid,'\n');
      fprintf(fid,['%% ydata, CT+D\n']);  fprintf(fid,'%g, ',dp);  fprintf(fid,'\n');
      fprintf(fid,['%% ydata, CT-D\n']);  fprintf(fid,'%g, ',dn);  fprintf(fid,'\n');
      fprintf(fid,['%% ydata, CT\n']);  fprintf(fid,'%g, ',ct);  fprintf(fid,'\n');
    case 2  % per fly, grouped
      l=exp_separators(1:2:(end-1));
      r=exp_separators(2:2:end);
      hh=patch(0.5+[l r r l l]',repmat([0 0 maxy*1.05 maxy*1.05 0]',1,floor(length(exp_separators)/2)),...
          [0.95 0.95 0.95],'parent',ha);
      set(hh,'edgecolor','none');
      set(ha,'children',flipud(get(ha,'children')));
      k=round(cumsum(k)-k/2);
  end

  fprintf(fid,'\n%% raw data\n');
  for g=1:length(handles.grouplist)
    if ischar(individual)
      idx=(cumsum_num_selexp_per_group(g)+1):(cumsum_num_selexp_per_group(g+1));
    else
      find(cumsum_num_exp_per_group<ggee,1,'last');
      if(ans~=g)  continue;  end
      idx=1;
    end
    fprintf(fid,['%% group ' handles.grouplist{g} '\n']);
    for e=1:length(idx)
      fprintf(fid,'%% experiment %s\n',handlesexperimentlist{selected_exp(idx(e))});
      for i=1:length(length_data{idx(e)})
        fprintf(fid,'%% individual %d\n',i);
        print_csv_data(fid,length_data{idx(e)}{i}./handles.fps);
        fprintf(fid,'\n');
      end
    end
  end

  if(isempty(k))  k=1:length(length_data);  end
  title(ha,tstr,'interpreter','none');
  ylabel(ha,ystr,'interpreter','none');
  set(ha,'xtick',k,'xticklabel',xticklabels);
  axis(ha,'tight');  vt=axis;
  axisalmosttight;  vat=axis;
  if(handles.boutstats_style==2)
    axis(ha,[vat(1) vat(2) 0 vt(4)]);
  else
    axis(ha,[vat(1) vat(2) 0 vat(4)]);
  end
  fprintf(fid,'\n');

end

uicontrol(hf,'style','pushbutton','string','Params','position',[5 5 60 20],...
    'callback',@figure_params_callback);
if(ischar(individual) && (length(handles.grouplist)>1))
  uicontrol(hf,'style','pushbutton','string','Stats','position',[70 5 50 20],...
      'callback',@figure_stats_callback);
  handles.statistics=calculate_statistics(table_data,handles.behaviorlist(bb),handles.grouplist,...
      fid,handles.pvalue);
end

fclose(fid);

guidata(hf,handles);

set(handles.Status,'string','Ready.','foregroundcolor','g');
set(handles.figure1,'pointer','arrow');
drawnow;


% --- 
function [during not_during raw_during raw_not_during]=...
    calculate_social(behavior_data,behavior_logic,behavior_data2,feature_data)

if(iscell(feature_data.data{1}))
  vals=unique([feature_data.data{:}]);
  if(length(vals)>2)  error('uhoh');  end
  feature_data.data=cellfun(@(x) strcmp(x,vals{1}),feature_data.data,'uniformoutput',false);
end

during={};  not_during={};
raw_during={};  raw_not_during={};
for i=1:length(behavior_data.allScores.t0s)  % individual
  tmp1=false(1,length(feature_data.data{i}));
  for j=1:length(behavior_data.allScores.t0s{i})
    tmp1((behavior_data.allScores.t0s{i}(j):(behavior_data.allScores.t1s{i}(j)-1))...
        -behavior_data.allScores.tStart(i)+1)=true;
  end
  if(behavior_logic>1)
    tmp2=false(1,length(feature_data.data{i}));
    for j=1:length(behavior_data2.allScores.t0s{i})
      tmp2((behavior_data2.allScores.t0s{i}(j):(behavior_data2.allScores.t1s{i}(j)-1))...
          -behavior_data2.allScores.tStart(i)+1)=true;
    end
  end
  switch(behavior_logic)
    case(1)
      partition_idx=tmp1;
    case(2)
      partition_idx=tmp1 & tmp2;
    case(3)
      partition_idx=tmp1 & ~tmp2;
    case(4)
      partition_idx=tmp1 | tmp2;
  end
  if length(partition_idx)==(length(feature_data.data{i}+1))  % n-1 features
    partition_idx=partition_idx(1:end-1);
  end
  partition_idx=[0 partition_idx 0];
  start=1+find(~partition_idx(1:(end-1)) &  partition_idx(2:end))-1;
  stop =  find( partition_idx(1:(end-1)) & ~partition_idx(2:end))-1;
  during{i}={};  not_during{i}={};
  raw_during{i}={};  raw_not_during{i}={};
  if(length(start)>0)
    for j=1:length(start)
      raw_during{i}{j}=feature_data.data{i}(start(j):stop(j));
      [m,f,c]=mode(raw_during{i}{j});
      during{i}{j}={c{1} f./length(raw_during{i}{j})*100};
    end
    if(length(start)>1)
      for j=1:(length(start)-1)
        raw_not_during{i}{j}=feature_data.data{i}(stop(j):start(j+1));
        [m,f,c]=mode(raw_not_during{i}{j});
        not_during{i}{j}={c{1} f./length(raw_not_during{i}{j})*100};
      end
    end
  end
end



% ---
function ret_val=print_stats2(arg)

sprintf('%d,',arg{1});
ret_val=[ans(1:(end-1)) ' (' num2str(arg{2},3) '%)']; %#ok<COLND>



% ---
function [table_data raw_table_data]=plot_social(experiment_value,experiment_list,...
    behavior_value,behavior_list,behavior_logic,behavior_value2,behavior_list2,...
    feature_value,feature_list,color)

during={};  not_during={};
parfor e=1:length(experiment_value)
%for e=1:length(experiment_value)
  behavior_data=load(fullfile(experiment_list{experiment_value(e)},...
        [behavior_list{behavior_value} '.mat']));
  if(behavior_logic>1)
    behavior_data2=load(fullfile(experiment_list{experiment_value(e)},...
        [behavior_list2{behavior_value2} '.mat']));
  else
    behavior_data2=[];
  end
  feature_data=load(fullfile(experiment_list{experiment_value(e)},'perframe',...
        [feature_list{feature_value} '.mat']));

  [during{e} not_during{e} raw_during{e} raw_not_during{e}]=...
      calculate_social(behavior_data,behavior_logic,behavior_data2,feature_data);
end

k=1;
table_data={};
raw_table_data={};
for e=1:length(experiment_value)
  for i=1:length(during{e})
    table_data{k,1}=['Exp #' color num2str(experiment_value(e)) ', Indi #' num2str(i)];
    if ~isempty(during{e}{i})
      cellfun(@(x) x{1},during{e}{i},'uniformoutput',false);
      raw_table_data{k,2}=cat(1,ans{:});
      [m,f,c]=mode(raw_table_data{k,2});
      sprintf('%d,',c{1});
      table_data{k,2}=[ans(1:(end-1)) ' (' num2str(f/length(raw_table_data{k,2})*100,3) '%)'];
    end
    if ~isempty(not_during{e}{i})
      cellfun(@(x) x{1},not_during{e}{i},'uniformoutput',false);
      raw_table_data{k,3}=cat(1,ans{:});
      [m,f,c]=mode(raw_table_data{k,3});
      sprintf('%d,',c{1});
      table_data{k,3}=[ans(1:(end-1)) ' (' num2str(f/length(raw_table_data{k,3})*100,3) '%)'];
    end
    if ~isempty(during{e}{i})
      tmp=cellfun(@print_stats2,during{e}{i},'uniformoutput',false);
      [table_data{k,4:2:(2+2*length(during{e}{i}))}]=deal(tmp{:});
      [raw_table_data{k,4:2:(2+2*length(during{e}{i}))}]=deal(raw_during{e}{i}{:});
    end
    if ~isempty(not_during{e}{i})
      tmp=cellfun(@print_stats2,not_during{e}{i},'uniformoutput',false);
      [table_data{k,5:2:(3+2*length(not_during{e}{i}))}]=deal(tmp{:});
      [raw_table_data{k,5:2:(3+2*length(not_during{e}{i}))}]=deal(raw_not_during{e}{i}{:});
    end
    k=k+1;
  end
end


% --- Executes on button press in SocialStats.
function SocialStats_Callback(hObject, eventdata, handles)
% hObject    handle to SocialStats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


set(handles.Status,'string','Thinking...','foregroundcolor','b');
set(handles.figure1,'pointer','watch');
drawnow;

experiment_value=get(handles.ExperimentList,'Value');
experiment_list=get(handles.ExperimentList,'String');
experiment_value2=get(handles.ExperimentList2,'Value');
experiment_list2=get(handles.ExperimentList2,'String');
behavior_value=get(handles.BehaviorList,'Value');
behavior_list=get(handles.BehaviorList,'String');
behavior_logic=get(handles.BehaviorLogic,'Value');
behavior_value2=get(handles.BehaviorList2,'Value');
behavior_list2=get(handles.BehaviorList2,'String');
feature_value=get(handles.FeatureList,'Value');
feature_list=get(handles.FeatureList,'String');

if(~strncmp(feature_list(feature_value),'closestfly_',11))
  feature_value=find(cellfun(@(x) strncmp(x,'closestfly_',11),feature_list),1);
  if(isempty(feature_value))
    error('no closestfly feature found');
  end
  set(handles.FeatureList,'Value',feature_value);
end

%axes(handles.Axes);  cla;  hold on;
figure;  hold on;

if(length(experiment_value)>0)
  [table_data raw_table_data]=plot_social(experiment_value,experiment_list,...
      behavior_value,behavior_list,behavior_logic,behavior_value2,behavior_list2,...
      feature_value,feature_list,'R');
end
if(length(experiment_value2)>0)
  [table_data2 raw_table_data2]=plot_social(experiment_value2,experiment_list2,...
      behavior_value,behavior_list,behavior_logic,behavior_value2,behavior_list2,...
      feature_value,feature_list,'B');
end

if((length(experiment_value)>0) && (length(experiment_value2)>0))
  tmp(1:2:2*size(table_data,1),1:size(table_data, 2))=table_data;
  tmp(2:2:2*size(table_data2,1),1:size(table_data2,2))=table_data2;
  raw_tmp(1:2:2*size(raw_table_data,1),1:size(raw_table_data, 2))=raw_table_data;
  raw_tmp(2:2:2*size(raw_table_data2,1),1:size(raw_table_data2,2))=raw_table_data2;
  ii=max(size(table_data,2),size(table_data2,2));
  set(handles.Table,'Data',tmp);
  set(handles.Table,'RowStriping','on','BackgroundColor',[1 0.8 0.8; 0.8 0.8 1]);
elseif(length(experiment_value)>0)
  tmp=table_data;
  raw_tmp=raw_table_data;
  ii=size(table_data,2);
  set(handles.Table,'Data',tmp);
  set(handles.Table,'RowStriping','off','BackgroundColor',[1 0.8 0.8]);
else
  tmp=table_data2;
  raw_tmp=raw_table_data2;
  ii=size(table_data2,2);
  set(handles.Table,'Data',tmp);
  set(handles.Table,'RowStriping','off','BackgroundColor',[0.8 0.8 1]);
end

handles.raw_table_data=raw_tmp;
tmp={'Individual' 'Bout Mode' 'Inter Bout Mode'};
cellstr(num2str((1:(((ii-3)/2)+1))','B #%d'))';
[tmp{4:2:ii}]=deal(ans{:});
cellstr(num2str((1: ((ii-3)/2)   )','IB #%d'))';
[tmp{5:2:ii}]=deal(ans{:});
%for i=1:(ii-3)/2
%  tmp{2*i+2}=['B #' num2str(i)];
%  tmp{2*i+3}=['IB #' num2str(i)];
%end
set(handles.Table,'ColumnName',tmp);
set(handles.Table,'ColumnWidth',{100 75 100});

handles.table='social_stats';
guidata(hObject,handles);

set(handles.Status,'string','Ready.','foregroundcolor','g');
set(handles.figure1,'pointer','arrow');
drawnow;


% --- 
function CellSelectionCallback(hObject,eventdata)

if(size(eventdata.Indices,1)==0)  return;  end

handles2=guidata(gcf);
handles=guidata(handles2.figure1);

if(strcmp(handles2.table,'social_stats'))
  if(eventdata.Indices(end,2)==1)  return;  end

  figure;  hold on;
  tmp={};
    experiment_value=get(handles.ExperimentList,'Value');
    experiment_value2=get(handles.ExperimentList2,'Value');
    if((length(experiment_value)>0) && (length(experiment_value2)>0))
      if(mod(eventdata.Indices(end,1),2))
        color='r';
      else
        color='b';
      end
    elseif(length(experiment_value)>0)
      color='r';
    else
      color='b';
    end

    data=handles.raw_table_data{eventdata.Indices(end,1),eventdata.Indices(end,2)};
    tmp=unique(handles.raw_table_data{eventdata.Indices(end,1),2});
    hist(data,min(tmp):max(tmp));

  xlabel('closest fly (#)');
  axis tight;

elseif(strcmp(handles2.table,'histogram'))
  switch(eventdata.Indices(end,2))
    case {1,3}
      group=handles2.table_data(eventdata.Indices(end,1),eventdata.Indices(end,2));
      switch(group)
        case(-1)
          questdlg('Remove all during / all-frames comparisons from table?','','Yes','No','No');
        case(-2)
          questdlg('Remove all during / not-during comparisons from table?','','Yes','No','No');
        otherwise
          questdlg(['Remove all ' handles.grouplist{group} ' groups from table?'],'','Yes','No','No');
      end
      if(strcmp(ans,'No'))  return;  end
      tmp=get(handles2.Table,'Data');
      idx=find((handles2.table_data(:,1)~=group)&(handles2.table_data(:,3)~=group));
      set(handles2.Table,'Data',tmp(idx,:));
      handles2.table_data=handles2.table_data(idx,:);
    case {2,4}
      thresh=handles2.table_data(eventdata.Indices(end,1),eventdata.Indices(end,2));
      questdlg(['Remove all rows for which n or n2 is less than or equal to ' num2str(thresh) ...
          ' from table?'],'','Yes','No','No');
      if(strcmp(ans,'No'))  return;  end
      tmp=get(handles2.Table,'Data');
      idx=find((handles2.table_data(:,2)>thresh)&(handles2.table_data(:,4)>thresh));
      set(handles2.Table,'Data',tmp(idx,:));
      handles2.table_data=handles2.table_data(idx,:);
    case {5}
      not_txt='';  if(handles2.table_data(eventdata.Indices(end,1),5)<0)  not_txt='not ';  end
      questdlg(['Remove all ' not_txt handles.behaviorlist{abs(handles2.table_data(eventdata.Indices(end,1),5))} ...
          ' behaviors from table?'],'','Yes','No','No');
      if(strcmp(ans,'No'))  return;  end
      tmp=get(handles2.Table,'Data');
      idx=find(handles2.table_data(:,5)~=handles2.table_data(eventdata.Indices(end,1),5));
      set(handles2.Table,'Data',tmp(idx,:));
      handles2.table_data=handles2.table_data(idx,:);
    case {6}
      questdlg(['Remove all ' handles.featurelist{handles2.table_data(eventdata.Indices(end,1),6)} ...
          ' features from table?'],'','Yes','No','No');
      if(strcmp(ans,'No'))  return;  end
      tmp=get(handles2.Table,'Data');
      idx=find(handles2.table_data(:,6)~=handles2.table_data(eventdata.Indices(end,1),6));
      set(handles2.Table,'Data',tmp(idx,:));
      handles2.table_data=handles2.table_data(idx,:);
    case {7}
      handles.analysis='feature_histogram';
      handles.behaviornot=round(0.5*(1+-sign(handles2.table_data(eventdata.Indices(end,1),5))));
      handles.behaviorvalue=max(1,abs(handles2.table_data(eventdata.Indices(end,1),5)));
      handles.behaviorlogic=1;
      handles.featurevalue=handles2.table_data(eventdata.Indices(end,1),6);
      handles.individual=1;
      handles.comparison=1;
      if(handles2.table_data(eventdata.Indices(end,1),3)<0)
        handles.comparison=-handles2.table_data(eventdata.Indices(end,1),3);
      end
      button_comparison_set(handles);
      %FeatureHistogram_Callback(hObject, eventdata, handles);
      handles=feature_histogram_plot(handles);
    case {8}
      if(isnan(handles2.table_data(eventdata.Indices(end,1),eventdata.Indices(end,2))))
        questdlg(['Remove all rows for which d''-AF is NaN?'],'','Yes','No','No');
        if(strcmp(ans,'No'))  return;  end
        tmp=get(handles2.Table,'Data');
        idx=find(~isnan(handles2.table_data(:,8)));
        set(handles2.Table,'Data',tmp(idx,:));
        handles2.table_data=handles2.table_data(idx,:);
      else
        crit=handles2.table_data(eventdata.Indices(end,1),7) ./ handles2.table_data(eventdata.Indices(end,1),8);
        questdlg(['Remove all rows for which the ratio of d'' to d''-AF is less than ' num2str(crit) '?'],...
            '','Yes','No','No');
        if(strcmp(ans,'No'))  return;  end
        tmp=get(handles2.Table,'Data');
        idx=find((handles2.table_data(:,7)./handles2.table_data(:,8))>=crit);
        set(handles2.Table,'Data',tmp(idx,:));
        handles2.table_data=handles2.table_data(idx,:);
      end
  end
  update_figure(handles);

end

guidata(handles.figure1,handles);
guidata(handles2.figure2,handles2);


% ---
function menu_classify_forcecompute_set(handles)

if(handles.classify_forcecompute)
  set(handles.MenuClassifyForceCompute,'Checked','on');
else
  set(handles.MenuClassifyForceCompute,'Checked','off');
end


% --------------------------------------------------------------------
function MenuClassifyForceCompute_Callback(hObject, eventdata, handles)
% hObject    handle to MenuClassifyForceCompute (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.classify_forcecompute=~handles.classify_forcecompute;
menu_classify_forcecompute_set(handles);
guidata(hObject,handles);


% --------------------------------------------------------------------
function MenuClassify_Callback(hObject, eventdata, handles)
% hObject    handle to MenuClassify (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in CloseAll.
function CloseAll_Callback(hObject, eventdata, handles)
% hObject    handle to CloseAll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

questdlg('Close all figures?','','Yes','No','No');
if(strcmp(ans,'No'))  return;  end
hh=findobj('type','figure');
for h=1:length(hh)
  if(hh(h)~=handles.figure1)  close(hh(h));  end
end


%[file,path]=uiputfile('*.txt','Select file to save table to');
%if (file==0)  return;  end
%fid=fopen(fullfile(path,file),'w');
%tmp=get(handles.Table,'columnname');
%tmp=transpose(tmp);
%tmp2=get(handles.Table,'data');
%tmp2(2:end+1,:)=tmp2;
%tmp2(1,:)=tmp;
%tmp2=cellfun(@char,tmp2,'uniformoutput',false);
%for i=1:size(tmp2,1)
%  for j=1:size(tmp2,2)
%    fprintf(fid,'%s\t',char(tmp2(i,j)));
%  end
%  fprintf(fid,'\r\n');
%end
%fclose(fid);


% --------------------------------------------------------------------
function MenuHelp_Callback(hObject, eventdata, handles)
% hObject    handle to MenuHelp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function MenuHelpDocumentation_Callback(hObject, eventdata, handles)
% hObject    handle to MenuHelpDocumentation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

html_file = 'http://jaaba.sourceforge.net/PlottingResults.html';
if isdeployed,
  [stat,msg] = myweb_nocheck(html_file);
  if stat ~= 0,
    errordlg({'Please see documentation at http://jaaba.sourceforge.net'
      'Error opening webpage within MATLAB:'
      msg});
  end
else
  web('-browser',html_file);
end


% --------------------------------------------------------------------
function MenuHelpBugReport_Callback(hObject, eventdata, handles)
% hObject    handle to MenuHelpBugReport (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

{...
'Having problems?  Before contacting us first try:'...
''...
'1. Going to the file menu and choosing update.  Does that fix it?'...
'2. Going to the file menu and choosing reset.  Add all of your experiments and try again.  Fixed?'...
'3. If you quit out of JAABAPlot and restart, does it still happen?  Does quiting out of Matlab help?  Rebooting your computer?'...
''...
'If it still doesn''t work, please e-mail our Google newsgroup with the following information:'...
''...
'1. Are you using JAABA on Windows, Mac, or Linux?'...
'2. Are you using the compiled version or running it within Matlab?'...
'3. What version of JAABA are you using?'...
'4. What exactly is the problem?  What did you expect JAABAPlot to do?'...
'5. What is the sequence of steps leading up to the problem?  That is, the order of buttons you push.'...
'6. Does it happen every time?'...
'7. JAABAPlot creates a file called most_recent_config.mat.  Please attach it.'...
'8. Please also cut and paste any red error messages on the Matlab command line.  Or a screenshot if you''re using the compiled version.'...
''...
'JAABAPlot is still a work in progress and we appreciate your feedback.'...
};
msgbox(ans,'Bug Report');


% --------------------------------------------------------------------
function MenuFileReset_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFileReset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

questdlg('Reset configuration to default?','','Yes','No','No');
if(strcmp(ans,'No'))  return;  end
handles=initialize(handles);
update_figure(handles);
set(handles.figure1,'pointer','arrow');
guidata(hObject, handles);


% --------------------------------------------------------------------
function MenuFileUpdate_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFileUpdate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.Status,'string','Thinking...','foregroundcolor','b');
set(handles.figure1,'pointer','watch');
drawnow;

handles=update_experiment_data(handles,true,true,true);
update_figure(handles);
guidata(hObject, handles);

set(handles.Status,'string','Ready.','foregroundcolor','g');
set(handles.figure1,'pointer','arrow');
drawnow;


% --------------------------------------------------------------------
function MenuFileLoad_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFileLoad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[file,path]=uigetfile('*.mat','Select configuration file to open');
if(isnumeric(file) && isnumeric(path) && (file==0) && (path==0))  return;  end
handles=load_configuration_file(fullfile(path,file),hObject,eventdata,handles);
update_figure(handles);
guidata(hObject, handles);


% --------------------------------------------------------------------
function MenuFileSave_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFileSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[file,path]=uiputfile('*.mat','Select file to save configuration to');
if(isnumeric(file) && isnumeric(path) && (file==0) && (path==0))  return;  end
save(fullfile(path,file),'handles');


% --------------------------------------------------------------------
function MenuFile_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% ---
function figure_stats_callback(src,evt)

handles=guidata(src);

len=max(cellfun(@(x) length(regexprep(num2str(x),'<[^<]*>','*')),handles.statistics));
hf=figure('menubar','none','toolbar','none','numbertitle','off',...
    'name',['statistics for Figure ' num2str(gcbf)]);
ht=uitable('data',handles.statistics,'columnwidth',num2cell(8*len),...
    'rowname',[],'columnname',[],'rowstriping','off');
extT=get(ht,'extent');
posF=get(hf,'position');
%set(hf,'position',[posF(1) posF(2) 16*(posF(4)<extT(4))+extT(3) posF(4)]);
set(hf,'position',[posF(1) posF(2) min(posF(3),16*(posF(4)<extT(4))+extT(3)) min(posF(4),extT(4))]);
%set(ht,'position',[0 0 16*(posF(4)<extT(4))+extT(3) posF(4)]);
set(ht,'units','normalized','position',[0 0 1 1]);



% ---
function figure_params_callback(src,evt)

handles=guidata(src);

CT={'Mean' 'Median' 'Mode'};
D={'Std. Dev.' 'Std. Err.' '5%-95%' '25%-75%'};
xoffset={'none', 'start', 'min(start)'};

tmp={};

switch handles.type
  case 'feature histogram'
    tmp{end+1}=['style = ' handles.featurehistogram_stylelist{handles.featurehistogram_style}];
    tmp{end+1}=['style2 = ' handles.featurehistogram_stylelist2{handles.featurehistogram_style2}];
    tmp{end+1}=['allframes=' num2str(handles.comparison==1)];
    tmp{end+1}=['notduring=' num2str(handles.comparison==2)];
    tmp{end+1}=['logbinsize=' num2str(handles.logbinsize)];
    tmp{end+1}=['nbins=' num2str(handles.nbins)];
    tmp{end+1}=['central tendency = ' CT{handles.centraltendency}];
    tmp{end+1}=['dispersion = '  D{handles.dispersion}];
  case 'feature time series'
    tmp{end+1}=['style = ' handles.featuretimeseries_stylelist{handles.featuretimeseries_style}];
    tmp{end+1}=['style2 = ' handles.featuretimeseries_stylelist2{handles.featuretimeseries_style2}];
    tmp{end+1}=['subtractmean=' num2str(handles.subtractmean)];
    tmp{end+1}=['conv. width = '  num2str(handles.convolutionwidth) ' sec'];
    tmp{end+1}=['radius=' num2str(handles.windowradius)];
    tmp{end+1}=['x-offset = '  xoffset{handles.xoffset}];
    tmp{end+1}=['central tendency = ' CT{handles.centraltendency}];
    tmp{end+1}=['dispersion = '  D{handles.dispersion}];
  case 'behavior bar chart'
    tmp{end+1}=['style = ' handles.behaviorbarchart_stylelist{handles.behaviorbarchart_style}];
    tmp{end+1}=['central tendency = ' CT{handles.centraltendency}];
    tmp{end+1}=['dispersion = '  D{handles.dispersion}];
    tmp2='';  if(handles.behaviornormalizenot)  tmp2='not ';  end
    tmp3='all frames';  if(handles.behaviorvalue3>1)  tmp3=handles.behaviorlist(handles.behaviorvalue3-1);  end
    tmp{end+1}=['normalize = ' tmp2 char(tmp3)];
  case 'behavior time series'
    tmp{end+1}=['style = ' handles.behaviortimeseries_stylelist{handles.behaviortimeseries_style}];
    tmp{end+1}=['conv. width = '  num2str(handles.convolutionwidth) ' sec'];
    tmp{end+1}=['x-offset = '  xoffset{handles.xoffset}];
    tmp{end+1}=['central tendency = ' CT{handles.centraltendency}];
    tmp{end+1}=['dispersion = '  D{handles.dispersion}];
    tmp2='';  if(handles.behaviornormalizenot)  tmp2='not ';  end
    tmp3='all frames';  if(handles.behaviorvalue3>1)  tmp3=handles.behaviorlist(handles.behaviorvalue3-1);  end
    tmp{end+1}=['normalize = ' tmp2 char(tmp3)];
  case 'bout stats'
    tmp{end+1}=['style = ' handles.boutstats_stylelist{handles.boutstats_style}];
    tmp{end+1}=['style2 = ' handles.boutstats_stylelist2{handles.boutstats_style2}];
    tmp{end+1}=['central tendency = ' CT{handles.centraltendency}];
    tmp{end+1}=['dispersion = '  D{handles.dispersion}];
end
tmp{end+1}=['minimum trajectory length = '  num2str(handles.minimumtrajectorylength)];
tmp{end+1}='';

for g=1:length(handles.grouplist)
  if(length(handles.experimentvalue{g})==0)  continue;  end
  tmp{end+1}=['group ' handles.grouplist{g}];
  for e=1:length(handles.experimentvalue{g})
    [~,tmp{end+1},~]=fileparts(handles.experimentlist{g}{handles.experimentvalue{g}(e)});
  end
  tmp{end+1}='';
end

hf=figure('menubar','none','toolbar','none','numbertitle','off',...
    'name',['parameters for Figure ' num2str(gcbf)]);
ht=uitable('data',tmp','columnwidth',num2cell(8*max(cellfun(@length,tmp))),...
    'rowname',[],'columnname',[],'rowstriping','off');
extT=get(ht,'extent');
posF=get(hf,'position');
set(hf,'position',[posF(1) posF(2) min(posF(3),16*(posF(4)<extT(4))+extT(3)) min(posF(4),extT(4))]);
set(ht,'units','normalized','position',[0 0 1 1]);



%BOUT STATISTICS -- not converted to multivariate yet
%
%Create a table of bout lengths (BL) and inter-bout lengths (IBL) using the
%"Bout Stats" button.  As for behavior stats, each row is a behavior;  the
%columns break it down by sex and individual;  the experiment groups are
%striped;  logical operators can be applied; and the Prefs contextual menu
%controls statistics.
%
%Selecting a cell plots a histogram, normalized to unit area, of bout and
%inter-bout lengths.  The "LogY" and "Stats" buttons can be used to scale
%the axis and overlay summary statistics, respectively.  Multiple selected
%cells are overlayed.
%
%
%SOCIAL STATISTICS -- not converted to multivariate yet
%
%Create a table of the closest individual using the "Social Stats" button.
%Choose which behavior to analysis using the pull-down menu in the Behavior
%panel, and which metric to use for closest using the pull-down menu in the
%Feature panel.  In the resulting table, each row is an individual, and the
%columns show the mode and percentage of frames of the closest individual
%for each bout, as well as the overall mode of the per-bout modes.
%
%Selecting a cell in the table plots a histogram of the closest indidividual.



function MinimumTrajectoryLength_Callback(hObject, eventdata, handles)
% hObject    handle to MinimumTrajectoryLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MinimumTrajectoryLength as text
%        str2double(get(hObject,'String')) returns contents of MinimumTrajectoryLength as a double

handles.minimumtrajectorylength=str2num(get(hObject,'String'));
handles.interestingfeaturehistograms_cache=[];
handles.interestingfeaturetimeseries_cache=[];
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function MinimumTrajectoryLength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MinimumTrajectoryLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ConvolutionWidth_Callback(hObject, eventdata, handles)
% hObject    handle to ConvolutionWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ConvolutionWidth as text
%        str2double(get(hObject,'String')) returns contents of ConvolutionWidth as a double

%inputdlg({'Convolution width:'},'',1,{num2str(handles.convolutionwidth)});
%if(isempty(ans))  return;  end
handles.convolutionwidth=str2num(get(hObject,'String'));
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function ConvolutionWidth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ConvolutionWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function PValue_Callback(hObject, eventdata, handles)
% hObject    handle to PValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PValue as text
%        str2double(get(hObject,'String')) returns contents of PValue as a double

%inputdlg({'P value:'},'',1,{num2str(handles.pvalue)});
%if(isempty(ans))  return;  end
handles.pvalue=str2num(get(hObject,'String'));
guidata(hObject,handles);




% --- Executes during object creation, after setting all properties.
function PValue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function NBins_Callback(hObject, eventdata, handles)
% hObject    handle to NBins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NBins as text
%        str2double(get(hObject,'String')) returns contents of NBins as a double

%handles.nbins=str2num(char(inputdlg({'Number of bins:'},'',1,...
%    {num2str(handles.nbins)})));
handles.nbins=str2num(get(hObject,'String'));
guidata(hObject,handles);




% --- Executes during object creation, after setting all properties.
function NBins_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NBins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function WindowRadius_Callback(hObject, eventdata, handles)
% hObject    handle to WindowRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of WindowRadius as text
%        str2double(get(hObject,'String')) returns contents of WindowRadius as a double

handles.windowradius=str2num(get(hObject,'String'));
guidata(hObject,handles);



% --- Executes during object creation, after setting all properties.
function WindowRadius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to WindowRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in StatisticsList.
function StatisticsList_Callback(hObject, eventdata, handles)
% hObject    handle to StatisticsList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns StatisticsList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from StatisticsList

handles.statistic=get(handles.Statistic,'value');


% --- Executes during object creation, after setting all properties.
function StatisticsList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to StatisticsList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in StyleList.
function StyleList_Callback(hObject, eventdata, handles)
% hObject    handle to StyleList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns StyleList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from StyleList

get(hObject,'Value');
switch(handles.analysis)
  case 'feature_histogram'
    handles.featurehistogram_style2=ans;
  case 'feature_timeseries'
    handles.featuretimeseries_style=ans;
  case 'behavior_barchart'
    handles.behaviorbarchart_style=ans;
  case 'behavior_timeseries'
    handles.behaviortimeseries_style=ans;
  case 'bout_stats'
    handles.boutstats_style=ans;
  case 'interesting_feature_histograms'
end
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function StyleList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to StyleList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in StyleList2.
function StyleList2_Callback(hObject, eventdata, handles)
% hObject    handle to StyleList2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns StyleList2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from StyleList2

get(hObject,'Value');
switch(handles.analysis)
  case 'feature_histogram'
    handles.featurehistogram_style=ans;
  case 'feature_timeseries'
    handles.featuretimeseries_style2=ans;
    update_figure(handles);
  case 'behavior_barchart'
  case 'behavior_timeseries'
  case 'bout_stats'
    handles.boutstats_style2=ans;
  case 'interesting_feature_histograms'
end
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function StyleList2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to StyleList2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% ---
function button_comparison_set(handles)

set(handles.AllFrames,'Value',0);
set(handles.NotDuring,'Value',0);
switch(handles.comparison)
  case(1), set(handles.AllFrames,'Value',1);
  case(2), set(handles.NotDuring,'Value',1);
end


% --- Executes on button press in AllFrames.
function AllFrames_Callback(hObject, eventdata, handles)
% hObject    handle to AllFrames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if(handles.comparison==1)
  handles.comparison=0;
else
  handles.comparison=1;
end
button_comparison_set(handles);
guidata(hObject,handles);


% --- Executes on button press in AllFrames.
function NotDuring_Callback(hObject, eventdata, handles)
% hObject    handle to AllFrames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if(handles.comparison==2)
  handles.comparison=0;
else
  handles.comparison=2;
end
button_comparison_set(handles);
guidata(hObject,handles);


% --- Executes on button press in OmitInf.
function OmitInf_Callback(hObject, eventdata, handles)
% hObject    handle to OmitInf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.omitinf=~handles.omitinf;
%button_omitinf_set(handles.omitinf);
%handles.interestingfeaturehistograms_cache=[];
guidata(hObject,handles);


% --- Executes on button press in OmitNaN.
function OmitNaN_Callback(hObject, eventdata, handles)
% hObject    handle to OmitNaN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.omitnan=~handles.omitnan;
%button_omitnan_set(handles.omitnan);
%handles.interestingfeaturehistograms_cache=[];
guidata(hObject,handles);


% --- Executes on button press in AbsDPrime.
function AbsDPrime_Callback(hObject, eventdata, handles)
% hObject    handle to AbsDPrime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.absdprime=~handles.absdprime;
%button_absdprime_set(handles.absdprime);
%handles.interestingfeaturehistograms_cache=[];
guidata(hObject,handles);


% --- Executes on button press in ForceCompute.
function ForceCompute_Callback(hObject, eventdata, handles)
% hObject    handle to ForceCompute (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in SubtractMean.
function SubtractMean_Callback(hObject, eventdata, handles)
% hObject    handle to SubtractMean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.subtractmean=~handles.subtractmean;
%button_subtractmean_set(handles.subtractmean);
guidata(hObject,handles);


% --- Executes on button press in LogBinSize.
function LogBinSize_Callback(hObject, eventdata, handles)
% hObject    handle to LogBinSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.logbinsize=~handles.logbinsize;
%button_logbinsize_set(handles.logbinsize);
guidata(hObject,handles);


% --- Executes on selection change in XOffset.
function XOffset_Callback(hObject, eventdata, handles)
% hObject    handle to XOffset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns XOffset contents as cell array
%        contents{get(hObject,'Value')} returns selected item from XOffset

handles.xoffset=get(handles.Xoffset,'value');
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function XOffset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to XOffset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Statistic.
function CentralTendency_Callback(hObject, eventdata, handles)
% hObject    handle to Statistic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Statistic contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Statistic

handles.centraltendency=get(handles.CentralTendency,'value');
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function CentralTendency_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Statistic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Statistic.
function Dispersion_Callback(hObject, eventdata, handles)
% hObject    handle to Statistic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Statistic contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Statistic

handles.dispersion=get(handles.Dispersion,'value');
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function Dispersion_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Statistic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
