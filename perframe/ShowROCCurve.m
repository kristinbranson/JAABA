function varargout = ShowROCCurve(varargin)
% SHOWROCCURVE MATLAB code for ShowROCCurve.fig
%      SHOWROCCURVE, by itself, creates a new SHOWROCCURVE or raises the existing
%      singleton*.
%
%      H = SHOWROCCURVE returns the handle to a new SHOWROCCURVE or the handle to
%      the existing singleton*.
%
%      SHOWROCCURVE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SHOWROCCURVE.M with the given input arguments.
%
%      SHOWROCCURVE('Property','Value',...) creates a new SHOWROCCURVE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ShowROCCurve_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ShowROCCurve_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ShowROCCurve

% Last Modified by GUIDE v2.5 30-Jan-2012 15:54:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
  'gui_Singleton',  gui_Singleton, ...
  'gui_OpeningFcn', @ShowROCCurve_OpeningFcn, ...
  'gui_OutputFcn',  @ShowROCCurve_OutputFcn, ...
  'gui_LayoutFcn',  [] , ...
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
end

% --- Executes just before ShowROCCurve is made visible.
function ShowROCCurve_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ShowROCCurve (see VARARGIN)

% Choose default command line output for ShowROCCurve
handles.output = hObject;
handles.labels = varargin{1};
handles.scores = varargin{2};
handles.JLDObj = varargin{3};
handles.JLabelObject = varargin{4};
handles.JLabelHandle = guidata(handles.JLabelObject);

numBins = 21;
handles.scoreNorm = handles.JLDObj.windowdata.scoreNorm;
if isnan(handles.scoreNorm); handles.scoreNorm = 1; end;
pos = handles.labels==1;
neg = ~pos;
bins = linspace(-handles.scoreNorm,handles.scoreNorm,numBins);
bins = [-inf bins(2:end-1) inf];
histPos = histc(handles.scores(pos),bins);
histNeg = histc(handles.scores(neg),bins);
histPos(end) = []; histNeg(end) = [];

handles.posColor = handles.JLabelHandle.labelcolors(1,:);
handles.negColor = handles.JLabelHandle.labelcolors(2,:);

handles.thres1 = round(handles.JLDObj.GetConfidenceThreshold(1)*(numBins-1)/2)/( (numBins-1)/2);
handles.thres2 = -round(handles.JLDObj.GetConfidenceThreshold(2)*(numBins-1)/2)/( (numBins-1)/2);


xLocs = linspace(-1+1/(numBins-1),1-1/(numBins-1),numBins-1);
hBar = bar(handles.axes1, xLocs,[histPos histNeg],'BarWidth',1.5);
set(hBar(1),'FaceColor',handles.posColor);
set(hBar(2),'FaceColor',handles.negColor);
ylim = get(handles.axes1,'ylim');
if(ylim(1)<0); 
  ylim(1) = 0; 
  set(handles.axes1,'ylim',ylim);
end
set(handles.axes1,'xlim',[-1 1])
shortylim = ylim;
shortylim(1) = ylim(1) + 0.001;
shortylim(2) = ylim(2) - 0.001;
linepos = imline(handles.axes1,[handles.thres1 handles.thres1]+0.005,shortylim);
lineneg = imline(handles.axes1,[handles.thres2 handles.thres2]-0.005,shortylim);
setColor(linepos,handles.posColor);
setColor(lineneg,handles.negColor);

possibleXLocs = linspace(-1,1,numBins);
possibleXLocs = possibleXLocs(2:end-1);
api = iptgetapi(linepos);
fcn = GetLineConstraintFcn(get(handles.axes1,'XLim'),...
  get(handles.axes1,'YLim'),[[handles.thres1 handles.thres1]; ylim],...
  possibleXLocs,hObject,1);
api.setPositionConstraintFcn(fcn);

api = iptgetapi(lineneg);
fcn = GetLineConstraintFcn(get(handles.axes1,'XLim'),...
  get(handles.axes1,'YLim'),[[handles.thres2 handles.thres2]; ylim],...
  possibleXLocs,hObject,2);
api.setPositionConstraintFcn(fcn);
UpdateText(handles);
% % Update handles structure
guidata(hObject, handles);

% UIWAIT makes ShowROCCurve wait for user response (see UIRESUME)
% uiwait(handles.figure1);
end

function UpdateText(handles)
pos = handles.labels>0;
posCorrect = sum( (handles.scores(pos)/handles.scoreNorm)>=handles.thres1);
posIncorrect = sum((handles.scores(pos)/handles.scoreNorm)<=handles.thres2);
posAbstain = sum( (handles.scores(pos)/handles.scoreNorm)<handles.thres1 & ...
  (handles.scores(pos)/handles.scoreNorm)>handles.thres2);
negCorrect = sum(  (handles.scores(~pos)/handles.scoreNorm)<=handles.thres2);
negIncorrect = sum((handles.scores(~pos)/handles.scoreNorm)>=handles.thres1);
negAbstain = sum(  (handles.scores(~pos)/handles.scoreNorm)<handles.thres1 & ...
                   (handles.scores(~pos)/handles.scoreNorm)>handles.thres2);
textStr = '';
textStr = sprintf('                      Predicted Positive      Abstained        Predicted Negative\n');
textStr = sprintf('%s Labeled Pos      %5d                 %5d                %5d      \n',...
  textStr,posCorrect, posAbstain,posIncorrect);
textStr = sprintf('%s Labeled Neg      %5d                 %5d                %5d      \n',...
  textStr,negIncorrect, negAbstain,negCorrect);
set(handles.text1,'String',textStr);
end


% --- Outputs from this function are returned to the command line.
function varargout = ShowROCCurve_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

end

function fcn = GetLineConstraintFcn(xlim,ylim,origPosition,xLocs,hObject,thresNum)

fcn = @constrainLineToRect;

%Store previous position matrix for use in constrainLineToRect.
line_pos_last = origPosition;

%-----------------------------------------
  function new_pos = constrainLineToRect(pos)

    handles = guidata(hObject);
    otherLineThres = handles.(sprintf('thres%d',3-thresNum));
    previous_position_cached = ~isempty(line_pos_last);
    
    is_end_point_drag = previous_position_cached &&...
      (any(pos(:,1) == line_pos_last(:,1)) &&...
      any(pos(:,2) == line_pos_last(:,2)));
    
    if is_end_point_drag
      new_pos = line_pos_last;
      %             new_pos = [constrainPointToRect(pos(1,:));constrainPointToRect(pos(2,:))];
    else
      %Apply correction made to first end point to both end points
      constrained_p1 = constrainPointToRect(pos(1,:));
      v1 = constrained_p1 - pos(1,:);
      temp_pos = pos + [v1; v1];
      
      % Now reconstrain both end points according to correction made
      % to second endpoint of partially constrained line.
      constrained_p2 = constrainPointToRect(temp_pos(2,:));
      v2 = constrained_p2 - temp_pos(2,:);
      new_pos = temp_pos + [v2; v2];
      
      if thresNum ==1
        if new_pos(1)<0 %otherLineThres)
          new_pos(:,1) = 0; %otherLineThres;
        end
      else
        if new_pos(1)>0 %otherLineThres)
          new_pos(:,1) = 0; %otherLineThres;
        end        
      end
      
      [~,closestXlocs] = min( abs(new_pos(1)-xLocs));
      new_pos(:,1) = xLocs(closestXlocs);
      handles.(sprintf('thres%d',thresNum)) = new_pos(1);
        
      UpdateText(handles);
      guidata(hObject,handles);
      if thresNum ==1
        new_pos(:,1) = new_pos(:,1)+0.005;
      else
        new_pos(:,1) = new_pos(:,1)-0.005;
      end
    end
    
    line_pos_last = new_pos;
  end

  function new_pos = constrainPointToRect(pos)
    
    x_candidate = pos(1);
    y_candidate = pos(2);
    
    x_new = min( xlim(2), max(x_candidate, xlim(1)) );
    y_new = min( ylim(2), max(y_candidate, ylim(1)) );
    
    new_pos = [x_new y_new];
    
  end

  end


% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure1);
end

% --- Executes on button press in pushbutton_apply.
function pushbutton_apply_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_apply (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.JLDObj.SetConfidenceThreshold(handles.thres1,1);
handles.JLDObj.SetConfidenceThreshold(-handles.thres2,2);
handles.JLabelHandle = guidata(handles.JLabelObject);
handles.JLabelHandle = JLabel('UpdatePrediction',handles.JLabelHandle);
guidata(handles.JLabelObject,handles.JLabelHandle);
end

% --- Executes on button press in pushbutton_ok.
function pushbutton_ok_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pushbutton_apply_Callback(hObject, eventdata, handles)
delete(handles.figure1);
end
