function varargout = JAABASplashScreen(varargin)
% JAABASPLASHSCREEN MATLAB code for JAABASplashScreen.fig
%      JAABASPLASHSCREEN, by itself, creates a new JAABASPLASHSCREEN or raises the existing
%      singleton*.
%
%      H = JAABASPLASHSCREEN returns the handle to a new JAABASPLASHSCREEN or the handle to
%      the existing singleton*.
%
%      JAABASPLASHSCREEN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in JAABASPLASHSCREEN.M with the given input arguments.
%
%      JAABASPLASHSCREEN('Property','Value',...) creates a new JAABASPLASHSCREEN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before JAABASplashScreen_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to JAABASplashScreen_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help JAABASplashScreen

% Last Modified by GUIDE v2.5 13-May-2012 20:22:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @JAABASplashScreen_OpeningFcn, ...
                   'gui_OutputFcn',  @JAABASplashScreen_OutputFcn, ...
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


% --- Executes just before JAABASplashScreen is made visible.
function JAABASplashScreen_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to JAABASplashScreen (see VARARGIN)

% Choose default command line output for JAABASplashScreen
handles.output = hObject;

try
  imname = 'janelia_logo.png';
  versionfile = 'version.txt';
  if exist(imname,'file'),
    im = imread(imname);
  else
    im = 0;
  end
  image(im,'Parent',handles.axes1);
  axis(handles.axes1,'image','off');
  if exist(versionfile,'file'),
    fid = fopen(versionfile,'r');
    while true,
      version = fgetl(fid);
      if ~ischar(version),
        break;
      end
      version = strtrim(version);
      if isempty(version),
        continue;
      end
      break;
    end
  else
    version = '???';
  end

  set(handles.text_info,'String',sprintf('JAABA: The Janelia Automatic Animal Behavior Annotator\nversion %s',version));
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes JAABASplashScreen wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = JAABASplashScreen_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
varargout{2} = handles.text_status;
