function handles = make_jlabel_results_movie(handles,t0,t1,varargin)

% compression: scheme for compression
% outavi_fps: output frames per second
% useVideoWriter: whether to use videowriter class

hselection_visible = get(handles.guidata.hselection,'Visible');
set(handles.guidata.hselection,'Visible','off');
clipsdir = handles.guidata.data.GetFile('clipsdir',handles.guidata.expi);
if ~exist(clipsdir,'dir'),
  mkdir(clipsdir);
end
flystr = sprintf('%02d_',handles.guidata.flies);
flystr = flystr(1:end-1);
if iscell(handles.guidata.configparams.behaviors.names),
  behaviorname = sprintf('%s',handles.guidata.configparams.behaviors.names{:});
else
  behaviorname = handles.guidata.configparams.behaviors.names;
end

avibasename = sprintf('%s_Target%s_Frames%05dto%05d.avi',behaviorname,flystr,t0,t1);
aviname = fullfile(clipsdir,avibasename);

% open avi file for writing
if handles.guidata.useVideoWriter,
  if strcmpi(handles.guidata.outavi_compression,'None') || strcmpi(handles.guidata.outavi_compression,'Uncompressed AVI'),
    profile = 'Uncompressed AVI';
  else
    profile = handles.guidata.outavi_compression;
  end
  aviobj = VideoWriter(aviname,profile);
  set(aviobj,'FrameRate',handles.guidata.outavi_fps);
  switch lower(handles.guidata.outavi_compression),
    case 'motion jpeg avi',
      set(aviobj,'Quality',handles.guidata.outavi_quality);
    case {'uncompressed avi','none'},
    case {'archival','motion jpeg 2000'},
      set(aviobj,'CompressionRatio',handles.guidata.outavi_quality);      
  end
  open(aviobj);
else
  if isempty(handles.guidata.avifileTempDataFile),
    aviobj = avifile(aviname,'fps',handles.guidata.outavi_fps,'quality',100,'compression',handles.guidata.outavi_compression);
  else
    aviobj = myavifile(aviname,'fps',handles.guidata.outavi_fps,'quality',100,'compression',handles.guidata.outavi_compression,...
      'TempDataFile',avifileTempDataFile);
    fprintf('Temporary data file for avi writing: %s\n',aviobj.TempDataFile);
  end
end

tsave = handles.guidata.ts;

axi = 1;
preview_pos = get(handles.guidata.panel_previews(axi),'Position');
fig_pos = get(handles.figure_JLabel,'Position');

rect = ...
  [handles.guidata.guipos.leftborder_leftpanels,...
  handles.guidata.guipos.bottomborder_bottompanels,...
  preview_pos(3),...
  fig_pos(4) - handles.guidata.guipos.bottomborder_bottompanels - handles.guidata.guipos.topborder_toppanels];

for t = t0:t1,
  handles = JLabel('SetCurrentFrame',handles,axi,t,handles.figure_JLabel);

  if ~all(get(handles.figure_JLabel,'Position') == fig_pos),
    set(handles.figure_JLabel,'Position',fig_pos);
  end
  fr = getframe(handles.figure_JLabel,rect);
  if t == t0,
    [height,width,~] = size(fr.cdata);
  else
    [height1,width1,~] = size(fr.cdata);
    if height1 ~= height || width1 ~= width,
      fr.cdata = imresize(fr.cdata,[height,width]);
    end
  end
  if handles.guidata.useVideoWriter,
    writeVideo(aviobj,fr);
  else
    aviobj = addframe(aviobj,fr);
  end
end
handles = JLabel('SetCurrentFrame',handles,axi,tsave,handles.figure_JLabel);

if handles.guidata.useVideoWriter,
  close(aviobj);
else
  aviobj = close(aviobj); %#ok<NASGU>
end

for i = 1:numel(handles.guidata.hselection),
  set(handles.guidata.hselection(i),'Visible',hselection_visible{i});
end

