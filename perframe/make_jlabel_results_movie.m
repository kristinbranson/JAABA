function handles = make_jlabel_results_movie(handles,t0,t1,varargin)

% compression: scheme for compression
% outavi_fps: output frames per second
% useVideoWriter: whether to use videowriter class

hselection_visible = get(handles.hselection,'Visible');
set(handles.hselection,'Visible','off');
clipsdir = handles.data.GetFile('clipsdir',handles.expi);
if ~exist(clipsdir,'dir'),
  mkdir(clipsdir);
end
flystr = sprintf('%02d_',handles.flies);
flystr = flystr(1:end-1);
avibasename = sprintf('Target%s_Frames%05dto%05d.avi',flystr,t0,t1);
aviname = fullfile(clipsdir,avibasename);

% open avi file for writing
if handles.useVideoWriter,
  if strcmpi(handles.outavi_compression,'None') || strcmpi(handles.outavi_compression,'Uncompressed AVI'),
    profile = 'Uncompressed AVI';
  else
    profile = handles.outavi_compression;
  end
  aviobj = VideoWriter(aviname,profile);
  set(aviobj,'FrameRate',handles.outavi_fps);
  if ~strcmpi(profile,'Uncompressed AVI'),
    set(aviobj,'Quality',100);
  end
  open(aviobj);
else
  if isempty(handles.avifileTempDataFile),
    aviobj = avifile(aviname,'fps',handles.outavi_fps,'quality',100,'compression',handles.outavi_compression);
  else
    aviobj = myavifile(aviname,'fps',handles.outavi_fps,'quality',100,'compression',handles.outavi_compression,...
      'TempDataFile',avifileTempDataFile);
    fprintf('Temporary data file for avi writing: %s\n',aviobj.TempDataFile);
  end
end

tsave = handles.ts;

axi = 1;
preview_pos = get(handles.panel_previews(axi),'Position');
fig_pos = get(handles.figure_JLabel,'Position');

rect = ...
  [handles.guipos.leftborder_leftpanels,...
  handles.guipos.bottomborder_bottompanels,...
  preview_pos(3),...
  fig_pos(4) - handles.guipos.bottomborder_bottompanels - handles.guipos.topborder_toppanels];

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
  if handles.useVideoWriter,
    writeVideo(aviobj,fr);
  else
    aviobj = addframe(aviobj,fr);
  end
end
handles = JLabel('SetCurrentFrame',handles,axi,tsave,handles.figure_JLabel);

if handles.useVideoWriter,
  close(aviobj);
else
  aviobj = close(aviobj); %#ok<NASGU>
end

for i = 1:numel(handles.hselection),
  set(handles.hselection(i),'Visible',hselection_visible{i});
end

