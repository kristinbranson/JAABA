classdef JLabelGUIData < handle
  
  properties (Access=public)
    
    classifierfilename = '';
    configfilename = '';
    defaultpath = '';
    packageoutputdir = '';
    
    isgroundtruthmode = false;
    
    status_bar_text = '';
    idlestatuscolor = [0,1,0];
    busystatuscolor = [1,0,1];

    movie_depth = 1;
    movie_height = 100;
    movie_width = 100;
    movie_filename = [];
    
    configparams = struct;

    panel_previews = [];
    axes_previews = [];
    slider_previews = [];
    edit_framenumbers = [];
    pushbutton_playstops = [];
    axes_timelines = [];
    labels_timelines = [];
    axes_timeline_props = [];
    axes_timeline_labels = [];
    text_timeline_props = [];
    text_timelines = [];
    
    hsplash = [];
    hsplashstatus = [];
    henabled = [];
    
    enabled = true;

    guipos = struct;

    rcfilename = '';
    rc = struct;

    timeline_nframes = 250;
    nframes_jump_go = 30;

    label_shortcuts = [];

    outavi_compression = 'None';
    outavi_fps = 15;
    outavi_quality = 95;
    useVideoWriter = false;

    play_FPS = 2;
    traj_nprev = 25;
    traj_npost = 25;

    bottomAutomatic = 'None';

    needsave = false;
    
    data = [];

    nflies_label = 1;
    classifier = [];

    expi = 0;
    flies = [];
    ts = 0;
    label_state = 0;
    label_imp = [];
    nflies_curr = 0;

    labelcolors = [];
    labelunknowncolor = [0,0,0];
    flies_extra_markersize = 12;
    flies_extra_marker = 'o';
    flies_extra_linestyle = '-';

    scorecolor = [];
    
    correctcolor = [0,.7,0];
    incorrectcolor = [.7,.7,0];
    suggestcolor = [0,.7,.7];
    selection_color = [1,.6,0];
    selection_alpha = .5;
    emphasiscolor = [.7,.7,0];
    unemphasiscolor = [1,1,1];

    togglebutton_label_behaviors = [];
    
    GUIAdvancedMode = false;

    timeline_prop_remove_string = '<html><body><i>Remove</i></body></html>';
    timeline_prop_help_string = '<html><body><i>Help</i></body></html>';
    timeline_prop_options = {};
    
    d = '';
    perframepropis = 1;

    timeline_data_ylims = [];
    
    max_click_dist_preview = .025^2;
    
    preview_zoom_mode = 'center_on_fly';
    zoom_fly_radius = nan(1,2);
    
    menu_view_zoom_options = [];

    selection_t0 = nan;
    selection_t1 = nan;
    selected_ts = nan(1,2);
    buttondown_t0 = nan;
    buttondown_axes = nan;
    selecting = false;
    
    NJObj = [];
    
    hplaying = nan;

    bookmark_windows = [];
    
    plot_labels_manual = true;
    plot_labels_automatic = false;

    doFastUpdates = true;

    axes_preview_curr = 1;
    hslider_listeners = [];
    hflies = [];
    hflies_extra = [];
    hfly_markers = [];
    htrx = [];
    fly_colors = [];
    hlabel_curr = [];
    himage_previews = [];
    hlabels = [];
    hpredicted = [];
    hlabelstarts = [];
    hzoom = [];
    hpan = [];
    himage_timeline_manual = [];
    htimeline_label_curr = [];
    himage_timeline_auto = [];
    htimeline_data = [];
    htimeline_errors = [];
    htimeline_suggestions = [];
    htimeline_gt_suggestions = [];
    hcurr_timelines = [];
    hselection = [];
    
    callbacks = struct;

    readframe = [];
    nframes = nan;
    movie_fid = [];
    movieheaderinfo = struct;
    t0_curr = nan;
    t1_curr = nan;
    labels_plot = struct;
    labels_plot_off = nan;

    current_interval = [];
    didclearselection = false;
    
    open_peripherals = [];

    cache_thread = [];
    cacheSize = 4000;
        
  end
     
  methods (Access=public)
    
    
    % Constructor
    function obj = JLabelGUIData(varargin)
      
    end
    
  end
  
end
