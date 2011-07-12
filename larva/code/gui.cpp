#include "common.h"
#include "gui.h"
#include "train.h"

#include <cv.h>
#include <highgui.h>
#include <ml.h>  


#ifdef WIN32
#pragma comment(lib, "cv210")
#pragma comment(lib, "cxcore210")
#pragma comment(lib, "highgui210")
#pragma comment(lib, "ml210")
#endif

#define ZOOM 20
#define DEFAULT_WIDTH 1024
#define DEFAULT_HEIGHT 768
#define CONTROL_WIDTH 260
#define TIMELINE_HEIGHT 40

#define BEHAVIOR_LABEL_WIDTH 60

#ifdef USE_OPENMP
#include <omp.h>
#endif

#define DEBUG_FIT_SKELETON 1

#ifdef __WXMAC__
  #include "ApplicationServices/ApplicationServices.h"
  #define TIMELINE_PANEL_HEIGHT 65
#else
#ifdef WIN32
  #define TIMELINE_PANEL_HEIGHT 75
#else
  #define TIMELINE_PANEL_HEIGHT 45
#endif
#endif


char g_processFile[FILENAME_LENGTH], g_processDir[FILENAME_LENGTH], g_meta_ext[FILENAME_LENGTH], g_detect_meta[FILENAME_LENGTH];
char g_train = 0; // CSC 20110209
//char g_trainFile[FILENAME_LENGTH], g_classifierFile[FILENAME_LENGTH];//CSC: moved to common.cpp

static const wxCmdLineEntryDesc g_cmdLineDesc [] =
{
     { wxCMD_LINE_SWITCH, wxT("h"), wxT("help"), wxT("Displays help on the command line parameters"),
          wxCMD_LINE_VAL_NONE, wxCMD_LINE_OPTION_HELP },
     { wxCMD_LINE_OPTION, wxT("o"), wxT("outline"), wxT("Load outline file"), wxCMD_LINE_VAL_STRING },
     { wxCMD_LINE_OPTION, wxT("d"), wxT("dir"), wxT("Specify the main data folder name (defaults to ../data)"), wxCMD_LINE_VAL_STRING },
     { wxCMD_LINE_SWITCH, wxT("s"), wxT("silent"), wxT("disables the GUI") },
     { wxCMD_LINE_OPTION, wxT("p"), wxT("process"), wxT("Process a .outline file, fitting both the skeleton and behaviors and producing a .anno file"), wxCMD_LINE_VAL_STRING },
     { wxCMD_LINE_OPTION, wxT("P"), wxT("folder"), wxT("Bulk processes all .outline files in a folder, fitting both the skeleton and behaviors to all .outline files"), wxCMD_LINE_VAL_STRING },
     { wxCMD_LINE_OPTION, wxT("m"), wxT("meta"), wxT("Detect and show meta behaviors.  Assumes a config file ../data/behaviors.txt.<str> exists"), wxCMD_LINE_VAL_STRING },
     { wxCMD_LINE_OPTION, wxT("M"), wxT("detect_meta"), wxT("Detect meta behaviors for all files in a directory, assuming a folder has already been processed using -P"), wxCMD_LINE_VAL_STRING },
     { wxCMD_LINE_SWITCH, wxT("l"), wxT("learn"), wxT("learn classifier") },
     { wxCMD_LINE_OPTION, wxT("t"), wxT("train file"), wxT("train file name"), wxCMD_LINE_VAL_STRING },
     { wxCMD_LINE_OPTION, wxT("c"), wxT("classifier file"), wxT("classifier file name"), wxCMD_LINE_VAL_STRING },
     { wxCMD_LINE_NONE }
};


LabelFrame *g_frame = NULL;

IMPLEMENT_APP(LabelApp)

BEGIN_EVENT_TABLE(LabelApp,wxApp)
EVT_KEY_DOWN(LabelApp::OnKeyDown)
END_EVENT_TABLE() 

BEGIN_EVENT_TABLE(LabelFrame, wxFrame)
EVT_CLOSE(LabelFrame::OnClose)

    EVT_BUTTON(Label_Quit, LabelFrame::OnExit)
    EVT_BUTTON(Label_Prev_Frame, LabelFrame::OnPrevFrame)
    EVT_BUTTON(Label_Next_Frame, LabelFrame::OnNextFrame)
    EVT_BUTTON(Label_Play, LabelFrame::OnPlay)
    EVT_BUTTON(Label_Fit_Frame, LabelFrame::OnFitFrame)
    EVT_BUTTON(Label_Fit_All, LabelFrame::OnFitAll)
    EVT_BUTTON(Label_Fit_Behaviors, LabelFrame::OnFitBehaviors)
    EVT_BUTTON(Label_Fit_Skeleton, LabelFrame::OnFitSkeleton)
    EVT_BUTTON(Label_Flip, LabelFrame::OnFlip)
    EVT_BUTTON(Label_Save, LabelFrame::OnSave)
    EVT_BUTTON(Label_Add, LabelFrame::OnAdd)
    EVT_BUTTON(Label_Delete, LabelFrame::OnDelete)
    EVT_BUTTON(Label_Validate, LabelFrame::OnValidate)
    EVT_BUTTON(Label_Train, LabelFrame::OnTrain)
    EVT_BUTTON(Label_Backup, LabelFrame::OnBackup)
    EVT_BUTTON(Label_Zoom_Plus, LabelFrame::OnZoomPlus)
    EVT_BUTTON(Label_Zoom_Minus, LabelFrame::OnZoomMinus)
   
    EVT_LISTBOX(Label_Train_List, LabelFrame::OnTrainList)

    EVT_MENU(wxID_EXIT, LabelFrame::OnExit)
    EVT_MENU(wxID_OPEN, LabelFrame::OnOpen)
    EVT_MENU(wxID_SAVE, LabelFrame::OnSave)
    EVT_MENU(wxID_SAVEAS, LabelFrame::OnSaveAs)

    EVT_SIZE(LabelFrame::OnSize)

//EVT_KEY_DOWN(LabelFrame::OnKeyDown)

    EVT_TIMER(Label_Play_Timer, LabelFrame::OnPlayTimer)

END_EVENT_TABLE()

BEGIN_EVENT_TABLE(BehaviorDialog, wxDialog)
    EVT_BUTTON(Behavior_Ok, BehaviorDialog::OnOk)
    EVT_BUTTON(Behavior_Cancel, BehaviorDialog::OnCancel)
END_EVENT_TABLE()

BEGIN_EVENT_TABLE(ImageWindow,wxScrolledWindow)
    EVT_PAINT(ImageWindow::OnPaint)
	EVT_ERASE_BACKGROUND(ImageWindow::OnEraseBackground)
END_EVENT_TABLE()

BEGIN_EVENT_TABLE(BlobWindow,wxScrolledWindow)
    EVT_PAINT(BlobWindow::OnPaint)
    EVT_MOUSE_EVENTS(BlobWindow::OnMouse)
    EVT_SIZE(BlobWindow::OnSize)
END_EVENT_TABLE()

BEGIN_EVENT_TABLE(TimelineWindow,wxScrolledWindow)
    EVT_PAINT(TimelineWindow::OnPaint)
    EVT_MOUSE_EVENTS(TimelineWindow::OnMouse)
    EVT_SIZE(TimelineWindow::OnSize)
END_EVENT_TABLE()


BEGIN_EVENT_TABLE(GraphWindow,wxScrolledWindow)
    EVT_PAINT(GraphWindow::OnPaint)
    EVT_MOUSE_EVENTS(GraphWindow::OnMouse)
    EVT_SIZE(GraphWindow::OnSize)
END_EVENT_TABLE()




// Temporary Hack: these functions should be in blob.cpp, but instead have put it here because it uses wx
#include "wx/regex.h"

const char *g_charMap = "abcdefghijklmnopqrstuvqxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789";
int g_charMapLen = strlen(g_charMap);

char *BehaviorCharString(BlobSequence *b, BehaviorGroups *beheviors, int beh) {
  char *str = (char*)malloc(1);
  unsigned int len = 0, alloc = 0;
  char tmp[400];

  // Build a string representation of this behavior sequence
  str[0] = '\0';
  for(int i = 0; i < b->num_frames; i++) {
    //if(i && b->frames[i].behaviors[beh] == b->frames[i-1].behaviors[beh])
    //continue;
    assert(b->frames[i].behaviors[beh] < g_charMapLen);
    sprintf(tmp, "%c", g_charMap[b->frames[i].behaviors[beh]]);
    if(len+strlen(tmp)+1 >= alloc) {
      alloc = (int)(alloc*1.1 + 1024);
      str = (char*)realloc(str, alloc);
    }
    strcat(str+len, tmp);
    len += strlen(tmp);
  }

  return str;
}

void DetectRegularExpression(const char *regExp, char *str, BlobSequence *b, int beh, int beh_val) {
  // Find all matches for this regular expression, append the corresponding behaviors to the list of detected meta_behaviors,
  // then remove it from str
  wxRegEx reg(wxString(regExp, wxConvUTF8));
  size_t start, len;
  while(1) {
    wxString wstr = wxString(str, wxConvUTF8);
    if(reg.Matches(wstr)) {
      reg.GetMatch(&start, &len); 
      for(int k = start; k < (int)(start+len); k++) {
	if(!(start > 0 && k == (int)start && b->frames[start-1].meta_behaviors[beh] == beh_val))
	  b->frames[k].meta_behaviors[beh] = beh_val;
	str[k] = '@';
      }
    } else
      break;
  }
}

void RegularExpressionCharString(const char *regExp, char *regA, BehaviorGroup *behaviors) {
  // Convert the regular expression for this behavior, into a string representation
  int inBracks = 0, num = 0, tagInd=0;
  char tag[1000];
  for(int i = 0; i < (int)strlen(regExp); i++) {
    if(regExp[i] == '<') {
      if(inBracks) {
	fprintf(stderr, "ERROR: nested '<' tags in %s\n", regExp);
	return;
      }
      inBracks = 1;
      tagInd = 0;
    } else if(regExp[i] == '>') {
      if(!inBracks) {
	fprintf(stderr, "ERROR: '>' not proceeded by a start tag in %s\n", regExp);
	return;
      }
      inBracks = 0;
      int ind = -1;
      tag[tagInd] = '\0';
      for(int j = 0; j < behaviors->num_values; j++) {
	if(tagInd == (int)strlen(behaviors->values[j].abbreviation) &&
	   !strncmp(behaviors->values[j].abbreviation, tag, tagInd)) {
	  ind = j;
	  break;
	}
      }
      tagInd = 0;
      if(ind < 0) {
	fprintf(stderr, "ERROR: tag %s not found in %s\n", tag, regExp);
	return;
      }
      regA[num++] = g_charMap[ind];
    } else if(!inBracks) {
      regA[num++] = regExp[i];
    } else
      tag[tagInd++] = regExp[i];
  }
  regA[num++] = '\0';
}


void ExtractMetaBehaviors(BlobSequence *b, BehaviorGroups *behaviors, BehaviorGroups *meta_behaviors) {
  if(meta_behaviors->num != behaviors->num) {
    fprintf(stderr, "ERROR: meta behaviors must have the same number of groups as original behaviors\n");
    return;
  }

  // Convert behaviors to a sequence of meta behaviors by parsing regular expressions
  char regA[1000];
  for(int beh = 0; beh < meta_behaviors->num; beh++) {
    char *str = BehaviorCharString(b, behaviors, beh);

    for(int k = 0; k < b->num_frames; k++) 
      b->frames[k].meta_behaviors[beh] = 0;
    for(int e = 2; e < meta_behaviors->behaviors[beh].num_values; e++) {
      RegularExpressionCharString(meta_behaviors->behaviors[beh].values[e].regExp, regA, &behaviors->behaviors[beh]);
      DetectRegularExpression(regA, str, b, beh, e);
    }
    free(str);
  }
}

// ============================================================================
// implementation
// ============================================================================

// ----------------------------------------------------------------------------
// app class
// ----------------------------------------------------------------------------

bool LabelApp::OnInit()
{
#ifdef __WXMAC__
    ProcessSerialNumber PSN;
    GetCurrentProcess(&PSN);
    TransformProcessType(&PSN,kProcessTransformToForegroundApplication);
#endif

    if ( !wxApp::OnInit() )
        return false;

    m_frame = new LabelFrame(_T("Larva Label Tool"), fname.length() ? (const char*)fname.mb_str() : NULL);
    m_frame->Show();

    //wxLog::AddTraceMask(_T("listbox"));
    //wxLog::AddTraceMask(_T("scrollbar"));
    //wxLog::AddTraceMask(_T("focus"));

    return true;
}

bool LabelApp::OnCmdLineParsed(wxCmdLineParser& parser)
{
  fname = wxT("");
  parser.Found(wxT("o"), &fname);

  wxString dir, fname, meta;


  strcpy(g_trainFile, "train.txt"); // default value
  strcpy(g_classifierFile, "classifier"); // default value


  if(parser.Found(wxT("d"), &dir)) {
    char str[1000];
    strcpy(str, wxString(dir).mb_str());
    DATA_DIR = StringCopy(str);
    fprintf(stderr, "Data dir is %s\n", str);
  }
  if(parser.Found(wxT("l"))) {
    g_train = 1;
    fprintf(stderr, "learn request detected\n");
  }
  strcpy(g_processFile, "");
  strcpy(g_processDir, "");
  strcpy(g_meta_ext, "");
  strcpy(g_detect_meta, "");
  if(parser.Found(wxT("p"), &fname)) {
    strcpy(g_processFile, wxString(fname).mb_str());
  }
  if(parser.Found(wxT("P"), &dir)) {
    strcpy(g_processDir, wxString(dir).mb_str());
  }
  if(parser.Found(wxT("m"), &meta)) {
    strcpy(g_meta_ext, wxString(meta).mb_str());
  }
  if(parser.Found(wxT("M"), &meta)) {
    strcpy(g_detect_meta, wxString(meta).mb_str());
    if(!strstr(g_detect_meta, "dir.multi"))
      strcat(g_detect_meta, "/dir.multi");
  }

  if(parser.Found(wxT("t"), &fname)) {
    strcpy(g_trainFile, wxString(fname).mb_str());
    fprintf(stderr, "train filename is %s\n", g_trainFile);
  }
  if(parser.Found(wxT("c"), &dir)) {
    strcpy(g_classifierFile, wxString(dir).mb_str());
    fprintf(stderr, "classifier directory is %s\n", g_classifierFile);
  }


  return true;
}

void LabelApp::OnInitCmdLine(wxCmdLineParser& parser)
{
    parser.SetDesc(g_cmdLineDesc);
    parser.SetSwitchChars(wxT("-"));
}


// ----------------------------------------------------------------------------
// LabelFrame construction
// ----------------------------------------------------------------------------

void LabelFrame::LoadBehaviors() {
  char bname[1000], classifier_dir[300];
  if(behaviors) free_behaviors(behaviors);
  sprintf(bname, "%s/behaviors.txt", DATA_DIR);
  sprintf(classifier_dir, "%s/%s", DATA_DIR, g_classifierFile); // CSC: add (by default empty) classifier appendix, for cross-validation classifier dir's are named classifier.<idx>, where idx is the index of the left out training resp. test file.
  behaviors = load_behaviors(bname, classifier_dir);
  assert(behaviors);

  if(strlen(g_meta_ext)) {
    sprintf(bname, "%s/behaviors.txt.%s", DATA_DIR, g_meta_ext);
    meta_behaviors = load_behaviors(bname, classifier_dir);
  }
}

LabelFrame::LabelFrame(const wxString& title, const char *fname)
            : wxFrame(NULL, wxID_ANY, title,
                      wxPoint(0, 50), wxSize(DEFAULT_WIDTH,DEFAULT_HEIGHT),
                      wxDEFAULT_FRAME_STYLE |
                      wxNO_FULL_REPAINT_ON_RESIZE |
                      wxCLIP_CHILDREN |
                      wxTAB_TRAVERSAL)
{
  int i;

  g_frame = this;
  timeline_zoom = 1;

    blobs = NULL;
    preview_frame = frame = 0;
    play = 0;
    mutex = new wxMutex();
    dirty = false;
    sel_train_file = -1;
    train_list_strings = NULL;
    sel_timeline = -1;
    timer = NULL;
    behaviors = NULL;
    meta_behaviors = NULL;

    /*
#if wxUSE_MENUS
    // create the menubar
    wxMenuBar *mbar = new wxMenuBar;
    wxMenu *menuWidget = new wxMenu;
#if wxUSE_TOOLTIPS
#endif // wxUSE_TOOLTIPS
    menuWidget->Append(wxID_OPEN, _T("&Open\tCtrl-O"));
    menuWidget->Append(wxID_SAVE, _T("&Save\tCtrl-S"));
    menuWidget->Append(wxID_SAVEAS, _T("&Save As"));
    menuWidget->Append(wxID_EXIT, _T("&Quit\tCtrl-Q"));
    mbar->Append(menuWidget, _T("&File"));
    //SetMenuBar(mbar);
#endif // wxUSE_MENUS
    */

    params = default_parameters();

//fprintf(stderr, "test\n");
    LoadBehaviors();
    
    char bname[400];
    sprintf(bname, "%s/%s", DATA_DIR, g_trainFile);
    train_list = load_train_list(bname, &num_train_files);
//fprintf(stderr, "bname = '%s', num_train_files = %d\n", bname, num_train_files);
    GenerateTrainListStrings();
    
    /*
    for(int j = 0; j < num_train_files; j++) {
      BlobSequence *b = load_blob_sequence(train_list[j], behaviors);
      int k, has_frames;
      for(k = 0; k < b->num_frames; k++) {
	if((b->frames[k].is_manual & 3) == 3 && b->frames[k].num_model_pts == params.num_spine_lines+1) {
	  has_frames = 1;
	  break;
	}
      }
      if(!has_frames) {
	free_blob_sequence(b);
	continue;
      }

      for(k = 0; k < b->num_frames; k++) {
	if(b->frames[k].behaviors[0] > 1 && b->frames[k].behaviors[0] < 12)
	  b->frames[k].behaviors[0] += 11;
	else if(b->frames[k].behaviors[0] > 12)
	  b->frames[k].behaviors[0] -= 11;
      }
      save_blob_sequence(b->fname, b, behaviors, params.num_spine_lines+1, params.num_orientations);
    }
    */

    // create controls
    timeline_splitter = new wxSplitterWindow(this, -1, wxDefaultPosition, wxDefaultSize, wxSP_3D);
    img_splitter = new wxSplitterWindow(timeline_splitter, -1, wxDefaultPosition, 
					wxDefaultSize, wxSP_3D);
    m_timeline_panel = new wxPanel(timeline_splitter, wxID_ANY, wxDefaultPosition, 
				   wxDefaultSize, wxCLIP_CHILDREN); 
    img_window = new BlobWindow(img_splitter, this);
    char rname[400];
    sprintf(rname, "%s/rotate.txt", DATA_DIR);
    FILE *fin = fopen(rname, "r");
    if(fin) {
      int rot, x, y;
      if(fscanf(fin, "rotate=%d flip_x=%d flip_y=%d", &rot, &x, &y))
	img_window->SetRotate(rot, x, y);
      fclose(fin);
    }

    m_control_panel = new wxPanel(img_splitter, wxID_ANY, wxDefaultPosition, 
				  wxDefaultSize, wxCLIP_CHILDREN);   

    m_timeline_sizer = new wxFlexGridSizer(2,1,0,0); 
    m_timeline_sizer2 = new wxFlexGridSizer((meta_behaviors ? 2 : 1)*behaviors->num+1,3,0,0);   
    for(i = 0; i < behaviors->num; i++) { 
      m_behavior_timelines[i] = new BehaviorTimelineControls(this, m_timeline_panel, i, m_timeline_sizer2, i);
    }
      

    m_graph = new GraphControls(this, m_timeline_panel, m_timeline_sizer2, behaviors->num);
    m_timeline_sizer->Add(m_timeline_sizer2, i, wxALIGN_TOP);

    if(meta_behaviors) {
      for(i = 0; i < meta_behaviors->num; i++) { 
	m_meta_behavior_timelines[i] = new BehaviorTimelineControls(this, m_timeline_panel, i, m_timeline_sizer2, i, true);
      }
    }

    m_frame_sizer = new wxFlexGridSizer(1,16,0,0);
    m_frame_label = new wxStaticText(m_timeline_panel, wxID_ANY, wxT("Frame: "));
    m_frame_number_label = new wxStaticText(m_timeline_panel, wxID_ANY, wxT("    "));
    m_frame_sizer->Add(m_frame_label, 0, wxLEFT | wxALIGN_CENTER_VERTICAL);
    m_frame_sizer->Add(m_frame_number_label, 1, wxLEFT | wxALIGN_CENTER_VERTICAL);
    m_frame_sizer->AddSpacer(30);

    m_prev_frame_button = new wxButton(m_timeline_panel, Label_Prev_Frame, wxT("Prev Frame"));
    m_play_button = new wxButton(m_timeline_panel, Label_Play, wxT("Play"));
    m_next_frame_button = new wxButton(m_timeline_panel, Label_Next_Frame, wxT("Next Frame"));
    m_frame_sizer->Add(m_prev_frame_button, 3, wxLEFT | wxALIGN_CENTER_VERTICAL);
    m_frame_sizer->Add(m_play_button, 4, wxLEFT | wxALIGN_CENTER_VERTICAL);
    m_frame_sizer->Add(m_next_frame_button, 5, wxLEFT | wxALIGN_CENTER_VERTICAL);
    m_frame_sizer->AddSpacer(30);

    //m_zoom_plus_button = new wxButton(m_timeline_panel, Label_Zoom_Plus, wxT("Zoom In"));
    //m_zoom_minus_button = new wxButton(m_timeline_panel, Label_Zoom_Minus, wxT("Zoom Out"));
    //m_frame_sizer->Add(m_zoom_plus_button, 7, wxLEFT | wxALIGN_CENTER_VERTICAL);
    //m_frame_sizer->Add(m_zoom_minus_button, 8, wxLEFT | wxALIGN_CENTER_VERTICAL);
    //m_frame_sizer->AddSpacer(30);

    m_fit_frame_button = new wxButton(m_timeline_panel, Label_Fit_Frame, wxT("Fit Frame"));
    m_fit_all_button = new wxButton(m_timeline_panel, Label_Fit_All, wxT("Fit All"));
    m_fit_skeleton_button = new wxButton(m_timeline_panel, Label_Fit_Skeleton, wxT("Fit Skeleton"));
    m_fit_behaviors_button = new wxButton(m_timeline_panel, Label_Fit_Behaviors, wxT("Fit Behaviors"));
    m_flip_button = new wxButton(m_timeline_panel, Label_Flip, wxT("Flip^"));
    m_validate_button = new wxButton(m_timeline_panel, Label_Validate, wxT("Validate^"));
    m_frame_sizer->Add(m_fit_frame_button, 10, wxLEFT | wxALIGN_CENTER_VERTICAL);
    m_frame_sizer->Add(m_fit_skeleton_button, 11, wxLEFT | wxALIGN_CENTER_VERTICAL);
    m_frame_sizer->Add(m_fit_all_button, 12, wxLEFT | wxALIGN_CENTER_VERTICAL);
    m_frame_sizer->Add(m_fit_behaviors_button, 13, wxLEFT | wxALIGN_CENTER_VERTICAL);
    m_frame_sizer->Add(m_flip_button, 14, wxLEFT | wxALIGN_CENTER_VERTICAL);
    m_frame_sizer->Add(m_validate_button, 15, wxLEFT | wxALIGN_CENTER_VERTICAL);
    m_timeline_sizer->Add(m_frame_sizer, behaviors->num, wxALIGN_TOP);
    
    


    m_control_sizer = new wxFlexGridSizer(5,1,0,0);

    m_control_sizer->AddSpacer(15);
    m_train_label = new wxStaticText(m_control_panel, wxID_ANY, wxT("Training Files:"));
    m_train_list = new wxListBox(m_control_panel, Label_Train_List, wxDefaultPosition, 
				 wxSize(CONTROL_WIDTH-10,-1), num_train_files, train_list_strings, wxLB_SINGLE);
    m_control_sizer->Add(m_train_label, 1, wxCENTER | wxALIGN_TOP);
    m_control_sizer->Add(m_train_list, 2, wxCENTER | wxALIGN_TOP);

    m_file_sizer = new wxFlexGridSizer(2,3,0,0);
    m_file_label = new wxStaticText(m_control_panel, wxID_ANY, wxT(""));
    m_add_button = new wxButton(m_control_panel, Label_Add, wxT("Add"));
    m_delete_button = new wxButton(m_control_panel, Label_Delete, wxT("Delete"));
    m_save_button = new wxButton(m_control_panel, Label_Save, wxT("Save"));
    m_file_sizer->Add(m_add_button, 0, wxLEFT | wxALIGN_TOP);
    m_file_sizer->Add(m_delete_button, 1, wxLEFT | wxALIGN_TOP);
    m_file_sizer->Add(m_save_button, 2, wxLEFT | wxALIGN_TOP);

    m_train_button = new wxButton(m_control_panel, Label_Train, wxT("Train"));
    m_backup_button = new wxButton(m_control_panel, Label_Backup, wxT("Backup"));
    m_control_sizer->Add(m_file_label, 3, wxLEFT | wxALIGN_TOP);
    m_control_sizer->Add(m_file_sizer, 4, wxLEFT | wxALIGN_TOP);
    m_file_sizer->Add(m_train_button, 5, wxLEFT | wxALIGN_TOP);
    m_file_sizer->Add(m_backup_button, 6, wxLEFT | wxALIGN_TOP);

    m_control_panel->Show(true);
    m_control_panel->SetSizerAndFit(m_control_sizer);
    m_timeline_panel->Show(true);
    m_timeline_panel->SetSizerAndFit(m_timeline_sizer);
    img_window->Show(true);
    img_splitter->SplitVertically(img_window, m_control_panel, DEFAULT_WIDTH-CONTROL_WIDTH);
    img_splitter->Show(true);
    timeline_splitter->SplitHorizontally(img_splitter, m_timeline_panel, DEFAULT_HEIGHT-TIMELINE_PANEL_HEIGHT-TIMELINE_HEIGHT*(behaviors->num+1));
    timeline_splitter->SetSashPosition(DEFAULT_HEIGHT-TIMELINE_PANEL_HEIGHT-TIMELINE_HEIGHT*(behaviors->num+1));
    img_splitter->SetSashPosition(DEFAULT_WIDTH-CONTROL_WIDTH);

    if(fname) 
      Open(fname);

    if(g_train) {
      fprintf(stderr, "training new classifier - this may take several hours...\n");
	wxCommandEvent event;
      OnTrain(event); // pass an WXUNUSED(event)
      fprintf(stderr, "All done!\n");
      exit(0);
    }
//    if(g_crossValidation) {
//      fprintf(stderr, "performing cross-validation analysis - this may take several dozens of hours...\n");
//      crossValidate();
//      fprintf(stderr, "All done!\n");
//      exit(0);      
//    }

    if(strlen(g_processFile)) {
      ProcessFile(g_processFile);
      fprintf(stderr, "All done!\n");
    }
    if(strlen(g_processDir)) {
      ProcessDirectory(g_processDir);
      fprintf(stderr, "All done!\n");
      exit(0);
    }
    if(strlen(g_detect_meta)) {
      DetectMetaBehaviors(g_detect_meta);
      fprintf(stderr, "All done!\n");
      exit(0);
    }

    /*for(int k = 0; k < num_train_files; k++) {
      Open(train_list[k]);
      Fit();
      Save(train_list[k]);
      }*/

    SynchEnable();
}

LabelFrame::~LabelFrame()
{
}

void LabelFrame::crossValidate() {// CSC
  char classifier_dir[FILENAME_LENGTH];
	char *curr_file;

// preferably do in parallel - fork?:

  for(int i = 0; i < num_train_files; i++) {

	// training step

//  	sprintf(classifier_dir, "%s/%s", DATA_DIR, g_classifierFile);
strcpy(classifier_dir, g_classifierFile);
	curr_file = train_list[i];
	train_list[i] = 0;
  	train_behavior_detector(classifier_dir, train_list, num_train_files, behaviors, &params);
	train_list[i] = curr_file;

	// prediction step

	char classifier_appendix[10]; // allows up to 99.999.999 test files
	sprintf(classifier_appendix, ".%d", i);
	LoadBehaviors();

	char orig_name[FILENAME_LENGTH], test_name[FILENAME_LENGTH], result_name[FILENAME_LENGTH];
	sprintf(orig_name, "%s.orig", curr_file);
	sprintf(test_name, "%s.test", curr_file);
	sprintf(result_name, "%s.result", curr_file);

	rename(curr_file, orig_name);
	ProcessFile(curr_file); // the rename steps make sure that we don't overwrite the manual annotation.
	rename(curr_file, test_name);
	rename(orig_name, curr_file);

	// compare prediction vs. manual annotation

	// execute compiled matlab script 'compareAnnoFiles curr_file test_name Turn > result_name';
  }
}



// ----------------------------------------------------------------------------
// LabelFrame event handlers
// ----------------------------------------------------------------------------

void LabelFrame::OnClose(wxCloseEvent& WXUNUSED(event)) {
  wxCommandEvent c;
  //this->Destroy();
  OnExit(c);
}

void LabelFrame::OnExit(wxCommandEvent& WXUNUSED(event))
{
  Cleanup();
  exit(0);

  Close();
}

void LabelFrame::Cleanup() {
  for(int i = 0; i < behaviors->num; i++)
    delete m_behavior_timelines[i];
  if(meta_behaviors)
    for(int i = 0; i < meta_behaviors->num; i++)
      delete m_meta_behavior_timelines[i];
  delete m_graph;

  if(blobs) free_blob_sequence(blobs);
  if(behaviors) free_behaviors(behaviors);
  if(train_list) free_train_list(train_list, num_train_files);

  if(img_window) img_window->Cleanup();
  //if(train_list_strings) delete [] train_list_strings;

  if(timer) delete timer;
  if(mutex) delete mutex;
}


void LabelFrame::OnOpen(wxCommandEvent& WXUNUSED(event)) {
  if(!PromptSave())
    return;

  
  wxFileDialog fd(this, wxT("Open a Blob File"), wxString(DATA_DIR, wxConvUTF8), wxT(""), 
		  wxT("Outline files (*.outline)|*.outline|Annotation files (*.anno)|*.anno"));
  if(fd.ShowModal() == wxID_OK) {
    wxString path = fd.GetPath();
    Open(path.mb_str());
  }
}


void LabelFrame::GenerateTrainListStrings() { 
  if(train_list_strings) delete [] train_list_strings;

  train_list_strings = new wxString[num_train_files];
  int day, sec, num;
  char str[400];

  for(int i = 0; i < num_train_files; i++) {
    strcpy(str, train_list[i]);
    if(!strncmp(train_list[i], DATA_DIR, strlen(DATA_DIR))) {
      if(sscanf(str+strlen(DATA_DIR)+1, "%d_%d", &day, &sec) && 
	 sscanf(str+strlen(str)-10, "%d", &num))
	sprintf(str, "%08d_%06d %05d", day, sec, num);
    }
      
    train_list_strings[i] = wxString(str, wxConvUTF8);
  }
}


void LabelFrame::OnAdd(wxCommandEvent& event) {
  if(!PromptSave())
    return;

  
  char *ptr;
  wxFileDialog fd(this, wxT("Open a Blob File"), wxString(DATA_DIR, wxConvUTF8), wxT(""), 
		  wxT("Outline files (*.outline)|*.outline|Annotation files (*.anno)|*.anno"));
  if(fd.ShowModal() == wxID_OK) {
    wxString path = fd.GetPath();
    
    train_list = (char**)realloc(train_list, sizeof(char*)*(num_train_files+1));
    train_list[num_train_files] = (char*)malloc(sizeof(char)*(strlen(path.mb_str())+100));
    strcpy(train_list[num_train_files], path.mb_str());
    if((ptr = strstr(train_list[num_train_files], ".outline")) != NULL)
      strcpy(ptr, ".anno");
    sel_train_file = num_train_files++;
    GenerateTrainListStrings();
    m_train_list->InsertItems(1, &train_list_strings[num_train_files-1], num_train_files-1);
    
    char bname[400];
    sprintf(bname, "%s/train.txt", DATA_DIR);
    assert(save_train_list(bname, train_list, num_train_files));

    Open(path.mb_str());
    if(!FileExists(train_list[num_train_files-1])) {
      Fit();
      Save(blobs->fname);
    } 
    

    SynchFileName();
  }
}

void LabelFrame::OnDelete(wxCommandEvent& event) {
  int sel = m_train_list->GetSelection();
  if(sel < 0 || sel >= num_train_files) return;
  char fname[1000];
  strcpy(fname, train_list[sel]);
  //wxString str = wxString::Format(_T("Are you sure remove %s from the list of training files?"), fname);
  wxString str = wxString(_T("Are you sure remove ")).Append(wxString(fname,wxConvUTF8)).Append(_T(" from the list of training files?"));
  wxMessageDialog dlg(this, str, wxT("Remove File?"), wxYES_NO);
  int retval = dlg.ShowModal();
  wxCommandEvent evt;
  if(retval == wxID_YES) {
    free(train_list[sel]);
    for(int i = sel; i < num_train_files-1; i++)
      train_list[i] = train_list[i+1];
    num_train_files--;
    
    GenerateTrainListStrings();
    m_train_list->Delete(sel);

    char bname[400];
    sprintf(bname, "%s/train.txt", DATA_DIR);
    assert(save_train_list(bname, train_list, num_train_files));
  }
}

void LabelFrame::OnSave(wxCommandEvent& event) {
  if(strlen(blobs->fname))
    Save(blobs->fname);
  else
    OnSaveAs(event);
}

void LabelFrame::OnSaveAs(wxCommandEvent& WXUNUSED(event)) {
  wxFileDialog fd(this, wxT("Save a Blob Annotation File"), wxT("DATA_DIR"), wxT(""), 
		  wxT("Annotation files (*.anno)|*.anno"));
  if(fd.ShowModal() == wxID_OK) {
    wxString path = fd.GetPath();
    Save(path.mb_str());
  }
}


void LabelFrame::Save(const char *fname) {
  mutex->Lock();
  char fnameB[400];
  save_blob_sequence(fname, blobs, behaviors, params.num_spine_lines+1, params.num_orientations);
  int has_b = 0;
  for(int i = 0; i < blobs->num_frames; i++) {
    if(blobs->frames[i].behaviors[0])
      has_b = 1;
  }
  if(has_b) {
    sprintf(fnameB, "%s.bouts", fname);
    save_blob_behavior_bouts(fnameB, blobs, behaviors, 0);
    if(meta_behaviors) {
      sprintf(fnameB, "%s.bouts.%s", fname, g_meta_ext);
      save_blob_behavior_bouts(fnameB, blobs, meta_behaviors, 1);
    }
  }
  mutex->Unlock();
  SetDirty(false);
}

 
bool LabelFrame::PromptSave() {
  if(dirty) {
    wxString str = wxString::Format(_T("Save changes?"));
    wxMessageDialog dlg(this, str, wxT("Save changes?"), wxYES_NO);
    int retval = dlg.ShowModal();
    wxCommandEvent evt;
    if(retval == wxID_YES)
      OnSave(evt);
    //else if(retval == wxID_NO)
      //  return false;
  }
  return true;
}



void LabelFrame::Open(const char *fname) {
  int is_import=0;

  BlobSequence *b;
  if(strstr(fname, ".anno")) {
    b = load_blob_sequence(fname, behaviors);
    is_import = 1;
  } else {
    b = import_blob_sequence(fname, /*params.num_spine_lines+1*/0);
    is_import = 1;
  }

  assert(b && b->num_frames > 0);
  SetBlobs(b, is_import);
}

void LabelFrame::SetBlobs(BlobSequence *b, int is_import) {
  int i;

  mutex->Lock();
  if(blobs && blobs != b) {
    free_blob_sequence(blobs);
    img_window->ClearBlobs();
    for(i = 0; i < behaviors->num; i++)
      m_behavior_timelines[i]->Clear();
    if(meta_behaviors)
      for(i = 0; i < meta_behaviors->num; i++)
	m_meta_behavior_timelines[i]->Clear();
  }

  blobs = b;
  params.world_to_pixel_scale = guess_world_to_pixel_scale(blobs->frames[0].contour, 
							   blobs->frames[0].num_pts);
  
  for(i = 0; i < blobs->num_frames; i++) {
    if(blobs->frames[i].num_model_pts != params.num_spine_lines+1 && (blobs->frames[i].is_manual&1))
      blobs->frames[i].is_manual -= 1;
    if(!blobs->frames[i].features)
      blobs->frames[i].features = (double*)malloc(2*sizeof(double)*((2*F_NUM_FEATURES)*(params.num_spine_lines+1)+F_NUM_GLOBAL_FEATURES));
  }
  if(meta_behaviors)
    ExtractMetaBehaviors(b, behaviors, meta_behaviors);


  compute_deterministic_spine_attributes_blob_sequence(blobs, &params, is_import);
  compute_all_global_features(blobs, &params, 0);

  frame = preview_frame = 0;
  img_window->AddBlob(&blobs->frames[0]);
  img_window->CropToSequence(blobs);
  mutex->Unlock();
  SetTimelineZoom(1);
  SetFrame(0);
  img_window->Refresh(true);
  OnFeatureCombo();
  SetDirty(false);
  SynchEnable();
}

void LabelFrame::OnTrainList(wxCommandEvent& WXUNUSED(event)) {
  wxArrayInt a;
  int num = m_train_list->GetSelections(a);
  if(num) {
    if(!PromptSave())
      return;
    sel_train_file = a[0];
    Open(train_list[a[0]]);
    SynchFileName();
  }
}

void LabelFrame::OnTrain(wxCommandEvent& WXUNUSED(event)) {
  char classifier_dir[400];

  sprintf(classifier_dir, "%s/%s", DATA_DIR, g_classifierFile);
  train_behavior_detector(classifier_dir, train_list, num_train_files, behaviors, &params);
}

void LabelFrame::SynchFileName() {
  m_file_label->SetLabel((dirty ? wxT("*") : wxT("")) + (sel_train_file >= 0 ? train_list_strings[sel_train_file] : wxT("")));
}

void LabelFrame::SynchEnable() {
  bool busy = IsBusy();

  m_save_button->Enable(dirty && !busy && blobs);
  m_add_button->Enable(!busy);
  
  m_next_frame_button->Enable(!busy && blobs);
  m_play_button->Enable((!busy || play) && blobs);
  m_prev_frame_button->Enable(!busy && blobs);
  m_validate_button->Enable(!busy && blobs);

  m_fit_frame_button->Enable(!busy && blobs);
  m_fit_all_button->Enable(!busy && blobs);
  m_fit_behaviors_button->Enable(!busy && blobs);
  m_fit_skeleton_button->Enable(!busy && blobs);
  m_flip_button->Enable(!busy && blobs);

  for(int i = 0; i < behaviors->num; i++) 
    m_behavior_timelines[i]->SynchEnable(busy, blobs);

  if(meta_behaviors)
    for(int i = 0; i < meta_behaviors->num; i++) 
      m_meta_behavior_timelines[i]->SynchEnable(busy, blobs);

  m_train_list->Enable(!busy);
}


void LabelFrame::OnPrevFrame(wxCommandEvent& WXUNUSED(event)) {
  if(frame > 0)
    SetFrame(frame-1);
}

void LabelFrame::OnNextFrame(wxCommandEvent& WXUNUSED(event)) {
  if(frame < blobs->num_frames-1)
    SetFrame(frame+1);
}



void LabelFrame::OnPlay(wxCommandEvent& WXUNUSED(event)) {
  if(!play) { 
    mutex->Lock();
    for(int i = 0; i < behaviors->num; i++)
      m_behavior_timelines[i]->Clear();

    if(meta_behaviors)
      for(int i = 0; i < meta_behaviors->num; i++)
	m_meta_behavior_timelines[i]->Clear();

    play = 1;
    SynchEnable();
    start_frame = frame;
    mutex->Unlock();
    timer = new wxTimer(this, Label_Play_Timer);
    timer->Start((blobs->frames[blobs->num_frames-1].frame_time -
		  blobs->frames[frame].frame_time)*1000 / 
		   (blobs->num_frames-frame));
  } else {
    mutex->Lock();
    play = 0;
    mutex->Unlock();
  }
  SynchEnable();
}

void LabelFrame::Backup(const char *dir) {
  char fname[10000], tmp[10000];
  wxString dir_str = wxString(dir, wxConvUTF8);
  char **tlist = (char**)malloc(sizeof(char*)*num_train_files);
  int i, j;

  assert(::wxMkdir(dir_str));

  // Copy trained classifier files
  sprintf(fname, "%s/%s", dir, g_classifierFile);
  assert(::wxMkdir(wxString(fname, wxConvUTF8)));
  if(behaviors->classifier_method) {
    sprintf(fname, "%s/%s/All.%s", DATA_DIR, g_classifierFile, g_classifer_extensions[behaviors->classifier_method]);
    if(FileExists(fname)) {
      sprintf(tmp, "%s/%s/All.%s", dir, g_classifierFile, g_classifer_extensions[behaviors->classifier_method]);
      fprintf(stderr, "Copying %s...\n", fname);
      assert(::wxCopyFile(wxString(fname, wxConvUTF8), wxString(tmp, wxConvUTF8)));
    }
    sprintf(fname, "%s/%s/All.%s.constraints", DATA_DIR, g_classifierFile, g_classifer_extensions[behaviors->classifier_method]);
    if(FileExists(fname)) {
      sprintf(tmp, "%s/%s/All.%s.constraints", dir, g_classifierFile, g_classifer_extensions[behaviors->classifier_method]);
      fprintf(stderr, "Copying %s...\n", fname);
      assert(::wxCopyFile(wxString(fname, wxConvUTF8), wxString(tmp, wxConvUTF8)));
    }
  }
  sprintf(fname, "%s/rotate.txt", DATA_DIR);
  if(FileExists(fname)) {
    sprintf(tmp, "%s/rotate.txt", dir);
    fprintf(stderr, "Copying %s...\n", fname);
    assert(::wxCopyFile(wxString(fname, wxConvUTF8), wxString(tmp, wxConvUTF8)));
  }
  for(i = 0; i < behaviors->num; i++) {
    if(behaviors->behaviors[i].classifier_method) {
	sprintf(fname, "%s/%s/%s.%s", DATA_DIR, g_classifierFile, behaviors->behaviors[i].name, 
		g_classifer_extensions[behaviors->behaviors[i].values[j].classifier_method]);
	if(FileExists(fname)) {
	  sprintf(tmp, "%s/%s/%s.%s", dir, g_classifierFile, behaviors->behaviors[i].name,
		  g_classifer_extensions[behaviors->behaviors[i].values[j].classifier_method]);
	  fprintf(stderr, "Copying %s...\n", fname);
	  assert(::wxCopyFile(wxString(fname, wxConvUTF8), wxString(tmp, wxConvUTF8)));
	}
    }
    for(j = 0; j < behaviors->behaviors[i].num_values; j++) {
      if(behaviors->behaviors[i].values[j].classifier_method) {
	sprintf(fname, "%s/%s/%s_%s.%s", DATA_DIR, g_classifierFile, behaviors->behaviors[i].name, behaviors->behaviors[i].values[j].name, 
		g_classifer_extensions[behaviors->behaviors[i].values[j].classifier_method]);
	if(FileExists(fname)) {
	  sprintf(tmp, "%s/%s/%s_%s.%s", dir, g_classifierFile, behaviors->behaviors[i].name, behaviors->behaviors[i].values[j].name, 
		  g_classifer_extensions[behaviors->behaviors[i].values[j].classifier_method]);
	  fprintf(stderr, "Copying %s...\n", fname);
	  assert(::wxCopyFile(wxString(fname, wxConvUTF8), wxString(tmp, wxConvUTF8)));
	}
      }
    }
  }

  // Copy behaviors
  sprintf(fname, "%s/behaviors.txt", dir);
  fprintf(stderr, "Saving %s...\n", fname);
  assert(save_behaviors(fname, behaviors));

  // Copy training files
  for(i = 0; i < num_train_files; i++) {
    strcpy(tmp, train_list[i]);
    for(j = 0; j < (int)strlen(tmp); j++) {
      if(tmp[j] == '/' || tmp[j] == '\\')
	tmp[j] = '_';
    }
    char *ptr = tmp; 
    while(*ptr == '.') ptr++;
    sprintf(fname, "%s/%s", dir, ptr);
    tlist[i] = StringCopy(fname);
    fprintf(stderr, "Copying %s...\n", train_list[i]);
    assert(::wxCopyFile(wxString(train_list[i], wxConvUTF8), wxString(tlist[i], wxConvUTF8)));
  }
  sprintf(fname, "%s/train.txt", dir);
  fprintf(stderr, "Saving %s...\n", fname);
  assert(save_train_list(fname, tlist, num_train_files));
  for(i = 0; i < num_train_files; i++)
    free(tlist[i]);
  free(tlist);

  fprintf(stderr, "Finished backing up to %s\n", dir);
}

void LabelFrame::OnBackup(wxCommandEvent& WXUNUSED(event)) {
  if(IsBusy() || !PromptSave())
    return;
  wxDateTime now = wxDateTime::Now();
  char dir[1000], str[2000];
  sprintf(dir, "../backup/");
  wxString dstr = now.Format(_T("%c"), wxDateTime::CET);
  dstr.Replace(_T(" "), _T("_"));
  strcat(dir, dstr.mb_str()); 
  sprintf(str, "Backup all %d annotation files to %s?", num_train_files, dir);
  wxMessageDialog dlg(this, wxString(str, wxConvUTF8), wxT("Backup Files?"), wxYES_NO);
  int retval = dlg.ShowModal();
  wxCommandEvent evt;
  if(retval == wxID_YES) 
    Backup(dir);
}


void on_fit_frame(int t) {
  g_frame->SetFrame(t);
}

void lock_frame(int f) {
  g_frame->Lock(f);
}

void LabelFrame::OnPlayTimer(wxTimerEvent& event)
{
  if(play && frame+1 < blobs->num_frames)
    SetFrame(frame+1);
  else {
    timer->Stop();
    mutex->Lock();
    play = 0;
    SynchEnable();
    mutex->Unlock();
  }
}

void LabelFrame::SetFrame(int f) {
  mutex->Lock();
  img_window->UpdateBlob(preview_frame >= 0 ? &blobs->frames[preview_frame] : NULL, &blobs->frames[f]);
  preview_frame = frame = f;
  m_frame_number_label->SetLabel(wxString::Format(wxT("%i"), frame+1));
  
  for(int i = 0; i < behaviors->num; i++) 
    m_behavior_timelines[i]->SetFrame(frame);

  if(meta_behaviors)
    for(int i = 0; i < meta_behaviors->num; i++) 
      m_meta_behavior_timelines[i]->SetFrame(frame);

  m_graph->SetFrame(frame);
  mutex->Unlock();
}


void LabelFrame::PreviewFrame(int f) {
  mutex->Lock();
  img_window->UpdateBlob(&blobs->frames[preview_frame], &blobs->frames[f]);
  preview_frame = f;
  m_frame_number_label->SetLabel(wxString::Format(wxT("%i"), preview_frame+1));
  mutex->Unlock();
}

void LabelFrame::OnFitFrame(wxCommandEvent& WXUNUSED(event)) {
  mutex->Lock();
  Blob *prev = frame > 0 ? &blobs->frames[frame-1] : NULL;
  fit_model_one_frame(&blobs->frames[frame], prev, &params, 1, NULL, DEBUG_FIT_SKELETON);
  if(!blobs->frames[frame].features)
    blobs->frames[frame].features = (double*)malloc(2*sizeof(double)*((2*F_NUM_FEATURES)*(params.num_spine_lines+1)+F_NUM_GLOBAL_FEATURES));
  compute_features(&blobs->frames[frame],  prev, prev ? (1/(blobs->frames[frame].frame_time-prev->frame_time)) : 0, &params, 1);
  compute_global_features(blobs, frame, &params, params.expected_length, 1);
  compute_all_global_features(blobs, &params, 0);
  predict_behaviors(&blobs->frames[frame], behaviors, NULL);
  img_window->DrawBlobs();
  SetDirty(true);
  mutex->Unlock();
}

void LabelFrame::OnValidate(wxCommandEvent& event) {
  Lock(true);
  if(m_behavior_timelines[0]->Mode() == MODE_NONE) {
    int i;
    for(i = 0; i < behaviors->num; i++) {
      m_behavior_timelines[i]->SetMode(MODE_VALIDATE);
      m_behavior_timelines[i]->Timeline()->BeginSelection(frame);
      m_behavior_timelines[i]->Timeline()->RecomputeImage(false);
    }
    if(meta_behaviors) {
      for(i = 0; i < meta_behaviors->num; i++) {
	m_meta_behavior_timelines[i]->SetMode(MODE_VALIDATE);
	m_meta_behavior_timelines[i]->Timeline()->BeginSelection(frame);
	m_meta_behavior_timelines[i]->Timeline()->RecomputeImage(false);
      }
    }

    m_graph->SetMode(MODE_VALIDATE);
    m_graph->Graph()->BeginSelection(-1);
    m_graph->Graph()->RecomputeImage(false);
    m_flip_button->SetLabel(wxT("Cancel"));
  } else {
    for(int i = 0; i < behaviors->num; i++) {
      m_behavior_timelines[i]->SetMode(MODE_NONE);
      m_behavior_timelines[i]->Timeline()->BeginSelection(-1);
      m_behavior_timelines[i]->Timeline()->RecomputeImage(false);
    }
    if(meta_behaviors) {
      for(int i = 0; i < meta_behaviors->num; i++) {
	m_meta_behavior_timelines[i]->SetMode(MODE_NONE);
	m_meta_behavior_timelines[i]->Timeline()->BeginSelection(-1);
	m_meta_behavior_timelines[i]->Timeline()->RecomputeImage(false);
      }
    }

    m_graph->SetMode(MODE_NONE);
    m_graph->Graph()->BeginSelection(-1);
    m_graph->Graph()->RecomputeImage(false);
    m_validate_button->SetLabel(wxT("Validate^"));
  }
  Lock(false);
}

void LabelFrame::OnFinishValidate() {
  mutex->Lock();
  int s = frame, e = frame, i;
  m_behavior_timelines[0]->GetSelection(&s, &e);
  for(i = s; i <= e; i++) 
    blobs->frames[i].is_manual |= 3;
  SetDirty(true);
  img_window->DrawBlobs();

  for(i = 0; i < behaviors->num; i++) {
    m_behavior_timelines[i]->SetMode(MODE_NONE);
    m_behavior_timelines[i]->Timeline()->BeginSelection(-1);
    m_behavior_timelines[i]->RecomputeImage(false);
  }
  if(meta_behaviors) {
    for(i = 0; i < meta_behaviors->num; i++) {
      m_meta_behavior_timelines[i]->SetMode(MODE_NONE);
      m_meta_behavior_timelines[i]->Timeline()->BeginSelection(-1);
      m_meta_behavior_timelines[i]->RecomputeImage(false);
    }
  }
  m_graph->SetMode(MODE_NONE);
  m_graph->Graph()->BeginSelection(-1);
  m_graph->Graph()->RecomputeImage(false);

  m_validate_button->SetLabel(wxT("Validate^"));

  mutex->Unlock();
}

  
void LabelFrame::OnFitAll(wxCommandEvent& WXUNUSED(event)) {
  Fit(); 
  predict_all_behaviors(blobs, behaviors);
  if(meta_behaviors)
    ExtractMetaBehaviors(blobs, behaviors, meta_behaviors);
  SetDirty(true);
  SetFrame(frame);

  //return;
  //wxThread *_hThread = new FitThread(this);
  //_hThread->Create();
  //_hThread->Run();
}

void LabelFrame::OnFitBehaviors(wxCommandEvent& WXUNUSED(event)) {
  predict_all_behaviors(blobs, behaviors);
  if(meta_behaviors)
    ExtractMetaBehaviors(blobs, behaviors, meta_behaviors);
  SetDirty(true);
  SetFrame(frame);
}

void LabelFrame::OnFitSkeleton(wxCommandEvent& WXUNUSED(event)) {
  Fit();
  SetFrame(frame);
}


void LabelFrame::OnPredictBehaviors(wxCommandEvent& WXUNUSED(event)) {
  predict_all_behaviors(blobs, behaviors);
  if(meta_behaviors)
    ExtractMetaBehaviors(blobs, behaviors, meta_behaviors);
  mutex->Lock();
  for(int i = 0; i < behaviors->num; i++)
    m_behavior_timelines[i]->RecomputeImage(false);
  mutex->Unlock();
}

void LabelFrame::OnFlip(wxCommandEvent& WXUNUSED(event)) {
  Lock(true);
  if(m_behavior_timelines[0]->Mode() == MODE_NONE) {
    for(int i = 0; i < behaviors->num; i++) {
      m_behavior_timelines[i]->SetMode(MODE_FLIP);
      m_behavior_timelines[i]->Timeline()->BeginSelection(frame);
      m_behavior_timelines[i]->Timeline()->RecomputeImage(false);
    }
    if(meta_behaviors) {
      for(int i = 0; i < meta_behaviors->num; i++) {
	m_meta_behavior_timelines[i]->SetMode(MODE_FLIP);
	m_meta_behavior_timelines[i]->Timeline()->BeginSelection(frame);
	m_meta_behavior_timelines[i]->Timeline()->RecomputeImage(false);
      }
    }
    m_graph->SetMode(MODE_FLIP);
    m_graph->Graph()->BeginSelection(frame);
    m_graph->Graph()->RecomputeImage(false);
    m_flip_button->SetLabel(wxT("Cancel"));
  } else {
    for(int i = 0; i < behaviors->num; i++) {
      m_behavior_timelines[i]->SetMode(MODE_NONE);
      m_behavior_timelines[i]->Timeline()->BeginSelection(-1);
      m_behavior_timelines[i]->Timeline()->RecomputeImage(false);
    }
    if(meta_behaviors) {
      for(int i = 0; i < meta_behaviors->num; i++) {
	m_meta_behavior_timelines[i]->SetMode(MODE_NONE);
	m_meta_behavior_timelines[i]->Timeline()->BeginSelection(-1);
	m_meta_behavior_timelines[i]->Timeline()->RecomputeImage(false);
      }
    }
    m_graph->SetMode(MODE_NONE);
    m_graph->Graph()->BeginSelection(-1);
    m_graph->Graph()->RecomputeImage(false);
    m_flip_button->SetLabel(wxT("Flip^"));
  }
  Lock(false);
}

void LabelFrame::OnFinishFlip() {
  mutex->Lock();
  int s = frame, e = frame, i;
  m_behavior_timelines[0]->GetSelection(&s, &e);
  for(i = s; i <= e; i++) 
    flip_spine(&blobs->frames[i], &params);
  SetDirty(true);
  img_window->DrawBlobs();

  for(i = 0; i < behaviors->num; i++) {
    m_behavior_timelines[i]->SetMode(MODE_NONE);
    m_behavior_timelines[i]->Timeline()->BeginSelection(-1);
    m_behavior_timelines[i]->RecomputeImage(false);
  }
  if(meta_behaviors) {
    for(i = 0; i < meta_behaviors->num; i++) {
      m_meta_behavior_timelines[i]->SetMode(MODE_NONE);
      m_meta_behavior_timelines[i]->Timeline()->BeginSelection(-1);
      m_meta_behavior_timelines[i]->RecomputeImage(false);
    }
  }
  m_graph->SetMode(MODE_NONE);
  m_graph->Graph()->BeginSelection(-1);
  m_graph->Graph()->RecomputeImage(false);

  m_flip_button->SetLabel(wxT("Flip^"));

  mutex->Unlock();
}

void LabelFrame::OnZoomPlus(wxCommandEvent& WXUNUSED(event)) {
  SetTimelineZoom(timeline_zoom*1.1);
}
void LabelFrame::OnZoomMinus(wxCommandEvent& WXUNUSED(event)) {
  SetTimelineZoom(timeline_zoom/1.1);
}
void LabelFrame::SetTimelineZoom(double z) {
  timeline_zoom = z;
  for(int i = 0; i < behaviors->num; i++) {
    m_behavior_timelines[i]->SetZoom(timeline_zoom);
  }
  if(meta_behaviors) {
    for(int i = 0; i < meta_behaviors->num; i++) {
      m_meta_behavior_timelines[i]->SetZoom(timeline_zoom);
    }
  }
  m_graph->Graph()->SetZoom(timeline_zoom);
}
void LabelFrame::OnScrollLeft(wxCommandEvent& WXUNUSED(event)) {
  TimelineScroll(-1);
}
void LabelFrame::OnScrollRight(wxCommandEvent& WXUNUSED(event)) {
  TimelineScroll(1);
}
void LabelFrame::TimelineScroll(int delta) {
  for(int i = 0; i < behaviors->num; i++) {
    m_behavior_timelines[i]->SetScroll(delta);
  }
  if(meta_behaviors) {
    for(int i = 0; i < meta_behaviors->num; i++) {
      m_meta_behavior_timelines[i]->SetScroll(delta);
    }
  }
  m_graph->Graph()->SetScroll(delta);
}

void LabelFrame::Fit() {
  int old = frame;
  g_frame = this;
  mutex->Lock();
  play = 1;
  for(int i = 0; i < behaviors->num; i++)
    m_behavior_timelines[i]->Clear();
  if(meta_behaviors)
    for(int i = 0; i < meta_behaviors->num; i++)
      m_meta_behavior_timelines[i]->Clear();
  SynchEnable();
  mutex->Unlock();
  fit_model(blobs, &params, DEBUG_FIT_SKELETON, /*&on_fit_frame*/ NULL, &play);
  mutex->Lock();
  play = 0;
  SynchEnable();
  mutex->Unlock();
  SetFrame(old);
  SetDirty(true);
}

void LabelFrame::Clear() {
  mutex->Lock();
  if(blobs) {
    free_blob_sequence(blobs);
    img_window->ClearBlobs();
    if(behaviors && m_behavior_timelines){
      for(int i = 0; i < behaviors->num; i++){
	if(m_behavior_timelines[i]){
	  m_behavior_timelines[i]->Clear();
	}
      }
    }
    if(meta_behaviors && m_meta_behavior_timelines){
      for(int i = 0; i < meta_behaviors->num; i++){
	if(m_meta_behavior_timelines[i]){
	  m_meta_behavior_timelines[i]->Clear();
	}
      }
    }
    blobs = NULL;
  }
  mutex->Unlock();
  SetDirty(false);
  SynchEnable();
}



void LabelFrame::OnSize(wxSizeEvent &event) {
  wxSize sz = GetSize();
  timeline_splitter->SetSashPosition(sz.y-TIMELINE_PANEL_HEIGHT-TIMELINE_HEIGHT*((meta_behaviors ? 2 : 1)*behaviors->num+1));
  img_splitter->SetSashPosition(sz.x-CONTROL_WIDTH);
  wxPoint pos;
  for(int i = 0; i < behaviors->num; i++) {
    pos = m_behavior_timelines[i]->Timeline()->GetPosition();
    m_behavior_timelines[i]->Timeline()->SetMinSize(wxSize(sz.x-pos.x,TIMELINE_HEIGHT));
  }
  if(meta_behaviors) {
    for(int i = 0; i < meta_behaviors->num; i++) {
      pos = m_meta_behavior_timelines[i]->Timeline()->GetPosition();
      m_meta_behavior_timelines[i]->Timeline()->SetMinSize(wxSize(sz.x-pos.x,TIMELINE_HEIGHT));
    }
  }
  m_graph->Graph()->SetPosition(wxPoint(pos.x, -1));
  m_graph->Graph()->SetMinSize(wxSize(sz.x-pos.x,TIMELINE_HEIGHT));

  m_timeline_sizer->RecalcSizes();
  img_window->SetMinSize(wxSize(sz.x-CONTROL_WIDTH,sz.y-TIMELINE_PANEL_HEIGHT-TIMELINE_HEIGHT*((meta_behaviors ? 2 : 1)*behaviors->num+1)));
  mutex->Lock();
  if(blobs) {
    img_window->CropToSequence(blobs);
    img_window->DrawBlobs();
  }
  mutex->Unlock();
  event.Skip();
}


void LabelFrame::OnKeyDown(wxKeyEvent& event) {
  if(blobs && !IsBusy()) {
    if(event.GetKeyCode() == WXK_SPACE) {
      compute_features(&blobs->frames[frame], frame ? &blobs->frames[frame-1] : NULL, 
		       frame ? 1/(blobs->frames[frame].frame_time-blobs->frames[frame-1].frame_time) : 0, &params, 1);
      return;
    } else {
    }
  } 
  event.Skip();
}


void LabelFrame::SaveBehaviors() {
  char bname[400];
  sprintf(bname, "%s/behaviors.txt", DATA_DIR);
  assert(save_behaviors(bname, behaviors));
}


void BehaviorTimelineControls::SynchEnable(bool busy, BlobSequence *blobs) {
  m_edit_button->Enable(!m_is_meta);
  m_add_button->Enable(!m_is_meta && !busy && blobs);
}

BehaviorTimelineControls::BehaviorTimelineControls(LabelFrame *parent, wxWindow *panel, 
						   int g_ind, wxSizer *sizer, int ind, bool is_meta) {
  this->m_is_meta = is_meta;
  this->group_ind = g_ind;
  this->parent = parent;
  this->panel = panel;
  this->mode = MODE_NONE;

  BehaviorGroup *group = &parent->GetBehaviors()->behaviors[group_ind];
  m_menu = NULL;
  selected_behavior = -1;

  m_timeline = new TimelineWindow(panel, this, m_is_meta);

  m_label = new wxStaticText(panel, wxID_ANY, wxString(group->name, wxConvUTF8) + wxT(": "), wxDefaultPosition, wxSize(BEHAVIOR_LABEL_WIDTH,-1), wxALIGN_RIGHT);
  m_add_button = new wxButton(panel, parent->BehaviorID(group)+1000, wxT("Begin Add->"));
  m_edit_button = new wxButton(panel, parent->BehaviorID(group)+1000, wxT("Edit"));
  m_b_sizer = new wxFlexGridSizer(1,2,0,0);
  sizer->Add(m_label,         ind*3 + 0, wxRIGHT | wxALIGN_CENTER_VERTICAL);
  m_b_sizer->Add(m_edit_button,   0, wxLEFT | wxALIGN_CENTER_VERTICAL);
  m_b_sizer->Add(m_add_button,    1, wxLEFT | wxALIGN_CENTER_VERTICAL);
  sizer->Add(m_b_sizer, ind*3 + 1, wxLEFT | wxALIGN_CENTER_VERTICAL);
  
  m_timeline_sizer = new wxFlexGridSizer(1,3,0,0);
  //m_scroll_left = new wxButton(panel, Label_Scroll_Left, wxT("<"), wxDefaultPosition, wxSize(25,-1));
  //m_scroll_right = new wxButton(panel, Label_Scroll_Right, wxT(">"), wxDefaultPosition, wxSize(25,-1));
  //m_timeline_sizer->Add(m_scroll_left, 0, wxLEFT | wxALIGN_CENTER_VERTICAL);
  m_timeline_sizer->Add(m_timeline, 1, wxLEFT | wxALIGN_CENTER_VERTICAL);
  //m_timeline_sizer->Add(m_scroll_right, 2, wxLEFT | wxALIGN_CENTER_VERTICAL);

  sizer->Add(m_timeline_sizer,      ind*3 + 2, wxLEFT | wxALIGN_CENTER_VERTICAL);

  m_add_button->Connect(wxEVT_COMMAND_BUTTON_CLICKED, (wxObjectEventFunction)&LabelFrame::OnAddBehaviorButton, (wxObject*)this, parent);
  m_edit_button->Connect(wxEVT_COMMAND_BUTTON_CLICKED, (wxObjectEventFunction)&LabelFrame::OnEditBehaviorButton, (wxObject*)this, parent);
  //m_scroll_left->Connect(wxEVT_COMMAND_BUTTON_CLICKED, (wxObjectEventFunction)&LabelFrame::OnScrollLeft, (wxObject*)this, parent);
  //m_scroll_right->Connect(wxEVT_COMMAND_BUTTON_CLICKED, (wxObjectEventFunction)&LabelFrame::OnScrollRight, (wxObject*)this, parent);

  BuildMenus();
}

BehaviorTimelineControls::~BehaviorTimelineControls() {
  /*m_add_button->Disconnect(wxEVT_COMMAND_BUTTON_CLICKED, (wxObjectEventFunction)&LabelFrame::OnAddBehaviorButton, (wxObject*)this, parent);
  m_edit_button->Disconnect(wxEVT_COMMAND_BUTTON_CLICKED, (wxObjectEventFunction)&LabelFrame::OnEditBehaviorButton, (wxObject*)this, parent);
  if(group) {
    for(int i = 0; i < group->num_values; i++) {
      m_menu->Disconnect(wxEVT_COMMAND_MENU_SELECTED, (wxObjectEventFunction)&LabelFrame::OnBehaviorMenu, (wxObject*)this, parent);
    }
    }*/
  
  
  if(m_timeline) m_timeline->Cleanup();
  
}


GraphControls::GraphControls(LabelFrame *parent, wxWindow *panel, wxSizer *sizer, int ind) {
  this->parent = parent;
  this->panel = panel;

  m_graph = new GraphWindow(panel, this);
  
  names = new wxString[F_NUM_GOOD_GLOBAL_FEATURES+1];
  names[0] = wxT("None");
  for(int i = 0; i < F_NUM_GOOD_GLOBAL_FEATURES; i++)
    names[i+1] = wxString(g_global_feature_names[i], wxConvUTF8);

  
  m_label = new wxStaticText(panel, wxID_ANY, wxT("Graph: "), wxDefaultPosition, wxSize(BEHAVIOR_LABEL_WIDTH,-1), wxALIGN_RIGHT);
  m_feature_combo = new wxComboBox(panel, 1000, _T("None"), wxDefaultPosition, wxDefaultSize,
                                    F_NUM_GOOD_GLOBAL_FEATURES+1, names, wxCB_READONLY|wxCB_DROPDOWN);
  sizer->Add(m_label,         ind*3 + 0, wxRIGHT | wxALIGN_CENTER_VERTICAL);
  sizer->Add(m_feature_combo, ind*3 + 1, wxLEFT | wxALIGN_CENTER_VERTICAL);
  sizer->Add(m_graph,      ind*3 + 2, wxRIGHT | wxALIGN_CENTER_VERTICAL);

  m_feature_combo->Connect(wxEVT_COMMAND_COMBOBOX_SELECTED, (wxObjectEventFunction)&LabelFrame::OnFeatureCombo, (wxObject*)this, parent);
}

GraphControls::~GraphControls() {
  //m_feature_combo->Disconnect(wxEVT_COMMAND_COMBOBOX_SELECTED, (wxObjectEventFunction)&LabelFrame::OnFeatureCombo, (wxObject*)this, parent);
  //delete [] names;
  m_graph->Cleanup();
}



void LabelFrame::OnFeatureCombo() {
  m_graph->SetFeature(m_graph->GetSelection()-1);
  m_graph->Graph()->RecomputeImage(true);
}

void BehaviorTimelineControls::BuildMenus() {
  //if(m_menu) 
  //delete m_menu;

  m_menu = new wxMenu;
  BehaviorGroup *group = &parent->GetBehaviors()->behaviors[group_ind];
  for(int i = 0; i < group->num_values; i++) {
    m_menu->Append(1000*(parent->BehaviorID(group)+1)+i, wxString(group->values[i].abbreviation, wxConvUTF8) + wxT("(") + 
		   wxString(group->values[i].name, wxConvUTF8) + wxT(")"));
    m_menu->Connect(wxEVT_COMMAND_MENU_SELECTED, (wxObjectEventFunction)&LabelFrame::OnBehaviorMenu, (wxObject*)this, parent);

  }
}

void LabelFrame::OnBehaviorMenu(wxCommandEvent& evt) {
  m_behavior_timelines[evt.GetId()/1000-1]->OnMenu(evt);
}

void BehaviorTimelineControls::OnMenu(wxCommandEvent& evt) {
  parent->Lock(true);
  selected_behavior = evt.GetId()%1000;
  SetMode(MODE_ADD);
  m_timeline->BeginSelection(parent->GetFrame());
  m_timeline->RecomputeImage(false);
  m_add_button->SetLabel(wxT("Cancel"));
  parent->Lock(false);
}


void LabelFrame::OnAddBehaviorButton(wxCommandEvent& evt) {
  m_behavior_timelines[evt.GetId()%1000]->OnAddButton(evt);
}

void BehaviorTimelineControls::OnAddButton(wxCommandEvent& WXUNUSED(event)) {
  int s, e;

  parent->Lock(true);
  if(!m_timeline->GetSelection(&s, &e) || s < 0) {
    parent->Lock(false);
    panel->PopupMenu(m_menu, m_add_button->GetPosition().x+m_add_button->GetSize().x, m_add_button->GetPosition().y);
    return;
  } else if(m_timeline->GetSelection(&s, &e)) {
    selected_behavior = -1;
    m_add_button->SetLabel(wxT("Begin Add->"));
    SetMode(MODE_NONE);
    m_timeline->BeginSelection(-1);
    m_timeline->RecomputeImage(false);
  }
  parent->Lock(false);
}

void LabelFrame::OnEditBehaviorButton(wxCommandEvent& evt) {
  if(IsBusy() || !PromptSave())
    return;
  Clear();

  BehaviorManagerControls dlg(m_behavior_timelines[evt.GetId()%1000]);
  dlg.ShowModal();
}

void LabelFrame::RecomputeTimelines() {
  int i;
  for(i = 0; i < behaviors->num; i++) 
    m_behavior_timelines[i]->Timeline()->RecomputeImage(false);
  if(meta_behaviors) {
    for(i = 0; i < meta_behaviors->num; i++) 
      m_meta_behavior_timelines[i]->Timeline()->RecomputeImage(false);
  }
}


void BehaviorTimelineControls::OnFinishAdd() {
  int s, e, i;

  BlobSequence *blobs = parent->GetBlobs();
  if(selected_behavior >= 0 && m_timeline->GetSelection(&s, &e)) {
    for(i = s; i <= e; i++) {
      blobs->frames[i].is_manual |= 2;
      blobs->frames[i].behaviors[group_ind] = selected_behavior;
    }
    if(parent->GetMetaBehaviors())
      ExtractMetaBehaviors(blobs, parent->GetBehaviors(), parent->GetMetaBehaviors());
    
    parent->SetDirty(true);
    parent->GetBlobWindow()->DrawBlobs();
    parent->RecomputeTimelines();
    selected_behavior = -1;
    m_add_button->SetLabel(wxT("Begin Add->"));
  }
  SetMode(MODE_NONE);
  m_timeline->BeginSelection(-1);
  m_timeline->RecomputeImage(true);
}

void LabelFrame::SetSelection(int s, int e) {
  Lock(true);
  for(int i = 0; i < behaviors->num; i++) {
    m_behavior_timelines[i]->Timeline()->SetSelection(s,e);
    m_behavior_timelines[i]->Timeline()->RecomputeImage(false);
  }
  if(meta_behaviors) {
    for(int i = 0; i < meta_behaviors->num; i++) {
      m_meta_behavior_timelines[i]->Timeline()->SetSelection(s,e);
      m_meta_behavior_timelines[i]->Timeline()->RecomputeImage(false);
    }
  }
  m_graph->Graph()->SetSelection(s,e);
  m_graph->Graph()->RecomputeImage(false);
  
  Lock(false);
}

void LabelFrame::ProcessDirectory(const char *dirName) {
  MultiBlobSequence *all_blobs = import_multi_blob_sequence(dirName, params.num_spine_lines+1);
  save_multi_blob_sequence(all_blobs, behaviors, params.num_spine_lines+1, params.num_orientations, 0);

#ifdef USE_OPENMP
  #pragma omp parallel for 
#endif
  for (int i = 0; i < all_blobs->num_blobs; i++) {
    char anno[400];
    strcpy(anno, all_blobs->blobs[i]->fname);
    strip_extension(anno);
    strcat(anno, ".anno");
    if(FileExists(anno)) {
      fprintf(stderr, "Skipping %s because annotation file already exists\n", all_blobs->blobs[i]->fname);
      continue;
    }

    BlobSequence *blobs = all_blobs->blobs[i];
    params.world_to_pixel_scale = guess_world_to_pixel_scale(blobs->frames[0].contour, 
							     blobs->frames[0].num_pts);
    for(int j = 0; j < blobs->num_frames; j++) {
      if(blobs->frames[j].num_model_pts != params.num_spine_lines+1 && (blobs->frames[j].is_manual&1))
	blobs->frames[j].is_manual -= 1;
      if(!blobs->frames[j].features)
	blobs->frames[j].features = (double*)malloc(2*sizeof(double)*((2*F_NUM_FEATURES)*(params.num_spine_lines+1)+F_NUM_GLOBAL_FEATURES));
    }

    compute_deterministic_spine_attributes_blob_sequence(blobs, &params, 1);
    compute_all_global_features(blobs, &params, 0);
    fit_model(blobs, &params, DEBUG_FIT_SKELETON, NULL, NULL);
    predict_all_behaviors(blobs, behaviors);
    if(meta_behaviors)
      ExtractMetaBehaviors(blobs, behaviors, meta_behaviors);

    char fnameB[400];
    save_blob_sequence(blobs->fname, blobs, behaviors, params.num_spine_lines+1, params.num_orientations);
    sprintf(fnameB, "%s.bouts", blobs->fname);
    save_blob_behavior_bouts(fnameB, blobs, behaviors, 0);
    if(meta_behaviors) {
      sprintf(fnameB, "%s.bouts.%s", blobs->fname, g_meta_ext);
      save_blob_behavior_bouts(fnameB, blobs, meta_behaviors, 1);
    }
  }
}

void LabelFrame::DetectMetaBehaviors(const char *multiName) {
  MultiBlobSequence *all_blobs = load_multi_blob_sequence(multiName, behaviors);
  
  for (int i = 0; i < all_blobs->num_blobs; i++) {
    ExtractMetaBehaviors(all_blobs->blobs[i], behaviors, meta_behaviors);
    if(meta_behaviors) {
      char fnameB[400];
      sprintf(fnameB, "%s.bouts.%s", all_blobs->blobs[i]->fname, g_meta_ext);
      fprintf(stderr, "Extract meta behaviors %s...\n", fnameB);
      save_blob_behavior_bouts(fnameB, all_blobs->blobs[i], meta_behaviors, 1);
    }
  }
}

void LabelFrame::ProcessFile(const char *fname) {
  Open(fname);
  Fit(); 
  predict_all_behaviors(blobs, behaviors);
  if(meta_behaviors)
    ExtractMetaBehaviors(blobs, behaviors, meta_behaviors);
  SetDirty(true);
  SetFrame(frame);
  Save(blobs->fname);
}

void LabelFrame::OpenDirectory(wxCommandEvent& evt) {
  wxDirDialog fd(this, wxT("Open a Blob Directory"), wxString(DATA_DIR, wxConvUTF8));
  if(fd.ShowModal() == wxID_OK) {
    ProcessDirectory(fd.GetPath().mb_str());
  }
}


BehaviorManagerControls::BehaviorManagerControls(BehaviorTimelineControls *parent) : 
  wxDialog(parent->Parent(), wxID_ANY, wxString(wxT("Edit behavior group ")) + wxString(parent->Behavior()->name, wxConvUTF8))
{
  behavior_ind = parent->BehaviorInd();
  this->parent = parent;
  frame = parent->Parent();
  
  char bname[1000], classifier_dir[1000];
  sprintf(bname, "%s/behaviors.txt", DATA_DIR);
  sprintf(classifier_dir, "%s/%s", DATA_DIR, g_classifierFile);
  behaviors = load_behaviors(bname, classifier_dir);
  group = &behaviors->behaviors[behavior_ind];
  assert(behaviors && group);

  m_strings = NULL;
  GenerateBehaviorStrings();

  m_sizer = new wxFlexGridSizer(3,1,0,0);
  m_list = new wxListBox(this, wxID_ANY, wxDefaultPosition, 
			 wxSize(240, -1), group->num_values, m_strings, wxLB_SINGLE);
  m_sizer->Add(m_list, 0, wxCENTER | wxALIGN_TOP);

  m_b_sizer = new wxFlexGridSizer(1,3,0,0);
  m_b_sizer2 = new wxFlexGridSizer(1,2,0,0);
  m_add_button = new wxButton(this, wxID_ANY, wxT("Add"));
  m_edit_button = new wxButton(this, wxID_ANY, wxT("Edit"));
  m_delete_button = new wxButton(this, wxID_ANY, wxT("Delete"));
  m_save_button = new wxButton(this, wxID_ANY, wxT("Save"));
  m_cancel_button = new wxButton(this, wxID_ANY, wxT("Cancel"));
  m_b_sizer->Add(m_add_button, 0, wxLEFT | wxALIGN_TOP);
  m_b_sizer->Add(m_edit_button, 1, wxLEFT | wxALIGN_TOP);
  m_b_sizer->Add(m_delete_button, 2, wxLEFT | wxALIGN_TOP);
  m_sizer->Add(m_b_sizer, 1, wxCENTER | wxALIGN_TOP);
  m_b_sizer2->Add(m_save_button, 0, wxCENTER | wxALIGN_TOP);
  m_b_sizer2->Add(m_cancel_button, 1, wxCENTER | wxALIGN_TOP);
  m_sizer->Add(m_b_sizer2, 1, wxCENTER | wxALIGN_TOP);

  m_add_button->Connect(wxEVT_COMMAND_BUTTON_CLICKED, (wxObjectEventFunction)&BehaviorManagerControls::OnAddButton, NULL, this);
  m_delete_button->Connect(wxEVT_COMMAND_BUTTON_CLICKED, (wxObjectEventFunction)&BehaviorManagerControls::OnDeleteButton, NULL, this);
  m_edit_button->Connect(wxEVT_COMMAND_BUTTON_CLICKED, (wxObjectEventFunction)&BehaviorManagerControls::OnEditButton, NULL, this);
  m_save_button->Connect(wxEVT_COMMAND_BUTTON_CLICKED, (wxObjectEventFunction)&BehaviorManagerControls::OnSaveButton, NULL, this);
  m_cancel_button->Connect(wxEVT_COMMAND_BUTTON_CLICKED, (wxObjectEventFunction)&BehaviorManagerControls::OnCancelButton, NULL, this);

  SetSizerAndFit(m_sizer);
}

void BehaviorManagerControls::Cleanup() {
  //if(m_strings) delete [] m_strings;
  
  
}


void BehaviorManagerControls::OnAddButton(wxCommandEvent& event) {
  BehaviorDialog dlg(this, -1, wxT("Add ") + wxString(group->name, wxConvUTF8) + 
		     wxT(" Behavior"), NULL);
  if(dlg.ShowModal() == wxID_OK) {
    group->values = (BehaviorValue*)realloc(group->values, sizeof(BehaviorValue)*(group->num_values+1));
    group->values[group->num_values++] = dlg.GetBehavior();
    GenerateBehaviorStrings();
    m_list->InsertItems(1, &m_strings[group->num_values-1], group->num_values-1);
  }
}

void BehaviorManagerControls::GenerateBehaviorStrings() {
  if(m_strings) delete [] m_strings;

  m_strings = new wxString[group->num_values];
  for(int i = 0; i < group->num_values; i++)
    m_strings[i] = wxString(group->values[i].abbreviation, wxConvUTF8) +
      wxT(" (") + wxString(group->values[i].name, wxConvUTF8) + wxT(")");
}


void BehaviorManagerControls::OnEditButton(wxCommandEvent& event) {
  int sel = m_list->GetSelection();
  if(!frame->PromptSave() || sel <= 1)
    return;

  frame->Clear();

  BehaviorDialog dlg(this, -1, wxString(wxT("Edit Behavior ")) + 
		     wxString(group->values[sel].abbreviation, wxConvUTF8), &group->values[sel]);
  if(dlg.ShowModal() == wxID_OK) {
    group->values[sel] = dlg.GetBehavior();
  }
}
void BehaviorManagerControls::OnDeleteButton(wxCommandEvent& event) {
  int sel = m_list->GetSelection();
  if(!frame->PromptSave() || sel <= 1)
    return;

  frame->Clear();

  for(int i = sel; i < group->num_values-1; i++)
    group->values[i] = group->values[i+1];
  group->num_values--;
  GenerateBehaviorStrings();
  m_list->Delete(sel);
}
void BehaviorManagerControls::OnSaveButton(wxCommandEvent& event) { 
  wxString str = wxString(_T("Warning! Changing behaviors may invalidate existing behavior annotations or behavior detectors! Are you want to save all changes to the behavior list?"));
  wxMessageDialog dlg(this, str, wxT("Save Behaviors?"), wxYES_NO);
  int retval = dlg.ShowModal();
  wxCommandEvent evt;
  if(retval == wxID_YES) {
    char bname[1000];
    sprintf(bname, "%s/behaviors.txt", DATA_DIR);
    assert(save_behaviors(bname, behaviors));
    frame->LoadBehaviors();
    parent->BuildMenus();
    frame->SynchEnable();
    Close();
  }
}
void BehaviorManagerControls::OnCancelButton(wxCommandEvent& event) { 
  Close();
}

BehaviorDialog::BehaviorDialog (BehaviorManagerControls *parent, wxWindowID id, const wxString & title,
				BehaviorValue *b) : wxDialog( parent, id, title) {
  memset(&behavior, 0, sizeof(BehaviorValue));

  this->parent = parent;
  this->frame = parent->Parent()->Parent();


  edit_behavior = b;
  if(edit_behavior)
    behavior = *edit_behavior;

  char color[400];
  sprintf(color, "%x", behavior.color);
  m_behavior_sizer = new wxFlexGridSizer(4,2,0,0);
  m_behavior_name_label = new wxStaticText(this, wxID_ANY, wxT("Name: "));
  m_behavior_name = new wxTextCtrl(this, Behavior_Name, wxString(behavior.name, wxConvUTF8), wxDefaultPosition, wxSize(150,-1));
  m_behavior_abbrev_label = new wxStaticText(this, Behavior_Abbreviation, wxT("Abbreviation: "), wxDefaultPosition);
  m_behavior_abbrev = new wxTextCtrl(this, wxID_ANY, wxString(behavior.abbreviation, wxConvUTF8), wxDefaultPosition, wxSize(150,-1));
  m_behavior_color_label = new wxStaticText(this, Behavior_Abbreviation, wxT("Color: "), wxDefaultPosition);
  m_behavior_color = new wxTextCtrl(this, wxID_ANY, wxString(color, wxConvUTF8), wxDefaultPosition, wxSize(150,-1));
  
  m_behavior_sizer->Add(m_behavior_name_label);
  m_behavior_sizer->Add(m_behavior_name);
  m_behavior_sizer->Add(m_behavior_abbrev_label);
  m_behavior_sizer->Add(m_behavior_abbrev);
  m_behavior_sizer->Add(m_behavior_color_label);
  m_behavior_sizer->Add(m_behavior_color);

  m_ok_button = new wxButton(this, Behavior_Ok, wxT("Ok"));
  m_cancel_button = new wxButton(this, Behavior_Cancel, wxT("Cancel"));
  m_behavior_sizer->Add(m_ok_button);
  m_behavior_sizer->Add(m_cancel_button);
  
  SetSizerAndFit(m_behavior_sizer);
}

void BehaviorDialog::Cleanup() {
   
}


void BehaviorDialog::OnOk( wxCommandEvent & event ) { 
  char color[400];

  strcpy(behavior.name, m_behavior_name->GetValue().mb_str());
  strcpy(behavior.abbreviation, m_behavior_abbrev->GetValue().mb_str());
  strcpy(color, m_behavior_color->GetValue().mb_str());
  sscanf(color, "%x", &behavior.color);
  EndModal(wxID_OK);
}

void BehaviorDialog::OnCancel( wxCommandEvent & event ) {
  EndModal(wxID_CANCEL);
}


void BlobWindow::Crop(Blob *b) {
  int w, h;
  double min_x = INFINITY, max_x = -INFINITY, min_y = INFINITY, max_y = -INFINITY;
  blob_get_bounding_box(b ? &b : blobs, b ? 1 : num_blobs, &min_x, &min_y, &max_x, &max_y);
  offsets[0] = -min_x;
  offsets[1] = -min_y;

  double sx = GetSize().x, sy = GetSize().y;
  if(m_rotate == 90 || m_rotate == 270) { sx = GetSize().y, sy = GetSize().x; }
  double zx = sx / (params->world_to_pixel_scale*(max_x-min_x));
  double zy = sy / (params->world_to_pixel_scale*(max_y-min_y));
  zoom = my_min(zx, zy);
  w = sx;//ceil((max_x-min_x)*params->world_to_pixel_scale*zoom)+1;
  h = sy;//ceil((max_y-min_y)*params->world_to_pixel_scale*zoom)+1;
  
  img = cvCreateImage(cvSize(w, h), IPL_DEPTH_8U, 3);
  cvZero(img);
}
void BlobWindow::CropToSequence(BlobSequence *s) {
  int w, h;
  double min_x = INFINITY, max_x = -INFINITY, min_y = INFINITY, max_y = -INFINITY;
  for(int j = 0; j < s->num_frames; j++) {
    Blob *b = &s->frames[j];
    for(int i = 0; i < b->num_pts; i++) {
      if(b->contour[2*i] < min_x) min_x = b->contour[2*i];
      if(b->contour[2*i] > max_x) max_x = b->contour[2*i];
      if(b->contour[2*i+1] < min_y) min_y = b->contour[2*i+1];
      if(b->contour[2*i+1] > max_y) max_y = b->contour[2*i+1];
    }
  }
  offsets[0] = -min_x;
  offsets[1] = -min_y;

  double sx = GetSize().x, sy = GetSize().y;
  if(m_rotate == 90 || m_rotate == 270) { sx = GetSize().y, sy = GetSize().x; }
  double zx = sx / (params->world_to_pixel_scale*(max_x-min_x));
  double zy = sy / (params->world_to_pixel_scale*(max_y-min_y));
  zoom = my_min(zx, zy);

  w = sx;//ceil((max_x-min_x)*params->world_to_pixel_scale*zoom)+1;
  h = sy;//ceil((max_y-min_y)*params->world_to_pixel_scale*zoom)+1;
  
  img = cvCreateImage(cvSize(w, h), IPL_DEPTH_8U, 3);
  cvZero(img);
}

void BlobWindow::AddBlob(Blob *b) {
  blobs = (Blob**)realloc(blobs, sizeof(Blob*)*(num_blobs+1));
  blobs[num_blobs++] = b;
}

void BlobWindow::ClearBlobs() {
  if(blobs) {
    free(blobs);
    blobs = NULL;
    num_blobs = 0;
    cvReleaseImage(&img);
  }
}

void BlobWindow::UpdateBlob(Blob *old, Blob *b) {
  bool found = false;
  for(int i = 0; i < num_blobs; i++) {
    if(blobs[i] == old) {
      blobs[i] = b;
      found = true;
      break;
    }
  }
  if(!found)
    AddBlob(b);
  DrawBlobs();
}

void BlobWindow::DrawBlobs() {
  if(img)
    cvZero(img);
  else
    Crop(NULL);
  
  for(int i = 0; i < num_blobs; i++) {
    img = draw_blob(img, blobs[i], offsets, params, zoom, parent->GetBehaviors()); 
    char str[1000], tmp[400];
    strcpy(str, "");
    BehaviorGroups *behaviors = parent->GetBehaviors();
    for(int j = 0; j < behaviors->num; j++) {
      if(blobs[i]->behaviors[j] >= 0) {
	if(j) strcat(str, ", ");
	sprintf(tmp, "%s: %s", behaviors->behaviors[j].name, behaviors->behaviors[j].values[blobs[i]->behaviors[j]].name);
	strcat(str, tmp);
      }
    }
    strcpy(text, str);
  }
  UpdateImage();
}


void BlobWindow::OnSize(wxSizeEvent &event) {
  /*mutex->Lock();
  if(img) {
    DrawBlobs();
  }
  mutex->Unlock();*/
}


void BlobWindow::OnMouse(wxMouseEvent &event) {
  int xx,yy, tmp, i, j;

  CalcUnscrolledPosition(event.GetX(), event.GetY(), &xx, &yy);
  event.m_x = xx; event.m_y = yy;
  if(m_flip_x) 
    xx = GetSize().x-1-xx;
  if(m_flip_y) 
    yy = GetSize().y-1-yy;
  if(m_rotate == 270) { 
    tmp = xx;
    xx = GetSize().y-1-yy;
    yy = tmp;
  } else if(m_rotate == 90) { 
    tmp = xx;
    xx = yy;
    yy = GetSize().x-1-tmp;
  } else if(m_rotate == 180) { 
    xx = GetSize().x-1-xx;
    yy = GetSize().y-1-yy;
  } 


  if(parent->IsBusy() || !blobs)
    return;

  parent->Lock(true);
  //if(on_pt < 0 && event.ButtonDown()) {

  on_pt = -1;
  for(j = 0; j < num_blobs; j++) {
    for(i = 0; i <= params->num_spine_lines; i++) {
      if(blobs[j]->model_pts &&
	 sqrt(SQR(P2W(xx/zoom,0)-blobs[j]->model_pts[i].x) +
	      SQR(P2W(yy/zoom,1) - blobs[j]->model_pts[i].y)) <= 
	 5/zoom/params->world_to_pixel_scale) {
	on_pt = i;
	if(on_blob) on_blob->on_pt = 0;
	on_blob = blobs[j];
	break;
      }
    }
  }
  if(on_blob && on_blob->on_pt != on_pt+1) {
    on_blob->on_pt = on_pt+1;
    DrawBlobs();
  }
  if(on_pt >= 0 && event.Dragging() && on_blob->model_pts) {
    on_blob->model_pts[on_pt].x = P2W(xx/zoom,0);
    on_blob->model_pts[on_pt].y = P2W(yy/zoom,1);
    double *contour = on_blob->fixed_contour ? on_blob->fixed_contour : on_blob->contour;
    int num_pts = on_blob->fixed_contour ? on_blob->fixed_num_pts : on_blob->num_pts;
    
    
    if(on_pt < params->num_spine_lines) 
      compute_deterministic_spine_attributes(&on_blob->model_pts[on_pt], NULL, &on_blob->model_pts[on_pt+1], 
					     contour, num_pts, params, 1);
    else
      compute_deterministic_spine_attributes(&on_blob->model_pts[on_pt], &on_blob->model_pts[on_pt-1], NULL, 
					     contour, num_pts, params, 1);
    if(on_pt == params->num_spine_lines-1) 
      compute_deterministic_spine_attributes(&on_blob->model_pts[on_pt+1], &on_blob->model_pts[on_pt], NULL, 
					     contour, num_pts, params, 1);
    
          
    parent->SetDirty(true);
    on_blob->is_manual |= 1;
    DrawBlobs();
  } 
  parent->Lock(false);

  event.ResumePropagation(1); // Pass along mouse events (e.g. to parent)
  event.Skip();
}


TimelineWindow::TimelineWindow(wxWindow *w,  BehaviorTimelineControls *pa, bool is_meta, wxSize sz, long style) : 
  ImageWindow(w, pa->Parent(), sz, style) {
  cvInitFont(&font, CV_FONT_VECTOR0, 0.5f, 0.4f, 0, 2);
  start_frame = end_frame = -1;
  controls = pa;
  group_ind = pa->BehaviorInd();
  zoom = 1;
  scroll_x = 0;
  m_is_meta = is_meta;
}

void TimelineWindow::SetScroll(int d) {
  int w = GetSize().x;
  int w2 = w*zoom;

  scroll_x += 10;
  if(scroll_x + w2 > w) scroll_x = w-w2;
  if(scroll_x < 0) scroll_x = 0;
  EnableScrolling(true, false);
  Scroll(scroll_x,0);
}

void TimelineWindow::RecomputeImage(bool get_lock) { 
  if(get_lock) parent->Lock(true);
  BlobSequence *b = parent->GetBlobs();
  if(!b) {
    if(get_lock) parent->Lock(false);
    return;
  }

  BehaviorGroup *group = m_is_meta ? &parent->GetMetaBehaviors()->behaviors[group_ind] : &parent->GetBehaviors()->behaviors[group_ind];

  int w = GetSize().x*zoom;
  int h = GetSize().y;
  int hh, start;
  CvSize sz;
  int y[3], ymin;
  int last_pos = 2, last_last_pos = 1;

  if(img) cvReleaseImage(&img);
  img = cvCreateImage(cvSize(w, h), IPL_DEPTH_8U, 3);
  cvZero(img);
  double z = w / (double)b->num_frames;
  int i = 0, j, prev_prev_max_x = -100000, prev_max_x = -100000;
  while(i < b->num_frames) {
    // Group consecutive frames with the same behavior
    j = i;
    int beh = m_is_meta ? b->frames[i].meta_behaviors[group_ind] : b->frames[i].behaviors[group_ind];
    j = i;
    while(j < b->num_frames && beh == (m_is_meta ? b->frames[j].meta_behaviors[group_ind] : b->frames[j].behaviors[group_ind]))
      j++;

    if(i == 0) {
      cvGetTextSize("hello", &font, &sz, &ymin);
      hh = h;
      y[0] = my_min(hh-2, (hh+sz.height)/2);
      y[1] = h-6;
      y[2] = sz.height;
    }

    // Draw the group as a rectangle
    cvRectangle(img, cvPoint(z*i,0), cvPoint(z*j,hh), 
		CV_RGB((group->values[beh].color & 0xff0000)>>16,
		       (group->values[beh].color & 0x00ff00)>>8,
		       (group->values[beh].color & 0xff)), CV_FILLED);

    i = j;
  }

  // Draw lines outlining frames with verified skeletons
  i = 0;
  while(i < b->num_frames) {
    if(b->frames[i].is_manual) {
      start = i;
      while(i < b->num_frames && (b->frames[i].is_manual) && (b->frames[i].is_manual) == (b->frames[start].is_manual)) 
	i++;
      cvRectangle(img, cvPoint(z*start,0), cvPoint(z*i,3), CV_RGB(255, (b->frames[start].is_manual&1) ? 0 : 255, 
								  (b->frames[start].is_manual&2) ? 0 : 255), CV_FILLED);
    } else
      i++;
  }

  // Draw lines outlining frames with verified behaviors
  /*i = 0;
  while(i < b->num_frames) {
    if((b->frames[i].is_manual & 2)) {
      start = i;
      while(i < b->num_frames && (b->frames[i].is_manual & 2)) 
	i++;
      cvRectangle(img, cvPoint(z*start,0), cvPoint(z*i,3), CV_RGB(0,0,255), CV_FILLED);
    } else
      i++;
      }*/

  // Label each group
  i = 0;
  while(i < b->num_frames) {
    int beh = m_is_meta ? b->frames[i].meta_behaviors[group_ind] : b->frames[i].behaviors[group_ind];
    j = i;
    while(j < b->num_frames && beh == (m_is_meta ? b->frames[j].meta_behaviors[group_ind] : b->frames[j].behaviors[group_ind]))
      j++;

    cvGetTextSize(group->values[beh].abbreviation, &font, &sz, &ymin);
    int pos = 0;
    if((j+i-sz.width)/2 < prev_prev_max_x+4) 
      pos = 3 - last_last_pos - last_pos;
    else if((j+i-sz.width)/2 < prev_max_x+4) 
      pos = last_pos > 0 ? 0 : 1;
    else {
      pos = 0;
      last_pos = 2;
    }
    last_last_pos = last_pos;
    last_pos = pos;
    prev_prev_max_x = prev_max_x;
    prev_max_x = (j+i+sz.width)/2;
    cvPutText(img, group->values[beh].abbreviation, cvPoint((z*(i+j)-sz.width)/2, y[pos]), 
	      &font, CV_RGB(255,255,255));

    i = j;
  }

  DrawSelection();

  image_changed = true;
  Refresh(false);

  if(get_lock) parent->Lock(false);
}

void TimelineWindow::DrawSelection() {
  BlobSequence *b = parent->GetBlobs();
  assert(b && img);
  double z = img->width / (double)b->num_frames;
  int frame = parent->GetFrame(), i, j;
  cvRectangle(img, cvPoint(z*frame-1,0), cvPoint(z*frame+1,img->height), CV_RGB(0,255,0), CV_FILLED);

  if(end_frame >= 0) {
    int s = my_max(0,my_min(start_frame, end_frame)*z);
    int e = my_min(img->width,my_max(start_frame, end_frame)*z);
    unsigned char *ptr2 = (unsigned char*)img->imageData, *ptr;
    float alpha = .5;
    unsigned char c[3] = {0, 255, 0 };
    for(i = 0; i < img->height; i++, ptr2 += img->widthStep) {
      for(j = s, ptr = ptr2+s*3; j < e; j++, ptr += 3) {
	ptr[0] = (unsigned char)(ptr[0]*(1-alpha) + c[0]*alpha);
	ptr[1] = (unsigned char)(ptr[1]*(1-alpha) + c[1]*alpha);
	ptr[2] = (unsigned char)(ptr[2]*(1-alpha) + c[2]*alpha);
      }
    }
  }
}

void TimelineWindow::OnMouse(wxMouseEvent &event) {
  BlobSequence *b = parent->GetBlobs();

  if(parent->IsBusy() || !b)
    return;
  

  parent->Lock(true);
  int w = GetSize().x*zoom;
  double z = w / (double)b->num_frames;
  int f, xx, yy;
    
  CalcUnscrolledPosition(event.GetX(), event.GetY(), &xx, &yy);
  event.m_x = xx; event.m_y = yy;
  
  f = (int)(xx/z+.5);
  if(f < 0) f = 0;
  if(f >= b->num_frames) f =  b->num_frames-1;

  parent->Lock(false);

  if(event.Leaving()) 
    parent->PreviewFrame(parent->GetFrame());
  else if(event.ButtonDown()) {
    parent->SetFrame(f);
    if(end_frame >= 0) {
      if(controls->Mode() == MODE_ADD)
	controls->OnFinishAdd();
      else if(controls->Mode() == MODE_FLIP)
	parent->OnFinishFlip();
      else if(controls->Mode() == MODE_VALIDATE)
	parent->OnFinishValidate();
    }
  } else {
    if(end_frame >= 0) {
      end_frame = f;
      if(controls->Mode() == MODE_FLIP || controls->Mode() == MODE_VALIDATE) 
	parent->SetSelection(start_frame, end_frame);
    }
    parent->PreviewFrame(f);
    RecomputeImage(true);
  }
}

void TimelineWindow::OnSize(wxSizeEvent &event) {
  RecomputeImage(true);
}

void GraphWindow::RecomputeImage(bool get_lock) { 
  if(get_lock) parent->Lock(true);
  BlobSequence *b = parent->GetBlobs();
  if(!b) {
    if(get_lock) parent->Lock(false);
    return;
  }


  int i;
  int w = GetSize().x*zoom;
  int h = GetSize().y;
  double z = w / (double)b->num_frames;
  if(img) cvReleaseImage(&img);
  img = cvCreateImage(cvSize(w, h), IPL_DEPTH_8U, 3);
  cvZero(img);

  int y = h/2;
  if(feat_ind >= 0) {
    double min_f = INFINITY, max_f = -INFINITY, f;
    for(i = 0; i < b->num_frames; i++) {
      f = *(GLOBAL_FEATURES(&b->frames[i]) + feat_ind);
      if(f < min_f) min_f = f; 
      if(f > max_f) max_f = f; 
    }

    if(feat_ind > 0) {
      y = h-6-(h-6)*(0-min_f) / (max_f-min_f);
    }
  } 

  if(b->num_frames) {
    double st = b->frames[0].frame_time, et = b->frames[b->num_frames-1].frame_time;
    CvFont font;
    CvSize sz;
    int ymin;
    char str[400];

    char *ptr = b->fname;
    while(ptr && (ptr=strstr(ptr,"_")) != NULL) {
      float start, num, dur, inter;
      ptr++;
      if(sscanf(ptr, "%fs%fx%fs%fs", &start, &num, &dur, &inter) == 4) {
	for(int i = 0; i < num; i++) {
	  cvRectangle(img, cvPoint(w*(start+i*inter-st)/(et-st),0), cvPoint(w*(start+i*inter+dur-st)/(et-st),h), CV_RGB(0,255,255), CV_FILLED);
	}
	break;
      }
    }

    cvLine(img, cvPoint(0, y), cvPoint(w, y), CV_RGB(0,0,255), 1);
    cvInitFont(&font, CV_FONT_VECTOR0, 0.5f, 0.4f, 0, 2);
    
    double g, g2;
    if(et-st < 2) g = .1; 
    else if(et-st < 4) g = .2;
    else if(et-st < 10) g = .5;
    else if(et-st < 20) g = 1;
    else g = (double)((int)((et-st)/10));
    g2 = g / 5;
    double s = ((int)(st/g2))*g2;
    while(s < et) {
      int x = w*(s-st)/(et-st);
      cvLine(img, cvPoint(x, y-3), cvPoint(x, y+3), CV_RGB(0,0,255), 1);
      s += g2;
    }
    s = ((int)(st/g))*g;
    while(s < et) {
      int x = w*(s-st)/(et-st);
      if(g <= 1) sprintf(str, "%.1f", (float)s);
      else sprintf(str, "%d", (int)s);
      cvGetTextSize(str, &font, &sz, &ymin );
      cvLine(img, cvPoint(x, y-3), cvPoint(x, y+3), CV_RGB(0,0,255), 2);
      cvPutText(img, str, cvPoint(x-sz.width/2, y+3+sz.height), &font, CV_RGB(0,0,255));
      s += g;
    }
  }

  if(feat_ind >= 0) {
    double last_x = 0, last_y = 0, f, x, y;
    double min_f = INFINITY, max_f = -INFINITY;
    for(i = 0; i < b->num_frames; i++) {
      f = *(GLOBAL_FEATURES(&b->frames[i]) + feat_ind);
      if(f < min_f) min_f = f; 
      if(f > max_f) max_f = f; 
    }
    for(i = 0; i < b->num_frames; i++) {
      f = *(GLOBAL_FEATURES(&b->frames[i]) + feat_ind);
      x = z*i;
      y = h-6-(h-6)*(f-min_f) / (max_f-min_f);
      if(i) 
	cvLine(img, cvPoint(last_x, last_y), cvPoint(x, y), CV_RGB(255,0,255), 2);
      last_x = x; last_y = y;
    }
  }

  DrawSelection();

  image_changed = true;
  Refresh(false);

  if(get_lock) parent->Lock(false);
}

void GraphWindow::OnSize(wxSizeEvent &event) {
  RecomputeImage(true);
}

void GraphWindow::OnMouse(wxMouseEvent &event) {
  BlobSequence *b = parent->GetBlobs();

  if(parent->IsBusy() || !b)
    return;
  

  parent->Lock(true);
  int w = GetSize().x*zoom;
  double z = w / (double)b->num_frames;
  int f, xx, yy;
    
  CalcUnscrolledPosition(event.GetX(), event.GetY(), &xx, &yy);
  event.m_x = xx; event.m_y = yy;
  
  f = (int)(xx/z+.5);
  if(f < 0) f = 0;
  if(f >= b->num_frames) f =  b->num_frames-1;

  parent->Lock(false);

  if(event.Leaving()) 
    parent->PreviewFrame(parent->GetFrame());
  else if(event.ButtonDown()) {
    parent->SetFrame(f);
    if(end_frame >= 0) {
      if(controls->Mode() == MODE_FLIP)
	parent->OnFinishFlip();
      else if(controls->Mode() == MODE_VALIDATE)
	parent->OnFinishValidate();
    }
  } else {
    if(end_frame >= 0) {
      end_frame = f;
      if(controls->Mode() == MODE_FLIP || controls->Mode() == MODE_VALIDATE) 
	parent->SetSelection(start_frame, end_frame);
    }
    parent->PreviewFrame(f);
    RecomputeImage(true);
  }
}
