#ifndef __GUI_H
#define __GUI_H

#include "common.h"
#include "fit_model.h"


#include <cv.h>
#include <highgui.h>
#include <ml.h> 

// for compilers that support precompilation, includes "wx/wx.h".
#include "wx/wxprec.h"


#ifdef __BORLANDC__
    #pragma hdrstop
#endif

// for all others, include the necessary headers
#ifndef WX_PRECOMP
    #include "wx/app.h"
    #include "wx/log.h"
    #include "wx/frame.h"
    #include "wx/menu.h"

    #include "wx/button.h"
    #include "wx/checkbox.h"
    #include "wx/listbox.h"
    #include "wx/statbox.h"
    #include "wx/stattext.h"
    #include "wx/textctrl.h"
    #include "wx/msgdlg.h"
#endif

#include "wx/sysopt.h"
#include "wx/bookctrl.h"
#include "wx/sizer.h"
#include "wx/colordlg.h"
#include "wx/fontdlg.h"
#include "wx/textdlg.h"
#include "wx/panel.h"
#include "wx/slider.h"
#include "wx/image.h"
#include "wx/scrolwin.h"
#include "wx/dcclient.h"
#include "wx/filedlg.h"
#include "wx/splitter.h"
#include "wx/cmdline.h"
#include "wx/timer.h"
#include "wx/msgdlg.h"
#include "wx/combobox.h"
#include "wx/listbox.h"
#include "wx/dcbuffer.h"
#include "wx/textctrl.h"
#include "wx/dirdlg.h"


// control ids
enum
{
    Label_Quit = 100,
    Label_Add,
    Label_Save,
    Label_Open,
    Label_Delete,
    Label_Fit_Frame,
    Label_Fit_All,
    Label_Fit_Behaviors,
    Label_Fit_Skeleton,
    Label_Flip,
    Label_Prev_Frame,
    Label_Next_Frame,
    Label_Validate,
    Label_Train,
    Label_Backup,
    Label_Play,
    Label_Play_Timer,
    Label_Train_List,
    Label_Zoom_Minus,
    Label_Zoom_Plus,
    Label_Scroll_Left,
    Label_Scroll_Right,

    Behavior_Name,
    Behavior_Abbreviation,
    Behavior_Ok,
    Behavior_Cancel,
#if wxUSE_TOOLTIPS
    Label_SetTooltip,
#endif // wxUSE_TOOLTIPS
    Label_Enable
};

typedef enum {
  MODE_NONE=0,
  MODE_ADD,
  MODE_FLIP,
  MODE_VALIDATE,
} TimelineMode;


class BlobWindow;
class TimelineWindow;
class BehaviorTimelineControls;
class GraphControls;

// ----------------------------------------------------------------------------
// our classes
// ----------------------------------------------------------------------------

// Define a new frame type: this is going to be our main frame
class LabelFrame : public wxFrame
{
public:
    // ctor(s) and dtor
  LabelFrame(const wxString& title, const char *arg);
    virtual ~LabelFrame();

    void SetFrame(int f);
    void PreviewFrame(int f);
    void Play();
    void Fit();
    void SetDirty(bool d= true) { dirty = d; SynchEnable(); SynchFileName(); }
    bool PromptSave();
    BlobSequence *GetBlobs() { return blobs; }
    int GetFrame() { return frame; }
    BehaviorGroups *GetBehaviors() { return behaviors; }
    BehaviorGroups *GetMetaBehaviors() { return meta_behaviors; }
    void Lock(bool l) { 
      if(l) mutex->Lock(); 
      else mutex->Unlock(); 
    }
    void SynchEnable();
    void SynchFileName();
    bool IsBusy() { return play; }
    void GenerateTrainListStrings();
    void GenerateBehaviorStrings();
    void Clear();
    void SaveBehaviors();
    void OnKeyDown(wxKeyEvent& event); 
    FitParams *Params() { return &params; }
    BlobWindow *GetBlobWindow() { return img_window; }
    void OnBehaviorMenu(wxCommandEvent& event);
    void OnAddBehaviorButton(wxCommandEvent& event);
    void OnEditBehaviorButton(wxCommandEvent& event);
    void OnScrollLeft(wxCommandEvent& WXUNUSED(event));
    void OnScrollRight(wxCommandEvent& WXUNUSED(event));
    void OpenDirectory(wxCommandEvent& evt);
    void OnFinishValidate();
    void OnFinishFlip();
    void SetSelection(int s, int e);
    void OnFeatureCombo();
    void SetTimelineZoom(double z);
    void TimelineScroll(int s);
    void Cleanup();
    void Backup(const char *);
//  void LoadBehaviors(const char *classifier_appendix = ""); // CSC: parameter with default value added    
    void LoadBehaviors();
    void ProcessDirectory(const char *dirName);
    void ProcessFile(const char *fname);
    void crossValidate(); // CSC
    void SetBlobs(BlobSequence *b, int is_import);
    void RecomputeTimelines();
    void DetectMetaBehaviors(const char *multiName);

    int BehaviorID(BehaviorGroup *g) {
      if(g < behaviors->behaviors || g >= behaviors->behaviors+behaviors->num)
	return -1;
      return g-behaviors->behaviors;
    }

protected:
    // event handlers
    void OnClose(wxCloseEvent& event);
    void OnExit(wxCommandEvent& event);
    void OnOpen(wxCommandEvent& event);
    void OnSave(wxCommandEvent& event);
    void OnSaveAs(wxCommandEvent& event);
    void OnAdd(wxCommandEvent& event);
    void OnDelete(wxCommandEvent& event);
    void OnBackup(wxCommandEvent& event);
    void Open(const char *);  
    void Save(const char *);  

    void OnFitFrame(wxCommandEvent& event);
    void OnFitAll(wxCommandEvent& event);
    void OnFitBehaviors(wxCommandEvent& event);
    void OnFitSkeleton(wxCommandEvent& event);
    void OnFlip(wxCommandEvent& event);
    void OnPrevFrame(wxCommandEvent& event);
    void OnNextFrame(wxCommandEvent& event);
    void OnPlay(wxCommandEvent& event);
    void OnPlayTimer(wxTimerEvent& event);
    void OnTrainList(wxCommandEvent& event);
    void OnValidate(wxCommandEvent& event);
    void OnTrain(wxCommandEvent& event);
    void OnPredictBehaviors(wxCommandEvent& WXUNUSED(event));
    void OnZoomPlus(wxCommandEvent& event);
    void OnZoomMinus(wxCommandEvent& event);

    void OnSize(wxSizeEvent &event);


private:
  bool dirty;
  wxTimer *timer;
  wxMutex *mutex;
  double timeline_zoom;

  int sel_timeline;

  int start_frame;
  int play;

  int frame, preview_frame;
  BlobSequence *blobs;
  FitParams params;

  BehaviorGroups *behaviors, *meta_behaviors;

  char **train_list;
  wxString *train_list_strings;
  int num_train_files;
  int sel_train_file;

  // the panel containing everything
  wxPanel *m_control_panel;
  wxPanel *m_timeline_panel;
  BlobWindow *img_window;
  wxSplitterWindow *timeline_splitter;
  wxSplitterWindow *img_splitter;
  wxSizer *m_timeline_sizer;
  wxSizer *m_timeline_sizer2;
  wxSizer *m_control_sizer;

  BehaviorTimelineControls *m_behavior_timelines[MAX_BEHAVIOR_GROUPS];
  BehaviorTimelineControls *m_meta_behavior_timelines[MAX_BEHAVIOR_GROUPS];
  GraphControls *m_graph;

  wxSizer *m_frame_sizer;
  wxStaticText *m_frame_label, *m_frame_number_label;


  wxButton *m_next_frame_button;
  wxButton *m_play_button;
  wxButton *m_prev_frame_button;

  wxButton *m_fit_frame_button;
  wxButton *m_fit_all_button;
  wxButton *m_fit_behaviors_button;
  wxButton *m_fit_skeleton_button;
  wxButton *m_flip_button;
  wxButton *m_validate_button;
  wxButton *m_train_button;
  wxButton *m_backup_button;
  wxButton *m_zoom_plus_button;
  wxButton *m_zoom_minus_button;

  wxStaticText *m_train_label;
  wxListBox *m_train_list;
  wxSizer *m_file_sizer;
  wxStaticText *m_file_label;
  wxButton *m_save_button;
  wxButton *m_cancel_button;
  wxButton *m_add_button;
  wxButton *m_dir_button;
  wxButton *m_delete_button;

  // any class wishing to process wxWidgets events must use this macro
  DECLARE_EVENT_TABLE()
};



// Define a new application type, each program should derive a class from wxApp
class LabelApp : public wxApp
{
  wxString fname;
  LabelFrame *m_frame;

public:

    // override base class virtuals
    // ----------------------------

    // this one is called on application startup and is a good place for the app
    // initialization (doing it here and not in the ctor allows to have an error
    // return: if OnInit() returns false, the application terminates)
    virtual bool OnInit();    
    virtual void OnInitCmdLine(wxCmdLineParser& parser);
    virtual bool OnCmdLineParsed(wxCmdLineParser& parser);
    void OnKeyDown(wxKeyEvent& event) { m_frame->OnKeyDown(event); }
 
private:
DECLARE_EVENT_TABLE() 
};


inline wxImage wx_from_cv(IplImage *img, int rotate, int flip_x, int flip_y) {
  wxImage wx(rotate == 90 || rotate == 270 ? img->height : img->width, 
	     rotate == 90 || rotate == 270 ? img->width : img->height, (unsigned char*) malloc(img->imageSize), false);
  int dx=3, dy=img->width*3, sx=0, sy=0;

  if(rotate == 90 || rotate == 270) {
    if(rotate == 90) { dx = img->height*3; sx = 0; dy = -3; sy = (img->height-1)*3; }
    if(rotate == 270) { dx = -img->height*3; sx = (img->width-1)*img->height*3; dy = 3; sy = 0; }
    if(flip_x) { sy = dy > 0 ? (img->height-1)*3 : 0; dy = -dy; }
    if(flip_y) { sx = dx > 0 ? (img->width-1)*img->height*3 : 0; dx = -dx; }
  } else {
    if(rotate == 180) {
      dx = -3; sx = (img->width-1)*3; dy = -img->width*3; sy = (img->height-1)*img->width*3;
    }
    if(flip_x) { sx = dx > 0 ? (img->width-1)*3 : 0; dx = -dx; }
    if(flip_y) { sy = dy > 0 ? (img->height-1)*img->width*3 : 0; dy = -dy; }
  }


  int i, j;
  unsigned char *dst = (unsigned char*)wx.GetData(), *dst2;
  unsigned char *ptr, *ptr2 = (unsigned char*)img->imageData;
  for(i = 0, dst2 = dst+sy+sx; i < img->height; i++, ptr2 += img->widthStep, dst2 += dy) {
    for(j = 0, ptr = ptr2, dst = dst2; j < img->width; j++, dst += dx, ptr += 3) {
      dst[0] = ptr[2];
      dst[1] = ptr[1];
      dst[2] = ptr[0];
    }
  }
  return wx;
}

// A scrolled window for showing an image.
class ImageWindow: public wxScrolledWindow {
protected:
  wxImage wx_img;
  wxBitmap bmp;
  IplImage *img;
  bool image_changed;
  LabelFrame *parent;
  bool stretch;
  int m_rotate, m_flip_x, m_flip_y;
  char text[1000];

public:
 ImageWindow(wxWindow *w, LabelFrame *pa, wxSize sz=wxDefaultSize, long style=0)
   : wxScrolledWindow(w, wxID_ANY, wxDefaultPosition, sz, style), bmp(0,0)
    {
      parent = pa;
      img = NULL;
      image_changed = false;
      stretch = false;
      m_rotate = m_flip_x = m_flip_y = 0;
      strcpy(text, "");
    }
  void Cleanup() {
    if(img) cvReleaseImage(&img);
  }

  void SetRotate(int r, int x, int y) { m_rotate = r; m_flip_x = x; m_flip_y = y; }

    void Create(wxWindow *parent, wxWindowID id = -1)
    {
        wxScrolledWindow::Create(parent, id);
    }

    void LoadImage(wxImage &image) {
        bmp = wxBitmap(image);
	if(!stretch) {
	  SetVirtualSize(bmp.GetWidth(), bmp.GetHeight());
	  wxClientDC dc(this);
	  PrepareDC(dc);
	  dc.DrawBitmap(bmp, 0, 0);
	} 
    }

    void UpdateImage() {
      if(img) {
	image_changed = true;
	Refresh(false);
      }
    }

protected:

    void OnPaint(wxPaintEvent &event) {
      if(image_changed) {
	parent->Lock(true);
	wx_img = wx_from_cv(img, m_rotate, m_flip_x, m_flip_y);
	LoadImage(wx_img);
	image_changed = false;
	parent->Lock(false);
      }
      wxPaintDC dc(this);
      PrepareDC(dc);
      dc.DrawBitmap(bmp, 0,0, true);
      if(strlen(text)) {
	wxCoord w, h;
	dc.SetTextForeground( wxColour(0,255,0));
	dc.GetTextExtent(wxString(text, wxConvUTF8), &w, &h);
	dc.DrawText(wxString(text, wxConvUTF8), 0, GetSize().y-h);
      }
    }
	
    void OnEraseBackground(wxEraseEvent &event) {
    }
private:
    DECLARE_EVENT_TABLE()
};


class BlobWindow: public ImageWindow {
protected:
  Blob **blobs;
  int num_blobs;

  FitParams *params;
  double offsets[2];

  int on_pt;
  Blob *on_blob;

  double zoom;

public:
 BlobWindow(wxWindow *w, LabelFrame *pa, wxSize sz=wxDefaultSize, long style=0)
   : ImageWindow(w, pa, sz, style)
    {
      params = pa->Params();
      blobs = NULL;
      img = NULL;
      on_blob = NULL;
      on_pt = -1;
      zoom = 5;
      num_blobs = 0;
    }
  void Cleanup() { 
    ClearBlobs(); 
  }

  void Crop(Blob *b);
  void CropToSequence(BlobSequence *s);
  void AddBlob(Blob *b);
  void ClearBlobs();
  void UpdateBlob(Blob *old, Blob *b);
  void DrawBlobs();

protected:
    void OnMouse(wxMouseEvent &event);
    void OnSize(wxSizeEvent &event);

private:
    DECLARE_EVENT_TABLE()
};

class TimelineWindow : public ImageWindow {
  CvFont font;
  BehaviorTimelineControls *controls;
  int group_ind;
  bool m_is_meta;
 protected:
  int start_frame, end_frame;
  double zoom;
  int scroll_x;


public:
  TimelineWindow(wxWindow *w,  BehaviorTimelineControls *pa, bool is_meta=false, wxSize sz=wxDefaultSize, long style=wxBORDER_SUNKEN);

  void Clear() { start_frame = end_frame = -1; }
  virtual void RecomputeImage(bool get_lock=true);
  void SetScroll(int x);

  void BeginSelection(int f) { start_frame = end_frame = f; }
  void SetSelection(int s, int e) { start_frame = s; end_frame = e; }
  void DrawSelection();

  bool GetSelection(int *s, int *e) { 
    if(start_frame >= 0 && end_frame >= 0) {
      if(start_frame < end_frame) { 
	*s = start_frame;
	*e = end_frame;
      } else {
	*e = start_frame;
	*s = end_frame;
      }
      return true;
    } 
    return false;
  }

  void FrameChanged(int f) { 
    if(end_frame >= 0) 
      end_frame = f; 
  }
  void SetZoom(double z) { 
    zoom = z;
    SetVirtualSize(GetSize().x*zoom,GetSize().y); 
    RecomputeImage(true); 
  }

protected:
    void OnMouse(wxMouseEvent &event);
    void OnSize(wxSizeEvent &event);

private:
    DECLARE_EVENT_TABLE()
};

class GraphWindow : public TimelineWindow {
  CvFont font;
  GraphControls *controls;
  BehaviorGroup *group;
  int feat_ind;

public:
 GraphWindow(wxWindow *w,  GraphControls *pa, wxSize sz=wxDefaultSize, long style=wxBORDER_SUNKEN) :
  TimelineWindow(w, (BehaviorTimelineControls*)pa, false, sz,  style) { 
    feat_ind = -1; 
    controls = pa; 
  }
  
  void SetFeature(int f) { feat_ind = f; }
  void RecomputeImage(bool get_lock=true);

protected:
    void OnMouse(wxMouseEvent &event);
    void OnSize(wxSizeEvent &event);

private:
    DECLARE_EVENT_TABLE()
};

class BehaviorTimelineControls {
  int group_ind;
  bool m_is_meta;

  wxMenu *m_menu;
  wxButton *m_add_button;
  wxButton *m_edit_button;
  wxButton *m_scroll_left, *m_scroll_right;
  TimelineWindow *m_timeline;
  wxFlexGridSizer *m_b_sizer, *m_timeline_sizer;

  int selected_behavior;
  TimelineMode mode;

 protected:
  wxWindow *panel;
  LabelFrame *parent;
  wxStaticText *m_label;

 public:
  ~BehaviorTimelineControls();
  BehaviorTimelineControls(LabelFrame *parent, wxWindow *panel, int group_ind, wxSizer *sizer, int ind, bool is_meta=false);
  BehaviorTimelineControls() { m_is_meta = false; group_ind = -1; m_menu = NULL; m_add_button = m_edit_button = NULL; m_timeline = NULL; selected_behavior = -1; mode = MODE_NONE; }

  void SetFrame(int f) {
    m_timeline->FrameChanged(f);
    m_timeline->RecomputeImage(false);
  }
  TimelineMode Mode() { return mode; }
  void SetMode(TimelineMode m) { mode = m; }
  void RecomputeImage(bool f) { m_timeline->RecomputeImage(f); }
  void SynchEnable(bool busy, BlobSequence *blobs);
  bool GetSelection(int *s, int *e) { 
    return m_timeline->GetSelection(s, e);
  }
  void Clear() { m_timeline->Clear(); }
  void BuildMenus();
  TimelineWindow *Timeline() { return m_timeline; }
  BehaviorGroup *Behavior() { return &parent->GetBehaviors()->behaviors[group_ind]; }
  LabelFrame *Parent() { return parent; }
  void OnFinishAdd();
  int BehaviorInd() { return group_ind; }
  void Redraw() { m_timeline->Refresh(false); }
  void SetZoom(double z) { m_timeline->SetZoom(z); }
  void SetScroll(int x) { m_timeline->SetScroll(x); }

public:
    void OnMenu(wxCommandEvent& event);
    void OnAddButton(wxCommandEvent& event);
};

class GraphControls : public BehaviorTimelineControls{
  GraphWindow *m_graph;
  wxComboBox *m_feature_combo;
  wxString *names;

 public:
  ~GraphControls();
  GraphControls(LabelFrame *parent, wxWindow *panel, wxSizer *sizer, int ind);
  void SetFeature(int f) { m_graph->SetFeature(f); }
  int GetSelection() { return m_feature_combo->GetSelection(); }
  GraphWindow *Graph() { return m_graph; }
  void SetFrame(int f) {
    m_graph->FrameChanged(f);
    m_graph->RecomputeImage(false);
  }
};

class BehaviorManagerControls : public wxDialog {
  BehaviorGroup *group;
  BehaviorTimelineControls *parent;
  LabelFrame *frame;
  int behavior_ind;
  BehaviorGroups *behaviors;

  wxSizer *m_sizer;
  wxSizer *m_b_sizer, *m_b_sizer2;
  wxStaticText *m_label;
  wxListBox *m_list;
  wxButton *m_edit_button;
  wxButton *m_delete_button;
  wxButton *m_add_button;
  wxButton *m_save_button;
  wxButton *m_cancel_button;

  wxString *m_strings;

 public:
  BehaviorManagerControls(BehaviorTimelineControls *par);
  void Cleanup();
  void GenerateBehaviorStrings();
  BehaviorTimelineControls *Parent() { return parent; }
  

 protected:
  void OnAddButton(wxCommandEvent& event);
  void OnEditButton(wxCommandEvent& event);
  void OnDeleteButton(wxCommandEvent& event);
  void OnSaveButton(wxCommandEvent& event);
  void OnCancelButton(wxCommandEvent& event);
};

class BehaviorDialog : public wxDialog
{
  BehaviorValue *edit_behavior;
  BehaviorValue behavior;
  LabelFrame *frame;
  BehaviorManagerControls *parent;

  wxSizer *m_behavior_sizer;
  wxStaticText *m_behavior_name_label;
  wxTextCtrl *m_behavior_name;
  wxStaticText *m_behavior_abbrev_label;
  wxTextCtrl *m_behavior_abbrev;
  wxStaticText *m_behavior_color_label;
  wxTextCtrl *m_behavior_color;

  wxButton *m_ok_button;
  wxButton *m_cancel_button;

 public: 
  BehaviorDialog(BehaviorManagerControls *parent, wxWindowID id, const wxString & title,
		 BehaviorValue *b);
  void Cleanup();
  BehaviorValue GetBehavior() { return behavior; }

 private:
  void OnOk( wxCommandEvent & event );
  void OnCancel( wxCommandEvent & event );

  DECLARE_EVENT_TABLE();
};

class PlayThread : public wxThread {
  LabelFrame *f;
public:
  PlayThread(LabelFrame *f) : wxThread(wxTHREAD_JOINABLE){ this->f = f; }

  void *Entry(void){ f->Play(); return 0; }
};

class FitThread : public wxThread {
  LabelFrame *f;
public:
  FitThread(LabelFrame *f) : wxThread(wxTHREAD_JOINABLE){ this->f = f; }

  void *Entry(void){ f->Fit(); return 0; }
};

#endif

