
#if DEBUG > 0 
extern char *g_currFile; // CSC 20110420: hack to pass current filename for debug purposes
#endif

#ifndef  __SVM_STRUCT_API_BEHAVIOR_SEQUENCE__
#define __SVM_STRUCT_API_BEHAVIOR_SEQUENCE__


#include "structured_svm.h"


#define USE_DURATION_COST 1
#define USE_NEW_LOSS 0
#define SCALED_LOSS 1
#define HAMMING_LOSS 0

#define ALLOW_SAME_TRANSITIONS
#define DEBUG 1
#define MAX_FILENAME 1000
#define MAX_FEATURES 1000

#if DEBUG > 0
//#define MAX_FEATURES 1000
extern const char *bout_feature_names[]; // initialized in svm_struct_api_behavior_sequence.cpp; used to display feature names (in combination with g_feature_map)
//extern const char *g_feature_names[MAX_FEATURES];
#endif

#define NONE_BEHAVIOR 0



#define INTERPOLATE_INTEGRAL_FEATURES
#define MAX_REGIONS 10

// Assume no score (other than transition costs) for the "none" behavior.  Effectively, removes all
// bout features, unary score, and duration score for the none behavior
// This allows an exhaustive search for all bout durations for the "none" behavior
#define NONE_CLASS_HAS_NO_SCORE 1


// Single behavior class. 
typedef struct _Behavior {
  char name[400];
  char abbreviation[400];
  unsigned int color;
} Behavior;

// Group of mutually-exclusive behaviors that should be detected dependently. 
typedef struct _Behaviors {
  Behavior *values;
  int num_values;
} Behaviors;

typedef struct {
  char *name;
  double *histogram_thresholds;  
  int num_histogram_bins;
} FrameFeature;

typedef enum {
  B_SUM,   // sum of per-frame features in some temporal region
  B_AVE,   // average of per-frame features in some temporal region
  B_VAR,   // variance of per-frame features in some temporal region
  B_DEV,   // standard deviation of per-frame features in some temporal region
  B_RAW,   // raw difference in per-frame features at two different locations
  B_MIN,   // min of per-frame features in some temporal region
  B_MAX,   // max of per-frame features in some temporal region
  B_DUR,   // the duration in frames
} BoutOperator;

typedef struct {
  int frame_feature;  // index of the frame feature to operate over
  int hist_bin;       // if >= 0, use 1[bin(frame_feature)==hist_bin] instead of frame_feature
  BoutOperator op;    // take the sum, mean, min, max, etc.
  double t_start, t_end; // defines the start and end of the region, in the reference frame where 
                      // t=0 is the bout start and t=1 is the bout end
  bool is_global;     // override t_start and t_end and use the entire video sequence
  int frame_coords;   // normally this is 0.  If frame_coords=1 t_start and t_end are relative offsets (in frames) from the bout start. If    frame_coords=2 they are relative offsets from the bout end
} BoutRegionFeature;

typedef struct {
  BoutRegionFeature regions[MAX_REGIONS];  // the response of this bout feature is the weighted sum of region features
  int num_regions;
  double weights[MAX_REGIONS];
  bool time_normalize;                     // If true, divide by the bout duration
  bool absolute_value;                     // If true, take the absolute value of the weighted sum f
  double thresh;                           // If num_thresholds=1, then make this feature into a decision stump 1[f < thresh]
  int num_thresholds;                      // If num_thresholds>1, assumes bout_features[i,...,i+num_thresholds-1] are clones of this one and collectively define a set of histogram bins
  double mu, gamma;                        // for normalization
  char *name;
} BoutFeature;


#define REGION_FEATURE_FAST() (integral_features[f->frame_feature][it2] - \
			       integral_features[f->frame_feature][it1])
#define REGION_SQR_FEATURE_FAST() (integral_sqr_features[f->frame_feature][it2] - \
				   integral_sqr_features[f->frame_feature][it1])
#define FRAME_FEATURE_FAST(it) (safe_features[f->frame_feature][it])
#define REGION_HIST_FEATURE_FAST() (integral_histogram_features[f->frame_feature][f->hist_bin][it2] - \
				    integral_histogram_features[f->frame_feature][f->hist_bin][it1])
#define FRAME_HIST_FEATURE_FAST(it) (histogram_bins[f->frame_feature][it]==f)

#ifdef INTERPOLATE_INTEGRAL_FEATURES
#define REGION_FEATURE() (integral_features[f->frame_feature][it2] +	\
			  r2*safe_features[f->frame_feature][it2] -	\
			  integral_features[f->frame_feature][t] +	\
			  (1-r1)*safe_features[f->frame_feature][it1])
#define REGION_SQR_FEATURE() (integral_sqr_features[f->frame_feature][it2] +	\
			      r2*SQR(safe_features[f->frame_feature][it2]) - \
			      integral_sqr_features[f->frame_feature][t] + \
			      (1-r1)*SQR(safe_features[f->frame_feature][it1]))
#define FRAME_FEATURE(it,r) (safe_features[f->frame_feature][it]*(1-(r))+safe_features[f->frame_feature][(it)+1]*(r))
#define REGION_HIST_FEATURE() (integral_histogram_features[f->frame_feature][f->hist_bin][it2] + \
			       r2*(histogram_bins[f->frame_feature][it2]==f->hist_bin) - \
			       integral_histogram_features[f->frame_feature][f->hist_bin][t] + \
			       (1-r1)*(histogram_bins[f->frame_feature][it1]==f->hist_bin))
#define FRAME_HIST_FEATURE(it,r) ((histogram_bins[f->frame_feature][it]==f->hist_bin)*(1-(r)) + \
				  (histogram_bins[f->frame_feature][(it)+1]==f->hist_bin)*(r))
#else
#define REGION_FEATURE() REGION_FEATURE_FAST()
#define REGION_SQR_FEATURE() REGION_SQR_FEATURE_FAST()
#define FRAME_FEATURE(it,r) FRAME_FEATURE_FAST(it)
#define REGION_HIST_FEATURE() REGION_HIST_FEATURE_FAST()
#define FRAME_HIST_FEATURE(it,r) FRAME_HIST_FEATURE_FAST(it)
#endif



struct _DecisionStump;
struct _ExampleWeights;
class BehaviorBoutFeatures;
class SVMBehaviorSequence;

// When SVMBehaviorSequence is used, LABEL->data can be cast to a (BehaviorBoutSequence*)
// and SPATTERN->data can be cast to a BehaviorBoutFeatures*



/**
 * @struct _BehaviorBout
 * 
 * @brief A single labeled bout of behavior
 */
typedef struct _BehaviorBout {
  int start_frame;  /**< the starting frame number of a behavior bout in a tracked sequence */
  int end_frame;  /**< the end frame number of a behavior bout in a tracked sequence (non-inclusive).  The bout begins at frame start_frame and ends at frame end_frame-1*/
  int behavior;   /**< the index of the behavior */
  double bout_score;  /**< the score of the bout (the dot product <w_bout,f_bout>) */
  double transition_score;  /**< the component of the bout score due to transitioning from the previous behavior class to the class of this bout */
  double duration_score;  /** the compoenent of the bout score due to duration variying from the ideal duration for the behavior class of this bout */
  double unary_score;
  double loss_fn;  /**< the loss associated with missing detection of some behavior bout(s) that overlap with this bout */
  double loss_fp;  /**< the loss associated with predicting this bout incorrectly */
} BehaviorBout;


/**
 * @class BehaviorBoutSequence
 * 
 * @brief The set of all behavior bouts for a given trajectory
 */
class BehaviorBoutSequence : public StructuredLabel {
protected:
  char fname[400]; /**< The name of the file on disk from which this tracked sequence was loaded */ 
  BehaviorBoutFeatures *features;  /**< A pointer to the object used to compute all bout-level features */
  int num_bouts;  /**< the number of behavior bouts */
  BehaviorBout *bouts;  /**< A num_bouts array storing an assignment of behavior bouts to the trajectory sequence */
  double score;    /**< The score associated with the behavior predictions in bouts (the dot product <w,f>) */
  double loss;     /**< The loss associated with the behavior predictions with respect to the ground truth behavior labels */
  double slack;    /**< The error associated with the behavior predictions (score+loss-score_gt) */
  bool disable_checks;
  bool filled;

public:
  BehaviorBoutSequence(BehaviorBoutFeatures *x, SVMBehaviorSequence *svm);
  void merge_bouts_with_same_behavior();
  virtual ~BehaviorBoutSequence();
  virtual bool load(const Json::Value &x, StructuredSVM *s);
  virtual Json::Value save(StructuredSVM *s);
  virtual bool load(const char *fname) = 0;
  virtual bool save(const char *fname) { return true; };
  void Visualize(Behaviors *behaviors, const char *fname, char *html);
  virtual bool IsValid() { return num_bouts > 0; }

  friend class SVMBehaviorSequence;
};


/**
 * @class BehaviorBoutFeatures
 * 
 * @brief Features and pre-computed feature caches for an entire trajectory.  This
 * struct is used to compute and maintain bout-level features when fitting behaviors
 */
class BehaviorBoutFeatures : public StructuredData {
protected:
  double fps;
  BehaviorBoutSequence *partial_label;  /**< If non-null, contains a partial label where a subset of behavior bouts are pre-specified */
  unsigned char *memory_buffer;  /**< A malloc'd memory buffer containing dynamically allocated memory for this struct */

  /* Precomputed feature caches */
  double **features;             /**< A num_frame_features X T array of all features */
  double *frame_times;           /**< A num_frame_features array of all frame times */
  double **smoothed_features;    /**< A num_frame_features X T array of all features smoothed around some temporal window */
  double **safe_features;        /**< A num_frame_features X T array of all features (duplicate of features) */
  double **integral_features;    /**< A num_frame_features X T array encoding integral (sum) features */
  double **integral_sqr_features;/**< A num_frame_features X T array encoding integral (sum) of squared features */
  double ***feature_mins;        /**< A  num_frame_features X log2(T) X T array of the min feature response in windows of size 2^i */
  double ***feature_maxes;       /**< A  num_frame_features X log2(T) X T array of the max feature response in windows of size 2^i */
  double ***integral_histogram_features;  // A num_frame_features X num_thresholds X T encoding integral features for each histogram bin */
  int **histogram_bins;          /**< A num_frame_features X T array encoding the histogram index of each frame feature */
  //double *max_feature_responses; /**< A num_frame_features array encoding the global max response of each feature */
  //double *min_feature_responses; /**< A num_frame_features array encoding the global min response of each feature */
  //double *ave_feature_responses; /**< A num_frame_features array encoding the global ave response of each feature */

  /* Feature caches that are updated during dynamic programming */
  //double *bout_max_feature_responses; /**< A num_frame_features array encoding the max response of each feature in the current bout */
  //double *bout_min_feature_responses; /**< A num_frame_features array encoding the min response of each feature in the current bout */
  //int bout_start;  /**< The start frame of the last call to psi_bout() */
  //int bout_end;  /**< The end frame of the last call to psi_bout() */

  char fname[400];  /**< The name of the file on disk from which this tracked sequence was loaded */ 
  int num_frames;  /**< The total number of frames in this tracked sequence */

  SparseVector *fvec;  /**< segmentation-level features for this tracked (the return value of psi()) */
  int pad1, pad2;

 public:
  BehaviorBoutFeatures();
  virtual ~BehaviorBoutFeatures();
  virtual bool load(const Json::Value &x, StructuredSVM *s);
  virtual Json::Value save(StructuredSVM *s);
  virtual bool load(const char *fname, SVMBehaviorSequence *svm) = 0;

  void SetFileName(const char *fn) { strcpy(fname, fn); }
  const char *GetFileName() { return fname; }
  void AllocateBuffers(SVMBehaviorSequence *svm, bool full=true);
  void ComputeCaches(SVMBehaviorSequence *svm);
  void UpdateCaches(int t_start, int t_end, int c);
  void Clear();

  inline double ComputeRegionFeature(BoutRegionFeature *f, double t_bout_start, double t_bout_end) {
    double sum, ave, sum_sqr, t1, t2, r1, r2, m, dur, v1, v2;
    int i, it1, it2, t, d;

    if(f->is_global) {
      t1 = 0;
      t2 = num_frames;
    } else {
      if(!f->frame_coords) {
	// Normally, t_start and t_end are in a coordinate system where 0 is the start frame 
	// of the bout and 1 is the end frame of the bout
	t1 = t_bout_start + (t_bout_end-t_bout_start)*f->t_start; 
	t2 = t_bout_start + (t_bout_end-t_bout_start)*f->t_end;
      } else if(f->frame_coords == 1) {
	// t_start and t_end are assumed to be relative offsets from the start of the bout
	t1 = t_bout_start + f->t_start;
	t2 = t_bout_start + f->t_end;
      } else if(f->frame_coords == 2) {
	// t_start and t_end are assumed to be relative offsets from the start of the bout
	t1 = t_bout_end-1 + f->t_start;
	t2 = t_bout_end-1 + f->t_end;
      }
    }

#ifdef INTERPOLATE_INTEGRAL_FEATURES
    it1 = (int)t1;
    it2 = (int)t2;
    t = it1 + 1;
    r1 = t1-it1; 
    r2 = t2-it2;
    dur = t2-t1;
#else
    it1 = t = (int)(t1+.5);
    it2 = (int)(t2+.5);
    dur = it2-it1;
#endif

    if(f->hist_bin < 0) {
      // some operation on a frame feature in some region
      switch(f->op) {
      case B_SUM: 
	return REGION_FEATURE();
      case B_AVE: 
	sum = REGION_FEATURE();
	return sum/dur;
      case B_VAR: 
	sum = REGION_FEATURE();
	ave = sum/dur;
	sum_sqr = REGION_SQR_FEATURE();
	return my_max(0,sum_sqr - 2*ave*sum + SQR(ave)*dur);
      case B_DEV: 
	sum = REGION_FEATURE();
	ave = sum/dur;
	sum_sqr = REGION_SQR_FEATURE();
	return sqrt(my_max(0,sum_sqr - 2*ave*sum + SQR(ave)*dur)/dur); 
      case B_RAW:
	return FRAME_FEATURE(it2-1,r2)-FRAME_FEATURE(it1,r1);
      case B_DUR:
	return (t2-t1)/fps;
      case B_MIN:
	d = it2-t;
	i = 0;
	while(d >= (1<<(i+1))) 
	  i++;
	m = INFINITY;
	while(i >= 0) {
	  if(d & (1<<i)) {
	    m = my_min(m,feature_mins[f->frame_feature][i][t]);
	    t += 1<<i;
	  }
	  i--;
	}
#ifdef INTERPOLATE_INTEGRAL_FEATURES
	v1 = FRAME_FEATURE(it1,r1);
	v2 = FRAME_FEATURE(it2-1,r2); 
	return my_min(m,my_min(v1,v2));
#else
	return m;
#endif
      case B_MAX:
	d = it2-t;
	i = 0;
	while(d >= (1<<(i+1))) 
	  i++;
	m = -INFINITY;
	while(i >= 0) {
	  if(d & (1<<i)) {
	    m = my_max(m,feature_maxes[f->frame_feature][i][t]);
	    t += 1<<i;
	  }
	  i--;
	}
#ifdef INTERPOLATE_INTEGRAL_FEATURES
	v1 = FRAME_FEATURE(it1,r1);
	v2 = FRAME_FEATURE(it2-1,r2); 
	return my_max(m,my_max(v1,v2));
#else
	return m;
#endif
      }
    } else {
      // some operation on a histogrammed frame feature in some region
      switch(f->op) {
      case B_SUM: 
	return REGION_HIST_FEATURE();
      case B_AVE: 
	return REGION_HIST_FEATURE()/dur;
      case B_VAR: 
	sum = REGION_HIST_FEATURE();
	ave = sum/dur;
	sum_sqr = sum;
	return my_max(0,sum_sqr - 2*ave*sum + SQR(ave)/dur);
      case B_DEV: 
	sum = REGION_HIST_FEATURE();
	ave = sum/dur;
	sum_sqr = sum;
	return sqrt(my_max(0,sum_sqr - 2*ave*sum + SQR(ave)*dur)/dur); 
      case B_RAW:
	return FRAME_HIST_FEATURE(it2-1,r2)-FRAME_HIST_FEATURE(it1,r1);
      case B_DUR:
      case B_MIN:
      case B_MAX:
	fprintf(stderr, "ERROR: min, max, and dur region responses are invalid for histogram features\n");
      }
    }
    assert(0);
    return 0;
  }

  inline double ComputeFastRegionFeature(BoutRegionFeature *f, int it1, int it2, int dur, double inv_dur, double inv_fps) {
    double sum, ave, sum_sqr;

    switch(f->op) {
    case B_SUM: 
      return REGION_FEATURE_FAST();
    case B_AVE: 
      sum = REGION_FEATURE_FAST();
      return sum*inv_dur;
    case B_VAR: 
      sum = REGION_FEATURE_FAST();
      ave = sum*inv_dur;
      sum_sqr = REGION_SQR_FEATURE_FAST();
      return my_max(0,sum_sqr - 2*ave*sum + SQR(ave)*dur);
    case B_DEV: 
      sum = REGION_FEATURE_FAST();
      ave = sum*inv_dur;
      sum_sqr = REGION_SQR_FEATURE_FAST();
      return sqrt(my_max(0,sum_sqr - 2*ave*sum + SQR(ave)*dur)*inv_dur); 
    case B_RAW:
	return FRAME_FEATURE_FAST(it2-1)-FRAME_FEATURE_FAST(it1);
    case B_DUR:
      return dur*inv_fps;
    case B_MIN:
    case B_MAX:
      assert(0);
    }
    assert(0);
    return 0;
  }

  inline double ComputeBoutFeature(BoutFeature *f, double t_bout_start, double t_bout_end) {
    double r = 0, w;
    for(int i = 0; i < f->num_regions; i++) {
      w = f->weights[i];
      if(f->regions[i].is_global && f->regions[0].op == B_SUM)  // hack for certain global difference features
	w *= (t_bout_end-t_bout_start);
      r += w*ComputeRegionFeature(f->regions+i, t_bout_start, t_bout_end);
    }
    if(f->time_normalize) 
      r /= t_bout_end-t_bout_start;
    return ((f->absolute_value && r<0) ? -r : r);
  }


  friend class SVMBehaviorSequence;
  friend class BehaviorBoutSequence;
};


#define MAX_TEMPORAL_LEVELS 8
#define MAX_HARMONIC_LEVELS 8
#define MAX_BASE_FEATURES 1000


/**
 * @class SVMBehaviorSequence
 * 
 * @brief A class used for training and fitting behavior detectors/segmentors.  Typically, this class should be inherited
 * by a custom behavior class definition (e.g. SVMFlyBehaviorSequence, SVMBlobBehaviorSequence) which define domain
 * specific routines for reading in data files, features used, etc.
 *
 */
class SVMBehaviorSequence : public StructuredSVM {
 protected:
  Behaviors *behaviors;  /**< A pointer to the object defining all behavior classes */
  //int num_bouts_per_class[MAX_BEHAVIOR_GROUPS][10];  // SJB: computed in class_training_count
  double *false_negative_cost;   /**< A behaviors->num X num_classes[i] array of costs for missing detection of a given behavior class */
  double *false_positive_cost;  /**< A behaviors->num X num_classes[i] array of costs for incorrectly detecting a given behavior class */
#if USE_DURATION_COST > 0
  double *min_frame_duration;   /** shortest duration of each behavior class bout in the training set */
  double *max_frame_duration;   /** longest duration of each behavior class bout in the training set */
#endif

  int **class_training_transitions;  /**< A num_classes[i]Xclass_training_transitions_count[i][j] array of indices specifying indices of behavior classes that are allowed to proceed a given behavior class */
  int *class_training_transitions_count;  /**< A num_classes[i] array specifying the number of behavior classes that are allowed to proceed a given behavior class */
  int *class_training_count; /**< A num_classes[i] array specifying the number of times each class occurs in the training set */
  double search_all_bout_durations_up_to;
  int importance_sample_interval_size;
  double time_approximation; /**< When searching for behavior bouts, for computational purposes, the duration of bouts (in terms of # of frames) considered is a geometrically increasing series of size time_approximation,time_approximation^2,time_approximation^3...*/
  int min_bout_duration;  /**< The minimum length (in frames of a behavior bout) */
  bool **restrict_behavior_features;  /**< Untested: A num_classes[i]Xnum_features defining which bout-level features to use on a per-behavior basis.  Intended to allow different features to be used for different behaviors. */
  int max_inference_learning_frames;   // Speeds up Inference() during train-time.  Only consider starting and ending bouts at max_inference_learning_frames randomly chosen time frames

  char debugdir[400];
  bool debug_predictions, debug_weights, debug_features, debug_model;

  bool dontComputeFeatureMeanVarianceMedianStatistics;
  
  BoutFeature *bout_features; /**< The total number of bout-level features (not including class transition features) */
  int num_bout_features;
  FrameFeature *frame_features;  /**< For each frame feature, a set of parameters defining how frame-level features are expanded into bout-level features */
  int num_frame_features; /**< The total number of frame-level features */
  int feature_sample_smoothness_window;

  BoutFeature *bout_expansion_features; 
  int num_bout_expansion_features;

 public:

  Behaviors *GetBehaviors() { return behaviors; }

  void SetDebugDir(const char *dir, bool deb_pred=true, bool deb_wei=true, bool deb_feat=true, bool deb_mod=false) {
    strcpy(debugdir, dir);
    debug_predictions = deb_pred;
    debug_weights = deb_wei;
    debug_features = deb_feat;
    debug_model = deb_mod;
  }
  

  
  /**
   * @brief Constructor, assumes feature definitions will be defined later 
   *
   */
  SVMBehaviorSequence();

  /**
   * @brief Destructor
   */
  ~SVMBehaviorSequence();


  virtual bool Load(const Json::Value &root);
  virtual Json::Value Save();
  SparseVector Psi(StructuredData *x, StructuredLabel *y);
  double Inference(StructuredData *x, StructuredLabel *ybar, SparseVector *w, StructuredLabel *y_partial=NULL, StructuredLabel *y_gt=NULL, double w_scale=1);
  double Loss(StructuredLabel *y_gt, StructuredLabel *y_pred);
  void OnFinishedIteration(StructuredData *x, StructuredLabel *y, SparseVector *w=NULL, StructuredLabel *ybar=NULL);
  double ImportanceSample(StructuredData *x, SparseVector *w, StructuredLabel *y_gt, struct _SVM_cached_sample_set *set, double w_scale=1);
  void Debug();
  void DebugFeatures(const char *fname, BehaviorBoutSequence *y, double *ww);
  void DebugExample(StructuredData *x, StructuredLabel *y, StructuredLabel *ybar, SparseVector *w);
  void OnFinishedPassThroughTrainset();
  void AugmentFeatureSpace();

  virtual StructuredLabel *NewStructuredLabel(StructuredData *x) = 0;
  virtual StructuredData *NewStructuredData() = 0;

  virtual StructuredDataset *LoadDataset(const char *fname);
  virtual bool SaveDataset(StructuredDataset *d, const char *fname, int start_from);
  virtual char **load_examples(const char *fname, int *num) = 0;
  virtual void save_examples(const char *fname, StructuredDataset *dataset) = 0;

  void LoadBoutFeatureParams(const Json::Value &fe, BoutFeature *&bout_features, int &num_bout_features);
  void LoadBoutFeatureParams(const Json::Value &fe) { LoadBoutFeatureParams(fe, bout_features, num_bout_features); }
  Json::Value SaveBoutFeatureParams(BoutFeature *&bout_features, int &num_bout_features);
  void FreeBoutFeatureParams(BoutFeature *&bout_features, int &num_bout_features);
  void LoadFrameFeatureParams(const Json::Value &fe);
  Json::Value SaveFrameFeatureParams();
  void FreeFrameFeatureParams();
  bool LoadBehaviorDefinitions(const Json::Value &p);
  Json::Value SaveBehaviorDefinitions();
  void FreeBehaviorDefinitions();


  /**
   * @brief Helper function to initialize feature definitions are known before hand and passed to the constructor 
   *
   * @param num_feat The number of frame-level features
   * @param behaviors Definition of all behavior classes
   * @params sparams A num_feat array defining all bout-level features we want to use for each frame-level feature
   */
  void Init(Behaviors *behaviors, FrameFeature *frame_features = NULL, 
	    int num_frame_features = 0, BoutFeature *bout_features = NULL, int num_bout_features = 0);
  
  int NumBehaviors() { return behaviors ? behaviors->num_values : 0; }
  int NumFrameFeatures() { return num_frame_features; }
  int NumBoutFeatures() { return num_bout_features; }
  BoutFeature *GetBoutFeatureDefs() { return bout_features; }
  FrameFeature *GetFrameFeatureDefs() { return frame_features; }
  
  void print_features(FILE *fout, double *feat);
  void print_features(const char *fname, StructuredDataset *dataset, bool normalized);
  void print_weights(FILE *fout, double *w);
  void print_weights(const char *fname, double *w); 
  double loss2(StructuredLabel *y_gt,  StructuredLabel *y_pred, int debug);

  double *psi_bout(BehaviorBoutFeatures *b, int t_start, int t_end, int c, double *feat, bool normalize=true, 
		   bool fast_update=false, BoutFeature *bout_features=NULL, int num_bout_features=0);
  void saveBoutFeatures(StructuredDataset *dataset, const char *filename, bool sphered=true, bool addRandBouts=true, int window_size=32); 
  void compute_feature_mean_variance_median_statistics(StructuredDataset *dataset);
  int compute_feature_space_size();
  StructuredExample *read_struct_example(const char *label_fname, const char *features_fname, bool computeFeatures);

  // Currently, cost does not depend on the bout duration.  It only depends on the percent overlap
  // with the ground truth.  Therefore, short bouts have equal weight to long bouts with respect to  
  // the loss function.  If these were weighted linearly by duration, this would become like a per 
  // frame loss.  One could imagine doing something in the middle...
  double match_false_positive_cost(double duration, int behavior) {
    return false_positive_cost[behavior];
  }
  double match_false_negative_cost(double duration, int behavior) {
    return false_negative_cost[behavior];
  }
  void SetTimeApproximation(double a) { time_approximation = a; }

  friend class BehaviorBoutFeatures;
  friend class BehaviorBoutSequence;

  void init_bout_label(BehaviorBoutSequence *ybar, BehaviorBoutSequence *y);
 private:
  double compute_updated_bout_loss(BehaviorBoutFeatures *b, BehaviorBoutSequence *y, int T, int t_p, int t, int c_prev, double *fn, int *gt_bout, double *dur_gt, double &loss_fp, double &loss_fn);
  void update_transition_counts_with_partial_label(BehaviorBoutSequence *y_partial, int** &old_class_transitions, int* &old_class_transition_counts, int* &old_class_training_counts);
#if USE_DURATION_COST > 0
  void backtrack_optimal_solution(BehaviorBoutSequence *ybar, double **table, BehaviorBout **states, double *duration_weights, int T, int c_prev=-1, int t_p=-1);
  void sanity_check_dynamic_programming_solution(BehaviorBoutFeatures *b, BehaviorBoutSequence *ybar, BehaviorBoutSequence *y, SparseVector *w, double **class_weights, double **transition_weights, double *unary_weights, double *duration_weights, double **table, BehaviorBout **states, int T, BehaviorBoutSequence *y_partial);
#else
  void backtrack_optimal_solution(BehaviorBoutSequence *ybar, double **table, BehaviorBout **states, int T, int c_prev=-1, int t_p=-1);
  void sanity_check_dynamic_programming_solution(BehaviorBoutFeatures *b, BehaviorBoutSequence *ybar, BehaviorBoutSequence *y, SparseVector *w, double **class_weights, double **transition_weights, double *unary_weights, double **table, BehaviorBout **states, int T, BehaviorBoutSequence *y_partial);
#endif
  bool check_agreement_with_partial_label(BehaviorBoutSequence *y_partial, int t_p, int t, int *partial_label_bout, int &restrict_c_prev);
  void store_solution(BehaviorBout &state, int t_p, int t, int c_prev, double bout_score, double transition_score, double unary_score, double loss_fn, double loss_fp, double duration_score);
  void restore_transition_counts(BehaviorBoutSequence *y_partial, int** &old_class_transitions, int* &old_class_transition_counts, int* &old_class_training_counts);
  bool *get_allowable_frame_times(BehaviorBoutSequence *y_gt, BehaviorBoutSequence *y_partial, int T);
  int get_bout_start_time(int *duration, int &tt, int t_p, int t, int &next_duration, int &last_gt, int &last_partial, int *gt_bout, int *partial_label_bout, BehaviorBoutSequence *y, BehaviorBoutSequence *y_partial, int &restrict_c_prev, int &restrict_c_next, bool *allowable_time_frames);
  BehaviorBoutSequence *bout_sequence_remove_section(BehaviorBoutSequence *y_src, int t_start, int t_end);
  void print_bout_sequence_scores(BehaviorBoutSequence *y);
  void init_weight_pointers(double* &ptr, double** &class_weights, double** &transition_weights, double* &unary_weights, double* &duration_weights);
  double get_duration_score(int duration, int c, double *duration_weights);
  bool fill_unlabeled_gt_frames(BehaviorBoutSequence *&y_gt, BehaviorBoutSequence *&y_partial);
  void ComputeClassTrainsitionCounts();

  struct _DecisionStump *SelectBestFeature(struct _ExampleWeights *samples, int num_samples, double **f);
  void AppendNewFeature(struct _DecisionStump *s, struct _ExampleWeights *samples, int num_samples);
  struct _ExampleWeights *ComputeExampleWeights(StructuredDataset *trainset, int *numEx);
};

void free_behavior_bout_sequence(BehaviorBoutSequence *b, int num);

#endif
