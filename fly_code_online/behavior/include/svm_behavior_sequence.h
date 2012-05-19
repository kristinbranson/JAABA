
#if DEBUG > 0 
extern char *g_currFile; // CSC 20110420: hack to pass current filename for debug purposes
#endif

#ifndef  __SVM_STRUCT_API_BEHAVIOR_SEQUENCE__
#define __SVM_STRUCT_API_BEHAVIOR_SEQUENCE__



#include "blob.h"
#include "structured_svm.h"



#define ALLOW_SAME_TRANSITIONS
#define DEBUG 1
#define MAX_FILENAME 1000
#define MAX_FEATURES 1000

#if DEBUG > 0
//#define MAX_FEATURES 1000
extern const char *bout_feature_names[]; // initialized in svm_struct_api_behavior_sequence.cpp; used to display feature names (in combination with g_feature_map)
//extern const char *g_feature_names[MAX_FEATURES];
#endif

#define NUMFEAT 20

#define FEATURE_SAMPLE_SMOOTHNESS_WINDOW 1
#define NUM_TEMPORAL_LEVELS 2
#define NUM_BOUT_MAX_THRESHOLDS 3
#define NUM_BOUT_MIN_THRESHOLDS 4
#define NUM_BOUT_CHANGE_POINTS 5
#define NUM_HISTOGRAM_BINS 6
#define NUM_HISTOGRAM_TEMPORAL_LEVELS 7
#define NUM_DIFFERENCE_TEMPORAL_LEVELS 8
#define NUM_HARMONIC_FEATURES 9
#define USE_BOUT_SUM_FEATURES 10
#define USE_BOUT_AVE_FEATURES 11
#define USE_BOUT_SUM_ABSOLUTE_FEATURES 12
#define USE_BOUT_AVE_ABSOLUTE_FEATURES 13
#define USE_STANDARD_DEVIATION 14
#define USE_SUM_VARIANCE 15
#define USE_BOUT_MAX_FEATURE 16
#define USE_BOUT_MIN_FEATURE 17
#define USE_GLOBAL_DIFFERENCE_MAX_AVE_FEATURES 18
#define USE_GLOBAL_DIFFERENCE_MIN_AVE_FEATURES 19
#define USE_GLOBAL_DIFFERENCE_AVE_AVE_FEATURES 20 
#define USE_GLOBAL_DIFFERENCE_MAX_SUM_FEATURES 21
#define USE_GLOBAL_DIFFERENCE_MIN_SUM_FEATURES 22
#define USE_GLOBAL_DIFFERENCE_AVE_SUM_FEATURES 23
#define USE_BOUT_CHANGE 24
#define USE_BOUT_ABSOLUTE_CHANGE 25
#define USE_HISTOGRAM_SUM_FEATURES 26
#define USE_HISTOGRAM_AVE_FEATURES 27
#define USE_SUM_HARMONIC_FEATURES 28
#define USE_AVE_HARMONIC_FEATURES 29
#define USE_SUM_ABSOLUTE_HARMONIC_FEATURES 30
#define USE_AVE_ABSOLUTE_HARMONIC_FEATURES 31
#define USE_START_SUM_ABSOLUTE_DIFF_HAAR_FEATURES 32
#define USE_END_SUM_ABSOLUTE_DIFF_HAAR_FEATURES 33
#define USE_START_SUM_DIFF_HAAR_FEATURES 34
#define USE_END_SUM_DIFF_HAAR_FEATURES 35
#define USE_START_AVE_ABSOLUTE_DIFF_HAAR_FEATURES 36
#define USE_END_AVE_ABSOLUTE_DIFF_HAAR_FEATURES 37
#define USE_START_AVE_DIFF_HAAR_FEATURES 38
#define USE_END_AVE_DIFF_HAAR_FEATURES 39
#define USE_UNARY_FEATURE 40
#define USE_BOUT_DURATIONS 41


struct _BehaviorGroups;
class BehaviorBoutFeatures;
class SVMBehaviorSequence;

// When SVMBehaviorSequence is used, LABEL->data can be cast to a (BehaviorBoutSequence*)
// and SPATTERN->data can be cast to a BehaviorBoutFeatures*

/**
 * @struct _SVMFeatureParams
 * 
 * @brief Parameters defining how to construct bout-level features from frame-level
 * features
 *
 * TODO: comment this
 */
typedef struct _SVMFeatureParams {
  int feature_sample_smoothness_window, num_temporal_levels, num_bout_max_thresholds, num_bout_min_thresholds,
    num_bout_change_points, num_histogram_bins, num_histogram_temporal_levels, num_difference_temporal_levels, num_harmonic_features;   
  bool use_bout_sum_features, use_bout_ave_features, use_standard_deviation, use_sum_variance, use_bout_max_feature, use_bout_min_feature,
    use_global_difference_max_ave_features, use_global_difference_min_ave_features, use_global_difference_ave_ave_features,
    use_global_difference_max_sum_features, use_global_difference_min_sum_features, use_global_difference_ave_sum_features,
    use_bout_change, use_bout_absolute_change, use_histogram_sum_features,
    use_histogram_ave_features, use_sum_harmonic_features, use_ave_harmonic_features, use_sum_absolute_harmonic_features,
    use_ave_absolute_harmonic_features, use_start_sum_absolute_diff_haar_features,
    use_end_sum_absolute_diff_haar_features, use_start_sum_diff_haar_features, use_end_sum_diff_haar_features, 
    use_start_ave_absolute_diff_haar_features, use_end_ave_absolute_diff_haar_features, use_start_ave_diff_haar_features, 
    use_end_ave_diff_haar_features;
  bool use_bout_sum_absolute_features, use_bout_ave_absolute_features;

  int num_features;
} SVMFeatureParams;


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
  double loss_fn;  /**< the loss associated with missing detection of some behavior bout(s) that overlap with this bout */
  double loss_fp;  /**< the loss associated with predicting this bout incorrectly */
  double extreme_vals[2][NUMFEAT];
} BehaviorBout;


/**
 * @class BehaviorBoutSequence
 * 
 * @brief The set of all behavior bouts for a given trajectory
 */
class BehaviorBoutSequence : public StructuredLabel {
protected:
  char fname[400];  /**< The name of the file on disk from which this tracked sequence was loaded */ 
  struct _BehaviorGroups *behaviors;  /**< A pointer to the object defining all behavior classes */
  BehaviorBoutFeatures *features;  /**< A pointer to the object used to compute all bout-level features */
  int *num_bouts;  /**< A behaviors->num array storing the number or behavior bouts for each grouping of behaviors */
  BehaviorBout **bouts;  /**< A behaviors->numXnum_bouts[i] array storing an assignment of behavior bouts to the trajectory sequence */
  double score;  /**< The score associated with the behavior predictions in bouts (the dot product <w,f>) */
  double loss;  /**< The loss associated with the behavior predictions with respect to the ground truth behavior labels */
  double slack;  /**< The error associated with the behavior predictions (score+loss-score_gt) */
  double *scores;  /**< A behaviors->num array of scores for each behavior group */
  double *losses;  /**< A behaviors->num array of losses for each behavior group */

public:
  BehaviorBoutSequence(BehaviorBoutFeatures *x, SVMBehaviorSequence *svm);
  virtual ~BehaviorBoutSequence();
  virtual bool load(const Json::Value &x, StructuredSVM *s);
  virtual Json::Value save(StructuredSVM *s);
  virtual bool load(const char *fname) = 0;
  virtual bool save(const char *fname) { return true; };
  void Visualize(BehaviorGroups *groups, int beh, const char *fname, char *html);

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
  BehaviorBoutSequence *partial_label;  /**< If non-null, contains a partial label where a subset of behavior bouts are pre-specified */
  unsigned char *memory_buffer;  /**< A malloc'd memory buffer containing dynamically allocated memory for this struct */

  /* Precomputed feature caches */
  double **features;             /**< A num_base_features X T array of all features */
  double *frame_times;           /**< A num_base_features array of all frame times */
  double **smoothed_features;    /**< A num_base_features X T array of all features smoothed around some temporal window */
  double **integral_features;    /**< A num_base_features X T array encoding integral (sum) features */
  double **integral_sqr_features;/**< A num_base_features X T array encoding integral (sum) of squared features */
  double *max_feature_responses; /**< A num_base_features array encoding the global max response of each feature */
  double *min_feature_responses; /**< A num_base_features array encoding the global min response of each feature */
  double *ave_feature_responses; /**< A num_base_features array encoding the global ave response of each feature */
  double ***integral_histogram_features;  // A num_base_features X num_thresholds X T encoding integral features for each histogram bin */
  int **histogram_bins;          /**< A num_base_features X T array encoding the histogram index of each frame feature */

  /* Feature caches that are updated during dynamic programming */
  double *bout_max_feature_responses; /**< A num_base_features array encoding the max response of each feature in the current bout */
  double *bout_min_feature_responses; /**< A num_base_features array encoding the min response of each feature in the current bout */
  int bout_start;  /**< The start frame of the last call to psi_bout() */
  int bout_end;  /**< The end frame of the last call to psi_bout() */

  char fname[400];  /**< The name of the file on disk from which this tracked sequence was loaded */ 
  int num_frames;  /**< The total number of frames in this tracked sequence */
  int num_base_features;

  SparseVector *fvec;  /**< segmentation-level features for this tracked (the return value of psi()) */

 public:
  BehaviorBoutFeatures();
  virtual ~BehaviorBoutFeatures();
  virtual bool load(const Json::Value &x, StructuredSVM *s);
  virtual Json::Value save(StructuredSVM *s);
  virtual bool load(const char *fname, SVMBehaviorSequence *svm, BehaviorBoutSequence *y) = 0;

  void SetFileName(const char *fn) { strcpy(fname, fn); }
  const char *GetFileName() { return fname; }
  void AllocateBuffers(SVMBehaviorSequence *svm);
  void ComputeCaches(SVMBehaviorSequence *svm);
  void UpdateCaches(int t_start, int t_end, int c);

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
  struct _BehaviorGroups *behaviors;  /**< A pointer to the object defining all behavior classes */
  int num_classes[MAX_BEHAVIOR_GROUPS];  /**< The number of behavior classes for each group (same as behaviors->num_values) */
  int num_features;  /**< The total number of bout-level features (not including class transition features) */
  int num_base_features;  /**< The total number of frame-level features */
  int feature_diff_frames;  /**< DEPRECATED: I don't think this is used anymore; however deleting it might change the file format for the learned output files */
  int behavior;   /**< If not -1, only run the model on one behavior group */
  double *false_negative_cost[MAX_BEHAVIOR_GROUPS];   /**< A behaviors->num X num_classes[i] array of costs for missing detection of a given behavior class */
  double *false_positive_cost[MAX_BEHAVIOR_GROUPS];  /**< A behaviors->num X num_classes[i] array of costs for incorrectly detecting a given behavior class */
  int ***class_training_transitions;  /**< A behaviors->numXnum_classes[i]Xclass_training_transitions_count[i][j] array of indices specifying indices of behavior classes that are allowed to proceed a given behavior class */
  int **class_training_transitions_count;  /**< A behaviors->numXnum_classes[i] array specifying the number of behavior classes that are allowed to proceed a given behavior class */
  int **class_training_count; /**< A behaviors->numXnum_classes[i] array specifying the number of times each class occurs in the training set */
  double search_all_bout_durations_up_to;
  double time_approximation; /**< When searching for behavior bouts, for computational purposes, the duration of bouts (in terms of # of frames) considered is a geometrically increasing series of size time_approximation,time_approximation^2,time_approximation^3...*/
  SVMFeatureParams feature_params[MAX_BASE_FEATURES];  /**< For each frame feature, a set of parameters defining how frame-level features are expanded into bout-level features */
  double *features_mu;  /**< A num_features array defining the mean of each bout-level feature (used to normalize all features to be roughly on the same scale) */
  double *features_gamma;  /**< A num_features array defining the inverse of the standard deviation of each bout-level feature (used to normalize all features to be roughly on the same scale) */
  double **histogram_thresholds;  /**< A num_base_featuresXfeature_params[i]->num_histogram_bins defining a set of thresholds for histogram bins, such that each bin in the histogram is between two adjacent thresholds */
  double **min_thresholds;  /**< A num_base_featuresXfeature_params[i]->num_histogram_bins defining a set of thresholds for checking if the min of a given feature is below a given threshold */
  double **max_thresholds;  /**< A num_base_featuresXfeature_params[i]->num_histogram_bins defining a set of thresholds for checking if the max of a given feature is below a given threshold */
  int min_bout_duration;  /**< The minimum length (in frames of a behavior bout) */
  bool **restrict_behavior_features[MAX_BEHAVIOR_GROUPS];  /**< Untested: A behaviors->numXnum_classes[i]Xnum_features defining which bout-level features to use on a per-behavior basis.  Intended to allow different features to be used for different behaviors. */
  int max_inference_learning_frames;   // Speeds up Inference() during train-time.  Only consider starting and ending bouts at max_inference_learning_frames randomly chosen time frames

  char debugdir[400];
  bool debug_predictions, debug_weights, debug_features, debug_model;

 public:
  char **feature_names;  /**< A num_features array of strings defining a human-interpretable name for each feature */
  char base_feature_names[MAX_BASE_FEATURES][1001];

  /**
   * @brief Constructor, assumes feature definitions are known before hand and passed to the constructor 
   *
   * @param num_feat The number of frame-level features
   * @param behaviors Definition of all behavior classes
   * @param beh If -1 uses all groups of behaviors, otherwise uses just one group behaviors->behaviors[beh] (a group is a set of mutually exclusive behavior classes)
   * @params sparams A num_feat array defining all bout-level features we want to use for each frame-level feature
   */
  SVMBehaviorSequence(int num_feat, struct _BehaviorGroups *behaviors, int beh, SVMFeatureParams *sparams = NULL);
  
  /**
   * @brief Constructor, assumes feature definitions will be defined later 
   *
   * @param behaviors Definition of all behavior classes
   * @param beh If -1 uses all groups of behaviors, otherwise uses just one group behaviors->behaviors[beh] (a group is a set of mutually exclusive behavior classes)
   */
  SVMBehaviorSequence(struct _BehaviorGroups *behaviors, int beh);

  /**
   * @brief Destructor
   */
  ~SVMBehaviorSequence();


  virtual bool Load(const Json::Value &root);
  virtual Json::Value Save();
  SparseVector Psi(StructuredData *x, StructuredLabel *y);
  double Inference(StructuredData *x, StructuredLabel *ybar, SparseVector *w, StructuredLabel *y_partial=NULL, StructuredLabel *y_gt=NULL, double w_scale=1);
  double Loss(StructuredLabel *y_gt, StructuredLabel *y_pred);
  void OnFinishedIteration(StructuredData *x, StructuredLabel *y, StructuredLabel *ybar);

  virtual StructuredLabel *NewStructuredLabel(StructuredData *x) = 0;
  virtual StructuredData *NewStructuredData() = 0;

  virtual StructuredDataset *LoadDataset(const char *fname);
  virtual bool SaveDataset(StructuredDataset *d, const char *fname, int start_from);
  virtual const char *get_base_feature_name(int ind) = 0;
  virtual char **load_examples(const char *fname, int *num) = 0;
  virtual void save_examples(const char *fname, StructuredDataset *dataset) = 0;


  /**
   * @brief Helper function to initialize feature definitions are known before hand and passed to the constructor 
   *
   * @param num_feat The number of frame-level features
   * @param behaviors Definition of all behavior classes
   * @param beh If -1 uses all groups of behaviors, otherwise uses just one group behaviors->behaviors[beh] (a group is a set of mutually exclusive behavior classes)
   * @params sparams A num_feat array defining all bout-level features we want to use for each frame-level feature
   */
  void Init(int num_feat, struct _BehaviorGroups *behaviors, int beh, SVMFeatureParams *sparams = NULL);
  
  int NumBaseFeatures() { return num_base_features; }

  SVMFeatureParams DefaultParams();
  struct _BehaviorGroups *Behaviors() { return behaviors; }
  
  void print_features(FILE *fout, double *feat);
  void print_features(const char *fname, StructuredDataset *dataset, bool normalized);
  void print_weights(FILE *fout, double *w);
  void print_weights(const char *fname, double *w); 
  void set_feature_name(int feature_ind, int base_feature_ind, const char *name);
  double loss2(StructuredLabel *y_gt,  StructuredLabel *y_pred, int beh, int debug);


  double *psi_bout(BehaviorBoutFeatures *b, int t_start, int t_end, int beh, int c, double *feat, bool normalize=true, 
		   bool fast_update=false, double get_extreme_vals[2][NUMFEAT]=NULL, double use_extreme_vals[2][NUMFEAT]=NULL);
  void saveBoutFeatures(StructuredDataset *dataset, const char *filename, bool sphered=true, bool addRandBouts=true); 
  void compute_feature_mean_variance_median_statistics(StructuredDataset *dataset);
  int compute_feature_space_size();
  StructuredExample *read_struct_example(const char *label_fname, const char *features_fname, bool computeFeatures);

  // Currently, cost does not depend on the bout duration.  It only depends on the percent overlap
  // with the ground truth.  Therefore, short bouts have equal weight to long bouts with respect to  
  // the loss function.  If these were weighted linearly by duration, this would become like a per 
  // frame loss.  One could imagine doing something in the middle...
  double match_false_positive_cost(double duration, int group, int behavior) {
    return false_positive_cost[group][behavior];
  }
  double match_false_negative_cost(double duration, int group, int behavior) {
    return false_negative_cost[group][behavior];
  }

  friend class BehaviorBoutFeatures;
  friend class BehaviorBoutSequence;

 private:
  void init_bout_label(BehaviorBoutSequence *ybar, BehaviorBoutSequence *y);
  double compute_updated_bout_loss(BehaviorBoutFeatures *b, BehaviorBoutSequence *y, int beh, int T, int t_p, int t, int c_prev, double *fn, int *gt_bout, double *dur_gt, double &loss_fp, double &loss_fn);
  void update_transition_counts_with_partial_label(int beh, BehaviorBoutSequence *y_partial, int* &old_class_transition_counts, int* &old_class_training_counts);
  void backtrack_optimal_solution(BehaviorBoutSequence *ybar, int beh, double **table, BehaviorBout **states, double *unary_weights, int T);
  bool check_agreement_with_partial_label(BehaviorBoutSequence *y_partial, int beh, int t_p, int t, int *partial_label_bout, int &restrict_c_prev);
  void store_solution(BehaviorBout &state, int t_p, int t, int c_prev, double bout_score, double transition_score, double loss_fn, double loss_fp, double extreme_vals[2][NUMFEAT]);
  void restore_transition_counts(int beh, BehaviorBoutSequence *y_partial, int* &old_class_transition_counts, int* &old_class_training_counts);
  void sanity_check_dynamic_programming_solution(int beh, BehaviorBoutFeatures *b, BehaviorBoutSequence *ybar, BehaviorBoutSequence *y, SparseVector *w, double **class_weights, double **transition_weights, double *unary_weights, double **table, BehaviorBout **states, int T);
  bool *get_allowable_frame_times(BehaviorBoutSequence *y_gt, BehaviorBoutSequence *y_partial, int T);
  int get_bout_start_time(int beh, int *duration, int &tt, int t_p, int t, int &next_duration, int &last_gt, int &last_partial, int *gt_bout, int *partial_label_bout, BehaviorBoutSequence *y, BehaviorBoutSequence *y_partial, int &restrict_c_prev, int &restrict_c_next);
};

void free_behavior_bout_sequence(BehaviorBoutSequence *b, int num);

#endif
