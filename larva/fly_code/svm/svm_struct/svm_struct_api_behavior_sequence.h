#ifndef  __SVM_STRUCT_API_BEHAVIOR_SEQUENCE__
#define __SVM_STRUCT_API_BEHAVIOR_SEQUENCE__

#include "svm_struct_api_types.h"
#include "svm_struct_api.h"


#include "../../blob.h"


#define FORMAT__BOUT_FEATURE_PARAMS "feature_sample_smoothness_window=%d, num_temporal_levels=%d, num_bout_max_thresholds=%d, "\
		     "num_bout_min_thresholds=%d, num_bout_change_points=%d, num_histogram_bins=%d, "\
		     "num_histogram_temporal_levels=%d, num_difference_temporal_levels=%d, num_harmonic_features=%d, "\
		     "use_bout_sum_features=%d, use_bout_ave_features=%d, use_bout_sum_absolute_features=%d, use_bout_ave_absolute_features=%d, use_standard_deviation=%d, use_sum_variance=%d, "\
		     "use_bout_max_feature=%d, use_bout_min_feature=%d, "\
		     "use_global_difference_max_ave_features=%d, use_global_difference_min_ave_features=%d, use_global_difference_ave_ave_features=%d, "\
		     "use_global_difference_max_sum_features=%d, use_global_difference_min_sum_features=%d, use_global_difference_ave_sum_features=%d, "\
		     "use_bout_change=%d, use_bout_absolute_change=%d, use_histogram_sum_features=%d, "\
		     "use_histogram_ave_features=%d, use_sum_harmonic_features=%d, use_ave_harmonic_features=%d, use_sum_absolute_harmonic_features=%d, "\
		     "use_ave_absolute_harmonic_features=%d, use_start_sum_absolute_diff_haar_features=%d, "\
		     "use_end_sum_absolute_diff_haar_features=%d, use_start_sum_diff_haar_features=%d, use_end_sum_diff_haar_features=%d,  "\
		     "use_start_ave_absolute_diff_haar_features=%d, use_end_ave_absolute_diff_haar_features=%d, use_start_ave_diff_haar_features=%d,  "\
		     "use_end_ave_diff_haar_features=%d%*[^\n]\n"


struct _BehaviorGroups;
struct _BehaviorBoutFeatures;

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
} BehaviorBout;


/**
 * @struct _BehaviorBoutSequence
 * 
 * @brief The set of all behavior bouts for a given trajectory
 */
typedef struct _BehaviorBoutSequence {
  struct _BehaviorGroups *behaviors;  /**< A pointer to the object defining all behavior classes */
  struct _BehaviorBoutFeatures *features;  /**< A pointer to the object used to compute all bout-level features */
  int *num_bouts;  /**< A behaviors->num array storing the number or behavior bouts for each grouping of behaviors */
  BehaviorBout **bouts;  /**< A behaviors->numXnum_bouts[i] array storing an assignment of behavior bouts to the trajectory sequence */
  double score;  /**< The score associated with the behavior predictions in bouts (the dot product <w,f>) */
  double loss;  /**< The loss associated with the behavior predictions with respect to the ground truth behavior labels */
  double slack;  /**< The error associated with the behavior predictions (score+loss-score_gt) */
  double *scores;  /**< A behaviors->num array of scores for each behavior group */
  double *losses;  /**< A behaviors->num array of losses for each behavior group */
} BehaviorBoutSequence;


/**
 * @struct _BehaviorBoutFeatures
 * 
 * @brief Features and pre-computed feature caches for an entire trajectory.  This
 * struct is used to compute and maintain bout-level features when fitting behaviors
 */
typedef struct _BehaviorBoutFeatures {
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

  void *data;  /**< A pointer to a structure defining all data for a particular tracked sequence (the return value of load_training_example()) */
  char fname[400];  /**< The name of the file on disk from which this tracked sequence was loaded */ 
  int num_frames;  /**< The total number of frames in this tracked sequence */

  SVECTOR *fvec;  /**< segmentation-level features for this tracked (the return value of psi()) */
} BehaviorBoutFeatures;


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
class SVMBehaviorSequence : public SVMStructMethod {
 protected:
  struct _BehaviorGroups *behaviors;  /**< A pointer to the object defining all behavior classes */
  int num_classes[MAX_BEHAVIOR_GROUPS];  /**< The number of behavior classes for each group (same as behaviors->num_values) */
  int num_features;  /**< The total number of bout-level features (not including class transition features) */
  int num_base_features;  /**< The total number of frame-level features */
  int feature_diff_frames;  /**< DEPRECATED: I don't think this is used anymore; however deleting it might change the file format for the learned output files */
  int sizePsi;   /**< The size of the feature space (like num_features but includes class transition features) */
  int behavior;   /**< If not -1, only run the model on one behavior group */
  double *false_negative_cost[MAX_BEHAVIOR_GROUPS];   /**< A behaviors->num X num_classes[i] array of costs for missing detection of a given behavior class */
  double *false_positive_cost[MAX_BEHAVIOR_GROUPS];  /**< A behaviors->num X num_classes[i] array of costs for incorrectly detecting a given behavior class */
  int ***class_training_transitions;  /**< A behaviors->numXnum_classes[i]Xclass_training_transitions_count[i][j] array of indices specifying indices of behavior classes that are allowed to proceed a given behavior class */
  int **class_training_transitions_count;  /**< A behaviors->numXnum_classes[i] array specifying the number of behavior classes that are allowed to proceed a given behavior class */
  int **class_training_count; /**< A behaviors->numXnum_classes[i] array specifying the number of times each class occurs in the training set */
  double time_approximation; /**< When searching for behavior bouts, for computational purposes, the duration of bouts (in terms of # of frames) considered is a geometrically increasing series of size time_approximation,time_approximation^2,time_approximation^3...*/
  SVMFeatureParams feature_params[MAX_BASE_FEATURES];  /**< For each frame feature, a set of parameters defining how frame-level features are expanded into bout-level features */
  double *features_mu;  /**< A num_features array defining the mean of each bout-level feature (used to normalize all features to be roughly on the same scale) */
  double *features_gamma;  /**< A num_features array defining the inverse of the standard deviation of each bout-level feature (used to normalize all features to be roughly on the same scale) */
  double **histogram_thresholds;  /**< A num_base_featuresXfeature_params[i]->num_histogram_bins defining a set of thresholds for histogram bins, such that each bin in the histogram is between two adjacent thresholds */
  double **min_thresholds;  /**< A num_base_featuresXfeature_params[i]->num_histogram_bins defining a set of thresholds for checking if the min of a given feature is below a given threshold */
  double **max_thresholds;  /**< A num_base_featuresXfeature_params[i]->num_histogram_bins defining a set of thresholds for checking if the max of a given feature is below a given threshold */
  int min_bout_duration;  /**< The minimum length (in frames of a behavior bout) */
  char **feature_names;  /**< A num_features array of strings defining a human-interpretable name for each feature */
  bool **restrict_behavior_features[MAX_BEHAVIOR_GROUPS];  /**< Untested: A behaviors->numXnum_classes[i]Xnum_features defining which bout-level features to use on a per-behavior basis.  Intended to allow different features to be used for different behaviors. */
  
 public:
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

  /**
   * @brief Helper function to initialize feature definitions are known before hand and passed to the constructor 
   *
   * @param num_feat The number of frame-level features
   * @param behaviors Definition of all behavior classes
   * @param beh If -1 uses all groups of behaviors, otherwise uses just one group behaviors->behaviors[beh] (a group is a set of mutually exclusive behavior classes)
   * @params sparams A num_feat array defining all bout-level features we want to use for each frame-level feature
   */
  void Init(int num_feat, struct _BehaviorGroups *behaviors, int beh, SVMFeatureParams *sparams = NULL);
  
  /**
   * @brief Train a new behavior detector/segmentor.  This is a wrapper for the standalone command line program to train behaviors
   */
  int train (int argc, const char* argv[], STRUCT_LEARN_PARM *struct_parm, STRUCTMODEL *structmodel);
  
  
  /**
   * @brief Test a new behavior detector/segmentor.  This is a wrapper for the standalone command line program to test behavior detectors
   */
  int test (int argc, const char* argv[], STRUCT_LEARN_PARM *struct_parm, STRUCTMODEL *structmodel);

  SVMFeatureParams DefaultParams();
  void compute_bout_feature_caches(BehaviorBoutFeatures *blob_features);
  int compute_feature_space_size();
  void update_bout_feature_caches(BehaviorBoutFeatures *b, int t_start, int t_end, int c);

  bool ReadFeatureParam(FILE *modelfl, SVMFeatureParams *p);
  void print_features(FILE *fout, double *feat);
  void print_features(const char *fname, EXAMPLE *ex, int n, bool normalized);
  void print_weights(FILE *fout, double *w);
  void print_weights(const char *fname, double *w);
  void set_feature_name(int feature_ind, int base_feature_ind, const char *name);

  void on_finished_iteration(CONSTSET c, STRUCTMODEL *sm, 
						STRUCT_LEARN_PARM *sparm, int iter_num) ;
  void on_finished_find_most_violated_constraint(LABEL *ybar, LABEL *y, int iter, STRUCT_LEARN_PARM *sparm, const char *ename);


  BehaviorBoutFeatures *create_behavior_bout_feature_cache(void *d, bool compute_cache=true);
  EXAMPLE *find_example(SAMPLE s, const char *fname);


  void        svm_struct_learn_api_init(int argc, const char* argv[]);
  void        svm_struct_learn_api_exit();
  void        svm_struct_classify_api_init(int argc, const char* argv[]);
  void        svm_struct_classify_api_exit();
  SAMPLE      read_struct_examples(char *file, STRUCT_LEARN_PARM *sparm);
  void        init_struct_model(SAMPLE sample, STRUCTMODEL *sm, 
			      STRUCT_LEARN_PARM *sparm, LEARN_PARM *lparm, 
			      KERNEL_PARM *kparm);
  CONSTSET    init_struct_constraints(SAMPLE sample, STRUCTMODEL *sm, 
				    STRUCT_LEARN_PARM *sparm);
  LABEL       find_most_violated_constraint_slackrescaling(SPATTERN x, LABEL y, 
						     STRUCTMODEL *sm, 
						     STRUCT_LEARN_PARM *sparm);
  LABEL       find_most_violated_constraint_marginrescaling(SPATTERN x, LABEL y, 
						     STRUCTMODEL *sm, 
						     STRUCT_LEARN_PARM *sparm);
  LABEL       classify_struct_example(SPATTERN x, STRUCTMODEL *sm, 
				    STRUCT_LEARN_PARM *sparm);
  int         empty_label(LABEL y);
  SVECTOR     *psi(SPATTERN x, LABEL y, STRUCTMODEL *sm, 
			  STRUCT_LEARN_PARM *sparm);
  double      loss2(LABEL y, LABEL ybar, STRUCT_LEARN_PARM *sparm, int beh, int debug=0);
  double      loss(LABEL y, LABEL ybar, STRUCT_LEARN_PARM *sparm);
  int         finalize_iteration(double ceps, int cached_constraint,
			       SAMPLE sample, STRUCTMODEL *sm,
			       CONSTSET cset, double *alpha, 
			       STRUCT_LEARN_PARM *sparm);
  void        print_struct_learning_stats(SAMPLE sample, STRUCTMODEL *sm,
					CONSTSET cset, double *alpha, 
					STRUCT_LEARN_PARM *sparm);
  void        print_struct_testing_stats(SAMPLE sample, STRUCTMODEL *sm,
				       STRUCT_LEARN_PARM *sparm,
				       STRUCT_TEST_STATS *teststats);
  void        eval_prediction(long exnum, EXAMPLE ex, LABEL prediction, 
			    STRUCTMODEL *sm, STRUCT_LEARN_PARM *sparm,
			    STRUCT_TEST_STATS *teststats);
  void        write_struct_model(const char *file,STRUCTMODEL *sm, 
			       STRUCT_LEARN_PARM *sparm);
  STRUCTMODEL read_struct_model(const char *file, STRUCT_LEARN_PARM *sparm);
  void        write_label(FILE *fp, LABEL y);
  void        free_pattern(SPATTERN x);
  void        free_label(LABEL y);
  void        free_struct_model(STRUCTMODEL sm);
  void        free_struct_sample(SAMPLE s);
  void        print_struct_help();
  void        parse_struct_parameters(STRUCT_LEARN_PARM *sparm);
  void        print_struct_help_classify();
  void        parse_struct_parameters_classify(STRUCT_LEARN_PARM *sparm);

  BehaviorBoutSequence *read_bout_sequence(FILE *fin, char *fname, int *iter, int *num);
  void write_bout_sequence(BehaviorBoutSequence *b, FILE *fout, int iter, int num);
  LABEL read_label(FILE *fin, char *fname, int *iter, int *num);
  void write_label(LABEL y, FILE *fout, int iter, int num);
  
  double      *psi_bout(BehaviorBoutFeatures *b, int t_start, int t_end, int beh, int c, double *feat, bool normalize=true, bool fast_update=false);
  LABEL       inference_via_dynamic_programming(SPATTERN *x, STRUCTMODEL *sm, 
						STRUCT_LEARN_PARM *sparm, LABEL *yy);
  void compute_feature_mean_variance_median_statistics(EXAMPLE *ex, int num_examples);

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


  /******** Functions that should be overridden by a child class that inherits this class **************/

  /**
   * @brief Read a dataset file containing a list of training files.  Each training file is some tracked sequence that may contain multiple bouts of behaviors.
   *
   * @param fname The name of the dataset file
   * @param num The returned number of training files (a pointer set by this function)
   * @return A *num array, a list of all training file names
   */
  virtual char **load_examples(const char *fname, int *num) = 0;

  /**
   * @brief Read a particular training file.  Each training file is some tracked sequence that may contain multiple bouts of behaviors.
   *
   * @param fname The name of the training file
   * @return A pointer to an object (of custom type) that contains all loaded data applicaple to this training example
   */
  virtual void *load_training_example(const char *fname) = 0;

  /**
   * @brief Get a human-interpretable name of the ith frame-level feature
   *
   * @param i The index of the frame-level feature
   */
  virtual const char *get_base_feature_name(int i) = 0;

  virtual void load_from_bout_sequence(BehaviorBoutSequence *y, void *b) = 0;
  virtual BehaviorBoutSequence *create_behavior_bout_sequence(void *b, BehaviorGroups *behaviors, bool build_partial_label) = 0;
  virtual void load_behavior_bout_features(void *b, BehaviorBoutFeatures *feature_cache) = 0; 
  virtual void free_data(void *d) = 0;
  virtual int num_frames(void *d) = 0;
};

void free_behavior_bout_sequence(BehaviorBoutSequence *b, int num);

#endif
