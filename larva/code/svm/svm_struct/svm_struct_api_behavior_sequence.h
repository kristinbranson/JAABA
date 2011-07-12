#ifndef  __SVM_STRUCT_API_BEHAVIOR_SEQUENCE__
#define __SVM_STRUCT_API_BEHAVIOR_SEQUENCE__

#include "svm_struct_api_types.h"
#include "svm_struct_api.h"


#include "../../blob.h"

struct _BehaviorGroups;
struct _BehaviorBoutFeatures;

// When SVMBehaviorSequence is used, LABEL->data can be cast to a (BehaviorBoutSequence*)
// and SPATTERN->data can be cast to a BehaviorBoutFeatures*

typedef struct {
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

  int num_features;
} SVMFeatureParams;



typedef struct _BehaviorBout {
  int start_frame, end_frame;
  int behavior;
  double bout_score, transition_score, loss_fn, loss_fp;
} BehaviorBout;

typedef struct _BehaviorBoutSequence {
  struct _BehaviorGroups *behaviors;
  struct _BehaviorBoutFeatures *features;
  int *num_bouts;
  BehaviorBout **bouts;
  double score, loss, slack;
  double *scores, *losses;
} BehaviorBoutSequence;


typedef struct _BehaviorBoutFeatures {
  BehaviorBoutSequence *partial_label;
  unsigned char *memory_buffer;

  /* Precomputed feature caches */
  double **features;             // A num_base_features X T array of all features
  double *frame_times;           // A num_base_features array of all frame times
  double **smoothed_features;    // A num_base_features X T array of all features smoothed around some temporal window
  double **integral_features;    // A num_base_features X T array encoding integral (sum) features
  double **integral_sqr_features;// A num_base_features X T array encoding integral (sum) of squared features
  double *max_feature_responses; // A num_base_features array encoding the global max response of each feature
  double *min_feature_responses; // A num_base_features array encoding the global min response of each feature
  double *ave_feature_responses; // A num_base_features array encoding the global ave response of each feature
  double ***integral_histogram_features;  // A num_base_features X num_thresholds X T encoding integral features for each histogram bin
  int **histogram_bins;          // A num_base_features X T array encoding the histogram index of each frame feature

  /* Feature caches that are updated during dynamic programming */
  double *bout_max_feature_responses; // A num_base_features array encoding the max response of each feature in the current bout
  double *bout_min_feature_responses; // A num_base_features array encoding the min response of each feature in the current bout
  int bout_start, bout_end;

  void *data;
  char fname[400];
  int num_frames;

  SVECTOR *fvec;
} BehaviorBoutFeatures;


#define MAX_TEMPORAL_LEVELS 8
#define MAX_HARMONIC_LEVELS 8
#define MAX_BASE_FEATURES 1000

class SVMBehaviorSequence : public SVMStructMethod {
 protected:
  struct _BehaviorGroups *behaviors;
  int num_classes[MAX_BEHAVIOR_GROUPS];
  int num_features, num_base_features;
  int feature_diff_frames;
  int sizePsi;
  int behavior;   // if not -1, only run the model on one behavior group
  double *false_negative_cost[MAX_BEHAVIOR_GROUPS], *false_positive_cost[MAX_BEHAVIOR_GROUPS];
  int ***class_training_transitions, **class_training_transitions_count, **class_training_count;
  double time_approximation;
  SVMFeatureParams feature_params[MAX_BASE_FEATURES];
  double *features_mu, *features_gamma;
  double **histogram_thresholds, **min_thresholds, **max_thresholds;
  int min_bout_duration;
  char **feature_names;
  
 public:
  SVMBehaviorSequence(int num_feat, struct _BehaviorGroups *behaviors, int beh, SVMFeatureParams *sparams = NULL);
  ~SVMBehaviorSequence();

  int train (int argc, const char* argv[], STRUCT_LEARN_PARM *struct_parm, STRUCTMODEL *structmodel);
  int test (int argc, const char* argv[], STRUCT_LEARN_PARM *struct_parm, STRUCTMODEL *structmodel);

  SVMFeatureParams DefaultParams();
  void compute_bout_feature_caches(BehaviorBoutFeatures *blob_features);
  int compute_feature_space_size();
  void update_bout_feature_caches(BehaviorBoutFeatures *b, int t_start, int t_end, int c);

  void print_features(FILE *fout, double *feat);
  void print_features(const char *fname, EXAMPLE *ex, int n, bool normalized);
  void print_weights(FILE *fout, double *w);
  void print_weights(const char *fname, double *w);
  void set_feature_name(int feature_ind, int base_feature_ind, const char *name);

  void on_finished_iteration(CONSTSET c, STRUCTMODEL *sm, 
						STRUCT_LEARN_PARM *sparm, int iter_num) ;

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
  
  double      *psi_bout(BehaviorBoutFeatures *b, int t_start, int t_end, int c, double *feat, bool normalize=true, bool fast_update=false);
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

 
  virtual const char *get_base_feature_name(int ind) = 0;
  virtual void load_from_bout_sequence(BehaviorBoutSequence *y, void *b) = 0;
  virtual BehaviorBoutSequence *create_behavior_bout_sequence(void *b, BehaviorGroups *behaviors, bool build_partial_label) = 0;
  virtual void load_behavior_bout_features(void *b, BehaviorBoutFeatures *feature_cache) = 0; 
  virtual void *load_training_example(const char *fname, BehaviorGroups *behaviors) = 0;
  virtual void free_data(void *d) = 0;
  virtual int num_frames(void *d) = 0;
  virtual char **load_examples(const char *fname, int *num) = 0;
};

void free_behavior_bout_sequence(BehaviorBoutSequence *b, int num);

#endif
