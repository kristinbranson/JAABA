#ifndef  __SVM_STRUCT_API_FLY_BEHAVIOR_SEQUENCE__
#define __SVM_STRUCT_API_FLY_BEHAVIOR_SEQUENCE__

#include "svm_behavior_sequence.h"

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
                     "use_end_ave_diff_haar_features=%d"

#define FORMAT__BOUT_FEATURE_PARAMS_READ FORMAT__BOUT_FEATURE_PARAMS \
        "%*[^\n]\n"


class FlyBehaviorBoutFeatures : public BehaviorBoutFeatures {
 public:
  bool load(const char *fname, SVMBehaviorSequence *svm, BehaviorBoutSequence *y);
};

class FlyBehaviorBoutSequence : public BehaviorBoutSequence {
  double version;
  char labelname[1001];
  char moviename[1001];
  char matname[1001];
  char trxname[1001];
  int nflies;
  int fly_ids[100];
  int firstframe, lastframe;

  char sex;
  double fps;
  bool is_labeled;

public:
  FlyBehaviorBoutSequence(BehaviorBoutFeatures *x, SVMBehaviorSequence *svm);
  bool load(const char *fname);
  bool save(const char *fname);

  friend class FlyBehaviorBoutFeatures;
  friend class SVMFlyBehaviorSequence;
};

typedef struct {
  char name[1001];
  char units_numer[1001];
  char units_denom[1001];
} FlyBehaviorFeatures;


class SVMFlyBehaviorSequence : public SVMBehaviorSequence {
  FlyBehaviorFeatures feature_defs[MAX_BASE_FEATURES];

 public:
  SVMFlyBehaviorSequence(const char *fname, struct _BehaviorGroups *behaviors, int beh);
  StructuredLabel *NewStructuredLabel(StructuredData *x);
  StructuredData *NewStructuredData();


  bool read_features(FILE *fin,  FlyBehaviorFeatures *f, double *frames, int num_frames) ;

  FlyBehaviorFeatures *GetFeatureDefs() { return feature_defs; }

  int ReadFeatureParams(const char *fname, SVMFeatureParams *p);
  bool ReadFeatureParam(FILE *modelfl, SVMFeatureParams *p);
  const char *get_base_feature_name(int ind);
  void load_from_bout_sequence(BehaviorBoutSequence *y, void *b);
  BehaviorBoutSequence *create_behavior_bout_sequence(void *b, BehaviorGroups *behaviors, bool build_partial_label);
  void load_behavior_bout_features(void *b, BehaviorBoutFeatures *feature_cache);
  void *load_training_example(const char *fname);
  void save_example(void *b, void *d, const char *fname);
  void free_data(void *d);
  char **load_examples(const char *fname, int *num);
  void save_examples(const char *fname, StructuredDataset *dataset);
};

#endif
