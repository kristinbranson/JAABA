#ifndef  __SVM_STRUCT_API_FLY_BEHAVIOR_SEQUENCE__
#define __SVM_STRUCT_API_FLY_BEHAVIOR_SEQUENCE__

#include "svm_struct_api_behavior_sequence.h"

typedef struct {
  BehaviorGroups *behaviors;

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

  BehaviorBoutSequence *bouts;
} FlyBehaviorBoutSequence;

typedef struct {
  char name[1001];
  char units_numer[1001];
  char units_denom[1001];
} FlyBehaviorFeatures;


class SVMFlyBehaviorSequence : public SVMBehaviorSequence {
FlyBehaviorFeatures feature_defs[MAX_BASE_FEATURES];

 public:
  SVMFlyBehaviorSequence(const char *fname, struct _BehaviorGroups *behaviors, int beh);

  bool read_features(FILE *fin,  FlyBehaviorFeatures *f, double *frames, int num_frames) ;

  int ReadFeatureParams(const char *fname, SVMFeatureParams *p);
  const char *get_base_feature_name(int ind);
  int num_frames(void *b);
  void load_from_bout_sequence(BehaviorBoutSequence *y, void *b);
  BehaviorBoutSequence *create_behavior_bout_sequence(void *b, BehaviorGroups *behaviors, bool build_partial_label);
  void load_behavior_bout_features(void *b, BehaviorBoutFeatures *feature_cache);
  void *load_training_example(const char *fname);
  void save_example(void *b, void *d, const char *fname);
  void free_data(void *d);
  char **load_examples(const char *fname, int *num);
  void save_examples(const char *fname, SAMPLE sample);
  char *getLabelName(void* d);
};

#endif
