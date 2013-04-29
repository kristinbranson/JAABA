#ifndef  __SVM_STRUCT_API_FLY_BEHAVIOR_SEQUENCE__
#define __SVM_STRUCT_API_FLY_BEHAVIOR_SEQUENCE__

#include "svm_behavior_sequence.h"

class FlyBehaviorBoutSequence;
class FlyBehaviorBoutFeatures;


class FlyBehaviorBoutFeatures : public BehaviorBoutFeatures {
  int fly_id;
  char sex;
  int nflies;
  int fly_ids[100];
  double version;
  char moviename[1001];
  char matname[1001];
  int firstframe, lastframe;

 public:
  FlyBehaviorBoutFeatures();
  bool load(const char *fname, SVMBehaviorSequence *svm);
  bool load(const Json::Value &x, StructuredSVM *s);
  Json::Value save(StructuredSVM *s);
  int GetFirstFrame() { return firstframe; }
  int GetLastFrame() { return lastframe; }

  friend class FlyBehaviorBoutSequence;
};

class FlyBehaviorBoutSequence : public BehaviorBoutSequence {
  double version;
  int fly_id;
  int firstframe, lastframe;
  bool is_labeled;
  char labelname[1001];
  char moviename[1001];
  char matname[1001];
  char trxname[1001];
  int nflies;
  int fly_ids[100];
  SVMBehaviorSequence *svm;

public:
  FlyBehaviorBoutSequence(BehaviorBoutFeatures *x, SVMBehaviorSequence *svm);
  bool load(const char *fname);
  bool save(const char *fname);
  bool load(const Json::Value &x, StructuredSVM *s);
  Json::Value save(StructuredSVM *s);
  void AddBout(int behavior, int start_frame, int end_frame) {
    bouts = (BehaviorBout*)realloc(bouts, (num_bouts+1)*sizeof(BehaviorBout));
    bouts[num_bouts].behavior = behavior;
    bouts[num_bouts].start_frame = start_frame;
    bouts[num_bouts].end_frame = end_frame;
    num_bouts++;
  }

  friend class FlyBehaviorBoutFeatures;
  friend class SVMFlyBehaviorSequence;
};


class SVMFlyBehaviorSequence : public SVMBehaviorSequence {
 public:
  SVMFlyBehaviorSequence();
  StructuredLabel *NewStructuredLabel(StructuredData *x);
  StructuredData *NewStructuredData();

  char **load_examples(const char *fname, int *num);
  void save_examples(const char *fname, StructuredDataset *dataset);
};

#endif
