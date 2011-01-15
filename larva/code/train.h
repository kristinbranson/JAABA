#ifndef __TRAIN_H
#define __TRAIN_H

#include "blob.h"

#define DEFAULT_DATA_DIR "../data"

extern const char *DATA_DIR;

typedef enum {
  CLASSIFIER_NONE,
  CLASSIFIER_DECISION_TREE,
  CLASSIFIER_RANDOM_FOREST,
  CLASSIFIER_BOOSTING,
  CLASSIFIER_BOOSTED_TREES,
  CLASSIFIER_NEAREST_NEIGHBOR,
  CLASSIFIER_SVM,
  CLASSIFIER_NORMAL_BAYES,
  CLASSIFIER_SVM_LIGHT,
  CLASSIFIER_SVM_STRUCT_BEHAVIOR
} ClassifierMethod;

void predict_behaviors(Blob *b, BehaviorGroups *behaviors, int *preds);
void predict_all_behaviors(BlobSequence *b, BehaviorGroups *behaviors);

void train_behavior_detector(const char *dir, char **train_list, int num_files, BehaviorGroups *groups,
			     struct _FitParams *params);

extern const char *g_classifer_extensions[400];

char **load_train_list(const char *fname, int *num);
int save_train_list(const char *fname, char **tlist, int num);
char *get_classifier_name(char *classifier_name, const char *classifier_dir, const char *b_name, int method);
CvStatModel *allocate_classifier(BehaviorGroups *behaviors, int classifier_method);
void deallocate_classifier(CvStatModel *m, int classifier_method);
void free_train_list(char **list, int num);


#endif
