#include "train.h"
#include "fit_model.h"

#include <cv.h>
#include <highgui.h>
#include <ml.h>  

#define SVM_LIGHT_C 1000

#include "svm/svm_struct/svm_struct_api_multiclass.h"
#include "svm/svm_struct/svm_struct_api_blob_behavior_sequence.h"

const char *DATA_DIR = StringCopy(DEFAULT_DATA_DIR);

class CvSVMLight /*: public CvStatModel*/ {
  STRUCTMODEL model;
  STRUCT_LEARN_PARM sparm;
  SVMMulticlass m;

public:
  CvSVMLight() {
    memset(&model, 0, sizeof(STRUCTMODEL));
    memset(&sparm, 0, sizeof(STRUCT_LEARN_PARM));
  }
  ~CvSVMLight() {
    m.free_struct_model(model);
  }

  void save(const char *filename, const char* name=0 ) {
    m.write_struct_model(filename, &model, &sparm);
  }
  void load(const char* filename, const char* name=0 ) {
    model=m.read_struct_model(filename,&sparm);
  }
  void write( CvFileStorage* storage, const char* name ) { assert(0); }
  void read( CvFileStorage* storage, CvFileNode* node ) { assert(0); }

  float predict(CvMat *mat, float *scores=NULL) {
    int pred = 0;
    m.svm_struct_classify_from_array((float*)mat->data.fl, &pred, scores, 1, &model, &sparm);
    return (float)pred;
  }

  bool train( const CvMat* train_data, const CvMat* responses, const CvMat *var_idx, const CvMat *sample_idx, double c = 1) {
    double mi, num_classes;
    char c_str[400]; sprintf(c_str, "%lf", c);
    const char *ptr = c_str;
    const char *args[3] = {"train", "-c", ptr};
    int num_args = 3;

    cvMinMaxLoc(responses, &mi, &num_classes);
    m.svm_struct_learn_from_array(&sparm, &model, (float*)train_data->data.fl, (int*)responses->data.i, sample_idx ? (unsigned char*)sample_idx->data.ptr : NULL,
				train_data->rows, train_data->cols, (int)num_classes, num_args, args);

    return true;
  }
};


class CvSVMStructBehaviorSequence : public CvStatModel {
  STRUCTMODEL model;
  STRUCT_LEARN_PARM sparm;
  SVMBlobBehaviorSequence *svm_struct;
  int group;
  BehaviorGroups *behaviors;

public:
  CvSVMStructBehaviorSequence(FitParams *params, BehaviorGroups *behaviors, int group = -1) {
    FitParams p = params ? *params : default_parameters();
    svm_struct = new SVMBlobBehaviorSequence(&p, behaviors, -1);
    memset(&model, 0, sizeof(STRUCTMODEL));
    memset(&sparm, 0, sizeof(STRUCT_LEARN_PARM));
    this->group = group;
    this->behaviors = behaviors;
  }
  ~CvSVMStructBehaviorSequence() {
    delete svm_struct;
  }

  void save(const char *filename, const char* name=0 ) {
    svm_struct->write_struct_model(filename, &model, &sparm);
  }
  void load(const char* filename, const char* name=0 ) {
    model=svm_struct->read_struct_model(filename,&sparm);
  }
  void write( CvFileStorage* storage, const char* name ) { assert(0); }
  void read( CvFileStorage* storage, CvFileNode* node ) { assert(0); }

  float predict(BlobSequence *blobs) {
    LABEL yy;
    BehaviorBoutSequence *y;
    float retval = 0;
    SPATTERN x;
    BehaviorBoutFeatures *blob_features;

    x.data = blob_features = svm_struct->create_behavior_bout_feature_cache(blobs);
    blob_features->partial_label = svm_struct->create_behavior_bout_sequence(blobs, behaviors, true);
    yy = svm_struct->classify_struct_example(x, &model, &sparm);
    y = (BehaviorBoutSequence*)yy.data;
    retval = (float)y->score;
    svm_struct->load_from_bout_sequence(y, blobs);
    svm_struct->free_label(yy);
    blob_features->data = NULL;
    svm_struct->free_pattern(x);
    return retval;
  }

  bool train(const char *trainfile, const char *modelfile=NULL, const char *constraints_file=NULL, double c = 1000) {
    char c_str[400]; sprintf(c_str, "%lf", c);
    char eps_str[400]; sprintf(eps_str, "%lf", c/500);
    char beh_str[400]; sprintf(beh_str, "%d", group);
    const char *ptr = c_str;
    //const char *ptr_beh = beh_str;
    const char *ptr_eps = eps_str;
    const char *argv[20] = { "./svm_struct_train", "-c", ptr, "-e", ptr_eps};
    int argc = 5;
    //if(constraints_file) { argv[argc++] = "-i"; argv[argc++] = constraints_file; } 
    if(constraints_file) { argv[argc++] = "-j"; argv[argc++] = constraints_file; } 
    argv[argc++] = trainfile;
    if(modelfile)argv[argc++] = modelfile;

    svm_struct->train(argc, argv, &sparm, &model);

    return true;
  }

  bool test(const char *testfile, const char *resultsfile=NULL, const char *modelfile=NULL, 
	    const char *constraints_file=NULL, double c = 1000) {
    const char *argv[20] = { "./svm_struct_test"};
    int argc = 1;

    argv[argc++] = testfile;
    if(resultsfile) { argv[argc++] = "-o"; argv[argc++] = resultsfile; } 
    svm_struct->test(argc, argv, &sparm, &model);
  }
};




CvStatModel *allocate_classifier(BehaviorGroups *behaviors, int classifier_method) {
  switch(classifier_method) {
  case CLASSIFIER_SVM:
    return new CvSVM;
  case CLASSIFIER_SVM_LIGHT:
    return (CvStatModel*)new CvSVMLight;
  case CLASSIFIER_SVM_STRUCT_BEHAVIOR:
    return (CvStatModel*)new CvSVMStructBehaviorSequence(NULL, behaviors, -1);
  case CLASSIFIER_NORMAL_BAYES:
    return new CvNormalBayesClassifier;
  case CLASSIFIER_DECISION_TREE:
    return new CvDTree;
  default:
    fprintf(stderr, "Error classifier %d is unknown\n", classifier_method);
    exit(0);
  }
}

void deallocate_classifier(CvStatModel *m, int classifier_method) {
  if(m) delete m;
}

char *get_classifier_name(char *classifier_name, const char *classifier_dir, const char *b_name, int method) {
  sprintf(classifier_name, "%s/%s.%s", classifier_dir, method == CLASSIFIER_SVM_STRUCT_BEHAVIOR ? "All" : b_name, 
	  g_classifer_extensions[method]);
  return classifier_name;
}

const char *g_classifer_extensions[400] = {
  "none", "dec.tree", "rand.forest", "boost", "tree.boost", "near.neighbor", "svm", "nbayes", "svmlight", "svmbehavior"
};

void train_behavior_detector(const char *dir, char **train_list, int num_files, BehaviorGroups *behaviors, FitParams *params) {
  int i, j, k, l, num_alloc = 0, num_train = 0;
  BehaviorGroup *beh;
  BlobSequence *blobs;
  float *train_features = NULL, *curr;
  double *g_feat;
  CvMat train_mat, labels_mat, mask_mat;
  Blob *b;
  int **labels[MAX_BEHAVIOR_GROUPS];
  int *class_labels[MAX_BEHAVIOR_GROUPS];
  int *num_pos_examples[MAX_BEHAVIOR_GROUPS];
  unsigned char *mask[MAX_BEHAVIOR_GROUPS];
  char fname[400], constraints_file[400];
  ClassifierMethod method = params->classifier_method;
  
  if(method == CLASSIFIER_SVM_STRUCT_BEHAVIOR) {
    char train_file[1000];
    sprintf(fname, "%s/All.%s", dir, g_classifer_extensions[method]);
    sprintf(train_file, "train.txt.tmp");
    save_train_list(train_file, train_list, num_files); 
    sprintf(constraints_file, "%s/All.%s.constraints", dir, g_classifer_extensions[method]);
    CvSVMStructBehaviorSequence *bs = new CvSVMStructBehaviorSequence(params, behaviors, -1);
    bs->train(train_file, fname, constraints_file);
    return;
  }


  for(k = 0; k < behaviors->num; k++) {
    class_labels[k] = NULL;
    mask[k] = NULL;
    labels[k] = (int**)malloc(sizeof(int*)*behaviors->behaviors[k].num_values);
    num_pos_examples[k] = (int*)malloc(sizeof(int*)*behaviors->behaviors[k].num_values);
    memset(labels[k], 0, sizeof(int*)*behaviors->behaviors[k].num_values);
    memset(num_pos_examples[k], 0, sizeof(int)*behaviors->behaviors[k].num_values);
  }

  // Load in all manually verified blob frames as training examples
  for(j = 0; j < num_files; j++) {
    blobs = load_blob_sequence(train_list[j], behaviors);
    for(i = 0; i < blobs->num_frames; i++) 
      if(!blobs->frames[i].features)
	blobs->frames[i].features = (double*)malloc(2*sizeof(double)*((2*F_NUM_FEATURES)*(params->num_spine_lines+1)+F_NUM_GLOBAL_FEATURES));
    //compute_all_features(blobs, params, 0);
    compute_deterministic_spine_attributes_blob_sequence(blobs, params, 0);
    compute_all_global_features(blobs, params, 0);

    for(i = 0; i < blobs->num_frames; i++) {
      b = &blobs->frames[i];
      g_feat = NORMALIZED_GLOBAL_FEATURES(b);
      if(b->is_good && (b->is_manual & 3) == 3 && b->num_model_pts == params->num_spine_lines+1) {
	if(num_train >= num_alloc) {
	  num_alloc += 128;
	  train_features = (float*)realloc(train_features, sizeof(float)*F_NUM_GOOD_GLOBAL_FEATURES*num_alloc);
	  for(k = 0; k < behaviors->num; k++) {
	    for(l = 0; l < behaviors->behaviors[k].num_values; l++) {
	      labels[k][l] = (int*)realloc(labels[k][l], sizeof(int)*num_alloc);
	    }
	    class_labels[k] = (int*)realloc(class_labels[k], sizeof(int)*num_alloc);
	    mask[k] = (unsigned char*)realloc(mask[k], sizeof(unsigned char)*num_alloc*F_NUM_GOOD_GLOBAL_FEATURES);
	  }
	}
	curr = train_features + num_train*F_NUM_GOOD_GLOBAL_FEATURES;
	for(k = 0; k < F_NUM_GOOD_GLOBAL_FEATURES; k++)
	  curr[k] = (float)g_feat[k];
		  
	for(k = 0; k < behaviors->num; k++) {
	  for(l = 0; l < behaviors->behaviors[k].num_values; l++)
	    labels[k][l][num_train] = l == b->behaviors[k] ? 1 : -1;
	  
	  class_labels[k][num_train] = b->behaviors[k];
	  mask[k][num_train] = b->behaviors[k] != 0;
	  num_pos_examples[k][b->behaviors[k]]++;
	}
	num_train++;
      }
    }
  }


  cvInitMatHeader(&train_mat, num_train, F_NUM_GOOD_GLOBAL_FEATURES, CV_32FC1, train_features); 
  for(k = 0; k < behaviors->num; k++) {
    beh = &behaviors->behaviors[k];
    beh->classifier_method = method;
    beh->is_multiclass = params->is_multiclass;
    if(params->is_multiclass) {
      cvInitMatHeader(&labels_mat, 1, num_train, CV_32SC1, class_labels[k]);
      cvInitMatHeader(&mask_mat, 1, num_train, CV_8UC1, mask[k]);
      switch(method) {
      case CLASSIFIER_DECISION_TREE:
	beh->classifier = new CvDTree;
	((CvDTree*)beh->classifier)->train(&train_mat, CV_ROW_SAMPLE, &labels_mat, NULL, NULL, &labels_mat, NULL, 
					   CvDTreeParams(4, 5, 0, false, beh->num_values, 10, false, true, NULL));
	break;
      case CLASSIFIER_SVM:
	beh->classifier = new CvSVM;
	((CvSVM*)beh->classifier)->train(&train_mat, &labels_mat, NULL, &mask_mat, 
					 CvSVMParams(CvSVM::C_SVC, CvSVM::RBF, 0, 0.1, 0, 10.0, 0, 0, NULL, 
						     cvTermCriteria(CV_TERMCRIT_EPS, 10000000, FLT_EPSILON)));
	break;
      case CLASSIFIER_NORMAL_BAYES:
	beh->classifier = new CvNormalBayesClassifier;
	((CvNormalBayesClassifier*)beh->classifier)->train(&train_mat, &labels_mat, NULL, &mask_mat, false);
	break;
	
      case CLASSIFIER_SVM_LIGHT:
	beh->classifier = (CvStatModel*)new CvSVMLight;
	((CvSVMLight*)beh->classifier)->train(&train_mat, &labels_mat, NULL, &mask_mat, SVM_LIGHT_C);
	
	break;
      default:
	fprintf(stderr, "Classifier %s not supported yet for multi-class\n", g_classifer_extensions[method]);
	exit(0);
      }
      sprintf(fname, "%s/%s.%s", dir, beh->name, g_classifer_extensions[method]);
      beh->classifier->save(fname);
    } else {
      // train a bunch of 1-vs-all classifiers
      for(l = 0; l < beh->num_values; l++) {
	beh->values[l].classifier_method = 0;
	if(num_pos_examples[k][l] > 5 && l) {
	  beh->values[l].classifier_method = method;
	  cvInitMatHeader(&labels_mat, 1, num_train, CV_32SC1, labels[k][l]);
	  switch(method) {
	  case CLASSIFIER_SVM:
	    beh->values[l].classifier = new CvSVM;
	    ((CvSVM*)beh->values[l].classifier)->train(&train_mat, &labels_mat, NULL, NULL, 
				CvSVMParams(CvSVM::ONE_CLASS, CvSVM::LINEAR, 1.0, 0.0001, 1.0, 10.0, 0.5, 0.1, NULL, 
					    cvTermCriteria(CV_TERMCRIT_EPS, 10000000, FLT_EPSILON)));
	    break;
	  default:
	    fprintf(stderr, "Classifier %s not supported yet for 1-vs-all\n", g_classifer_extensions[method]);
	    exit(0);
	  }
	  sprintf(fname, "%s/%s_%s.%s", dir, beh->name, beh->values[l].name, g_classifer_extensions[method]);
	  beh->values[l].classifier->save(fname);
	}
      }
    }
  }

  for(k = 0; k < behaviors->num; k++) {
    for(l = 0; l < behaviors->behaviors[k].num_values; l++) 
      if(labels[k][l]) free(labels[k][l]);
    free(labels[k]);
    free(class_labels[k]);
  }
  free(train_features);
}

void predict_behaviors(Blob *b, BehaviorGroups *behaviors, int *preds) {
  int i, j;
  double margins[10000];
  CvMat mat;
  double mi, ma;
  float feat[F_NUM_GOOD_GLOBAL_FEATURES];
  double *g_feat = NORMALIZED_GLOBAL_FEATURES(b);
  float scores[10000];

  if(!preds)
    preds = b->behaviors;

  assert(b->features);
  for(i = 0; i < F_NUM_GOOD_GLOBAL_FEATURES; i++)
    feat[i] = (float)g_feat[i];
  cvInitMatHeader(&mat, 1, F_NUM_GOOD_GLOBAL_FEATURES, CV_32FC1, feat); 

  for(i = 0; i < behaviors->num; i++) {    
    if(behaviors->behaviors[i].is_multiclass) {
      if(!behaviors->behaviors[i].classifier)
	continue;

      switch(behaviors->behaviors[i].classifier_method) {
      case CLASSIFIER_DECISION_TREE:
	preds[i] = (int)((CvDTree*)behaviors->behaviors[i].classifier)->predict(&mat)->value;
	break;
      case CLASSIFIER_SVM:
	preds[i] = (int)((CvSVM*)behaviors->behaviors[i].classifier)->predict(&mat);
	break;
      case CLASSIFIER_SVM_LIGHT:
	preds[i] = (int)((CvSVMLight*)behaviors->behaviors[i].classifier)->predict(&mat, scores);
	break;
      case CLASSIFIER_NORMAL_BAYES:
	preds[i] = (int)((CvNormalBayesClassifier*)behaviors->behaviors[i].classifier)->predict(&mat);
	break;
      default:
	fprintf(stderr, "Classifier %s not supported yet for multi-class\n", 
		g_classifer_extensions[behaviors->behaviors[i].classifier_method]);
	exit(0);
      }
    } else {
      mi = INFINITY;
      ma = -INFINITY;
      for(j = 0; j < behaviors->behaviors[i].num_values; j++) {
	if(behaviors->behaviors[i].values[j].classifier) {
	  if(!behaviors->behaviors[i].values[j].classifier)
	    continue;
	  switch(behaviors->behaviors[i].classifier_method) {
	  case CLASSIFIER_SVM:
	    margins[j] = ((CvSVM*)behaviors->behaviors[i].values[j].classifier)->predict(&mat);
	    break;
	  default:
	    fprintf(stderr, "Classifier %s not supported yet for 1-vs-all\n", 
		    g_classifer_extensions[behaviors->behaviors[i].values[j].classifier_method]);
	    exit(0);
	  }
	  if(margins[j] < mi) mi = margins[j];
	  if(margins[j] > ma) { ma = margins[j]; preds[i] = j; }
	}
      }
    }
  }
}

void predict_all_behaviors(BlobSequence *b, BehaviorGroups *behaviors) {
  int i;

  for(i = 0; i < b->num_frames; i++) {
    if(!b->frames[i].model_pts) {
      fprintf(stderr, "ERROR: can't predict behaviors without detecting the skeleton first!\n");
      return;
    }
  }

  if(behaviors->classifier_method == CLASSIFIER_SVM_STRUCT_BEHAVIOR && behaviors->classifier) {
    ((CvSVMStructBehaviorSequence*)behaviors->classifier)->predict(b);
  }

  for(i = 0; i < b->num_frames; i++)
    predict_behaviors(&b->frames[i], behaviors, NULL);
}

char **load_train_list(const char *fname, int *num) {
  char **retval = NULL, line[1000];
  FILE *fin = fopen(fname, "r");
  *num = 0;

  if(!fin) return NULL;
  
  while(fgets(line, 999, fin)) {
    chomp(line);
    retval = (char**)realloc(retval, sizeof(char*)*(*num+1));
    retval[*num] = (char*)malloc(strlen(line)+1);
    strcpy(retval[*num], line);
    *num = *num+1;
  }
  fclose(fin);

  return retval;
}

void free_train_list(char **list, int num) {
  for(int i = 0; i < num; i++) {
    free(list[i]);
  }
  free(list);
}

int save_train_list(const char *fname, char **tlist, int num) {
  FILE *fout = fopen(fname, "w");
  if(!fout) return 0;
  
  for(int i = 0; i < num; i++)
    fprintf(fout, "%s\n", tlist[i]);
  
  fclose(fout);

  return 1;
}
