#include <stdio.h>
#include <string.h>
#include "svm_behavior_sequence.h"

//#define DEBUG 0
//#define ALLOW_SAME_TRANSITIONS 0
#if DEBUG > 0
char *g_currFile; // CSC 20110420: hack to pass current filename for debug purposes
#endif

#define KEEP_FEATURES_IN_MEMORY 0


#define CAP_BOUT_FEATURES 4.0  // cap bout-level features to be no more than 4 standard deviations

// Assume no score (other than transition costs) for the "none" behavior.  Effectively, removes all
// bout features, unary score, and duration score for the none behavior
// This allows an exhaustive search for all bout durations for the "none" behavior
#define NONE_CLASS_HAS_NO_SCORE 1

#define INST_VERSION "0.0.0"

#if USE_DURATION_COST > 0
#define getPsiSize(p_features, p_classes) ((p_features) + (p_classes) + 1 + 1) * (p_classes)
#else
#define getPsiSize(p_features, p_classes) ((p_features) + (p_classes) + 1) * (p_classes)
#endif



/*
* Routines to fit an optimal scoring sequence of behavior bouts given a tracking sequence
* and to learn the optimal parameters for that scoring function from training data
*
* The most important functions to look at are:
*
* inference_via_dynamic_programming(): find the optimal behavior segmentation at test time 
*   given a known score function, or find the most violated constraint (highest scoring 
*   segmentation with a loss term incorporated) at train time
*
* psi_bout(): compute bout-level features (such as bout average, sum, min, max, etc.) from 
*   frame-level features.  Uses cache structures created via compute_bout_feature_caches()
*
* loss2(): evaluate the loss of a proposed bout ybar with respect to ground truth bout y
*
* The other functions are mostly used to convert between different data formats or initialize
* parameters
*
*/

SVMBehaviorSequence::SVMBehaviorSequence() : StructuredSVM() {
  restrict_behavior_features = NULL;
#if USE_DURATION_COST > 0
  min_frame_duration = max_frame_duration = NULL;
#endif
  false_negative_cost = false_positive_cost = NULL;
  Init(NULL, NULL, 0, NULL, 0);
}

void SVMBehaviorSequence::Init(Behaviors *behaviors, FrameFeature *frame_features, int num_frame_features, 
			       BoutFeature *bout_features, int num_bout_features) {
  int i;

  eps = .015;
  C = 10.0;
  method = SPO_DUAL_MULTI_SAMPLE_UPDATE_WITH_CACHE;// SPO_DUAL_UPDATE_WITH_CACHE;
  debugLevel = 3;

  this->bout_features = bout_features;
  this->num_bout_features = num_bout_features;

  this->frame_features = frame_features;
  this->num_frame_features = num_frame_features;

  max_inference_learning_frames = -1;  // don't disregard any frames at train time to speedup inference

  // Speed up inference, only searching over bout durations of size time_approximation^k, for some integer k.
  // However, do consider merging bouts of the same class, such that we can still predict bouts of any length
  // time_approximation = 0;   // disable approximate inference
  time_approximation = 1.2;
  //time_approximation = -1; 
  search_all_bout_durations_up_to = 50; // Search all bout durations from 1 to 50.  Can be combined with time_approximation

  runMultiThreaded = 1; 

  importance_sample_interval_size = 200;  // break the sequence into sub-sections of size 200

  numCacheUpdatesPerIteration = 50;
  maxCachedSamplesPerExample = 500;
  numMultiSampleIterations = 10;
  keepAllEvictedLabels = true;
  updateFromCacheThread = true;
  
  strcpy(debugdir, "");
  debug_predictions = debug_weights = debug_features = debug_model = false;

  this->behaviors = behaviors;

  class_training_transitions = NULL;
  class_training_transitions_count = class_training_count = NULL;
  min_bout_duration = 1;

  sizePsi = 0;
  if(!behaviors)
    return;

  false_negative_cost = (double*)realloc(false_negative_cost, 2*sizeof(double)*behaviors->num_values);
  false_positive_cost = false_negative_cost + behaviors->num_values;
#if USE_DURATION_COST > 0
  min_frame_duration  = (double*)realloc(min_frame_duration, 2*sizeof(double)*behaviors->num_values);
  max_frame_duration  = min_frame_duration+behaviors->num_values;
#endif

  // It is assumed here that the label 0 is the unknown label.  The user could define custom
  // class-specific values for these constants .
  // In the future, we could implement this with a per class-pair confusion cost
#if SCALED_LOSS > 0 
  false_negative_cost[0] = 0;  // 0 when beh 0 represents "other"
  false_positive_cost[0] = 0;  // 0 when beh 0 represents "other"
  for(i = 1; i < behaviors->num_values; i++) { 
    false_negative_cost[i] = 10; // 100 for eric's videos
    false_positive_cost[i] = 10; // 100 for eric's videos
  }
#else
  false_negative_cost[0] = 0;
  false_positive_cost[0] = 0;  //
  for(i = 1; i < behaviors->num_values; i++) {  
    false_negative_cost[i] = 1; 
    false_positive_cost[i] = 1;
  }
#endif
    
#if USE_DURATION_COST > 0
  for(i = 0; i < behaviors->num_values; i++) {
    min_frame_duration[i] = INFINITY;
    max_frame_duration[i] = -INFINITY;
  }
#endif
  sizePsi = getPsiSize(this->num_bout_features, behaviors->num_values);

  restrict_behavior_features = (bool**)realloc(restrict_behavior_features, behaviors->num_values*sizeof(bool*));
  memset(restrict_behavior_features, 0, sizeof(bool*)*behaviors->num_values);

  if(sizePsi) 
    fprintf(stderr, "Psi is %d-dimensional, over %d frame features, %d bout features\n", sizePsi, num_frame_features, num_bout_features);
}

SVMBehaviorSequence::~SVMBehaviorSequence() {
  if(class_training_transitions) free(class_training_transitions);
  if(class_training_transitions_count) free(class_training_transitions_count);
  if(class_training_count) free(class_training_count);
  
  FreeFrameFeatureParams();
  FreeBoutFeatureParams();
  FreeBehaviorDefinitions();
}





int cmp_double(const void * a, const void * b) { double d = *(double*)a - *(double*)b;  return (d < 0 ? -1 : (d > 0 ? 1 : 0)); }

void SVMBehaviorSequence::saveBoutFeatures(StructuredDataset *dataset, const char *filename, bool sphered, bool addRandBouts, int window_size) {
  // EYRUN: Store normalized bout features of training set in a txt file (want to verify that bout features are good enough to separate the classes)
  // Count the total number of bouts in the training set
  int num_bouts = 0;
  //int num_rand_factor = 5;
  //double window_size = 4; // 8, 16, 32, 64, 128
  StructuredExample **ex = dataset->examples;
  BehaviorBoutSequence *y;
  BehaviorBoutFeatures *x;
  double *tmp_features = (double*)malloc(sizeof(double)*num_bout_features);

  int num_bouts_extra = 0;
  for(int n = 0; n < dataset->num_examples; n++) {
    num_bouts += ((BehaviorBoutSequence*)ex[n]->y)->num_bouts;
    num_bouts_extra += ((BehaviorBoutFeatures*)ex[n]->x)->num_frames - window_size;
  }

  if(addRandBouts)
    num_bouts = num_bouts + num_bouts_extra; //num_bouts = num_bouts * (num_rand_factor+1); //

  FILE *featureFile;//, *featureFileUnsphered;
  featureFile = fopen(filename,"w");
  fprintf(featureFile, "num_features = %d\nnum_bouts = %d", num_bout_features, num_bouts);
  for(int n = 0; n < dataset->num_examples; n++) {
    y = ((BehaviorBoutSequence*)ex[n]->y);
    x = ((BehaviorBoutFeatures*)ex[n]->x);
    int *frame_labels = new int[x->num_frames];
    if(!x->memory_buffer)
      x->ComputeCaches(this);
    
    for(int j = 0; j < y->num_bouts; j++) {
      psi_bout(x, y->bouts[j].start_frame, y->bouts[j].end_frame, y->bouts[j].behavior, tmp_features, sphered, false);
      fprintf(featureFile, "\n%d %d %d %d ", n, y->bouts[j].start_frame, y->bouts[j].end_frame, y->bouts[j].behavior);
      for(int i = 0; i < num_bout_features; i++) {
	fprintf(featureFile, "%f ", tmp_features[i]);
      }
      for(int i = y->bouts[j].start_frame; i < y->bouts[j].end_frame; i++) 
	frame_labels[i] = y->bouts[j].behavior;
    }

    if(addRandBouts) {
      // Generate random intervals and calculate the bout-features for those
      //int num_rand = y->num_bouts[0] * num_rand_factor;
      //int num_rand = x->num_frames - window_size;
      for(int r = window_size/2; r<x->num_frames-window_size/2; r++) {
	int f_start, f_end;
	// generate random start frame between 0 and x->num_frames - 5
	//	f_start = rand() % (x->num_frames - 5);
	// generate random duration between 4 and 50 (the number of frames we go back)
	//	f_dur = rand() % 47 + 4;
	//	f_end = my_min(f_start + f_dur - 1, x->num_frames - 1);
	
	f_start = r-window_size/2;
	f_end = f_start+window_size;
	
	// calculate bout score for that interval
	psi_bout(x, f_start, f_end, -1, tmp_features, sphered, false);
	fprintf(featureFile, "\n%d %d %d %d ", n, f_start, f_end, frame_labels[r]);
	for(int i = 0; i < num_bout_features; i++) {
	  fprintf(featureFile, "%f ", tmp_features[i]);
	}
      }	
    }
    delete [] frame_labels;
#if KEEP_FEATURES_IN_MEMORY == 0
    x->Clear();
#endif	
  }
  fclose(featureFile);
  free(tmp_features);
}

/*
* Compute mean, variance and median statistics on the base features over the training set.
* the mean and variance are used to normalize each feature.  The median statistics are used to
* construct histogram and threshold constants, by placing some fraction of the training set 
* between each threshold
*/
void SVMBehaviorSequence::compute_feature_mean_variance_median_statistics(StructuredDataset *dataset) {
  int i, j, k, n, ind, num_bouts = 0, num_frames = 0;
  double w, target;
  BehaviorBoutSequence *y;
  double *tmp_features = (double*)malloc(sizeof(double)*num_bout_features);
  BehaviorBoutFeatures *x;
  StructuredExample **ex = dataset->examples;

  if(!trainset) trainset = dataset;
  ComputeClassTrainsitionCounts();
      
  // Count the total number of bouts in the training set
  for(n = 0; n < dataset->num_examples; n++) {
    num_frames += ((BehaviorBoutFeatures*)ex[n]->x)->num_frames;
    y = ((BehaviorBoutSequence*)ex[n]->y);
    num_bouts += y->num_bouts;
  }
  if(!num_bouts) {
    fprintf(stderr, "ERROR: no labeled bouts detected\n");
    return;
  }
  	

#if SCALED_LOSS > 0
  for(i=0; i < behaviors->num_values; i++){
    if(class_training_count[i]) {
      false_negative_cost[i] /= class_training_count[i];
      false_positive_cost[i] /= class_training_count[i];
    }
  }
#endif

  // Compute per-frame histogram bin thresholds
  fprintf(stderr, "Computing frame histogram thresholds...\n");
  double *frame_feat = (double*)malloc(num_frames*sizeof(double));
  for(i = 0; i < num_frame_features; i++) {
    // Extract all per-frame features into an array
    int curr_frame = 0;
    for(n = 0; n < dataset->num_examples; n++) {
      y = ((BehaviorBoutSequence*)ex[n]->y);
      x = ((BehaviorBoutFeatures*)ex[n]->x);
      int c = 0;
      for(j = 0; j < y->num_bouts; j++) {
	for(k = y->bouts[j].start_frame; k < y->bouts[j].end_frame; k++) {
	  assert(k >= 0);  // CSC 20110324
	  frame_feat[curr_frame++] = x->features[i][k];
	  c++;
	}
	if (y->bouts[j].end_frame == y->bouts[j].start_frame)
	  printf("Warning: bout has 0 frames!");
      }
    }

    // Compute thresholds for assigning histogram bins.  The thresholds are chosen such that an equal number
    // of training examples (training example frames) would lie in each histogram bin 
    qsort(frame_feat, curr_frame, sizeof(double), cmp_double);
    for(j = 0; j < frame_features[i].num_histogram_bins-1; j++) {
      target = (j+1)*curr_frame / (double)frame_features[i].num_histogram_bins;
      ind = (int)target;
      w = 1-(target-ind);
      while(j && ind < curr_frame && frame_feat[ind] == frame_features[i].histogram_thresholds[j-1]) ind++;
      frame_features[i].histogram_thresholds[j] = frame_feat[my_min(ind,curr_frame-1)]*w + 
	frame_feat[my_min(ind+1,curr_frame-1)]*(1-w);
      assert(!isnan(frame_features[i].histogram_thresholds[j]));
    }
  }
  free(frame_feat);

  // Extract all bout features for all training bouts into an array
  fprintf(stderr, "Computing bout feature thresholds...\n");
  double **bout_feat = (double**)malloc(num_bout_features*sizeof(double*));
  memset(bout_feat, 0, num_bout_features*sizeof(double*));
  i = 0;
  while(i < num_bout_features) {
    if(bout_features[i].num_thresholds > 1) {  // check if this bout feature is a histogram-type
      bout_feat[i] = (double*)malloc(num_bouts*sizeof(double));

      int curr_bout = 0;
      for(n = 0; n < dataset->num_examples; n++) {
	y = ((BehaviorBoutSequence*)ex[n]->y);
	x = ((BehaviorBoutFeatures*)ex[n]->x);
	// Now that frame histogram thresholds are set, we can compute bout-wise features for the training set
	if(!x->memory_buffer)
	  x->ComputeCaches(this);

	for(j = 0; j < y->num_bouts; j++) 
	  bout_feat[i][curr_bout++] = x->ComputeBoutFeature(&bout_features[i], y->bouts[j].start_frame, 
							    y->bouts[j].end_frame);
	
#if KEEP_FEATURES_IN_MEMORY == 0
	x->Clear();
#endif
      }
      assert(curr_bout == num_bouts);
      i += bout_features[i].num_thresholds;
    } else
      i++;
  }

  // Compute bout-feature histogram bin thresholds	
  i = 0;
  while(i < num_bout_features) {
    if(bout_features[i].num_thresholds > 1) {  // check if this bout feature is a histogram-type
      // Compute thresholds for assigning histogram bins.  The thresholds are chosen such that an equal number
      // of training examples (training example bouts) would lie in each histogram bin 
      qsort(bout_feat[i], num_bouts, sizeof(double), cmp_double);
      for(j = 0; j < bout_features[i].num_thresholds; j++) {
	target = (j+1)*num_bouts / (double)bout_features[i].num_thresholds;
	ind = (int)target;
	w = 1-(target-ind);
	while(j && ind < num_bouts && bout_feat[i][ind] == bout_features[i+j-1].thresh) ind++;
	bout_features[i+j].thresh = bout_feat[i][my_min(ind,num_bouts-1)]*w + bout_feat[i][my_min(ind+1,num_bouts-1)]*(1-w);
	assert(!isnan(bout_features[i+j].thresh));
      }
      for(j = 0; j < bout_features[i].num_thresholds; j++) 
	assert(bout_features[i+j].num_thresholds == bout_features[i+j].num_thresholds);
      free(bout_feat[i]);
      i += bout_features[i].num_thresholds;
    } else
      i++;
  }
  free(bout_feat);

  for(i = 0; i < num_bout_features; i++) {
    bout_features[i].mu = bout_features[i].gamma = 0;
  }

  fprintf(stderr, "Computing bout feature means...\n");
  for(n = 0; n < dataset->num_examples; n++) {
    y = ((BehaviorBoutSequence*)ex[n]->y);
    x = ((BehaviorBoutFeatures*)ex[n]->x);

    // Now that histogram thresholds are set, we can compute bout-wise features for the training set
    if(!x->memory_buffer)
      x->ComputeCaches(this);

    for(j = 0; j < y->num_bouts; j++) {
      psi_bout(x, y->bouts[j].start_frame, y->bouts[j].end_frame, y->bouts[j].behavior, tmp_features, false, false);
      for(i = 0; i < num_bout_features; i++) {
	// Now that we have bout features, we can compute the mean of each bout-wise feature over the training set
	assert(!isnan(tmp_features[i]));
	bout_features[i].mu  += tmp_features[i];
	assert(!isnan(bout_features[i].mu));
      }
    }
#if KEEP_FEATURES_IN_MEMORY == 0
    x->Clear();
#endif
  }
  for(i = 0; i < num_bout_features; i++) 
    bout_features[i].mu /= num_bouts;

  // Compute the (inverse of the) standard deviation bout-average feature response
  fprintf(stderr, "Computing bout feature standard deviations...\n");
  for(n = 0; n < dataset->num_examples; n++) {
    y = ((BehaviorBoutSequence*)ex[n]->y);
    x = ((BehaviorBoutFeatures*)ex[n]->x);
    if(!x->memory_buffer)
      x->ComputeCaches(this);
    
    for(j = 0; j < y->num_bouts; j++) {
      psi_bout(x, y->bouts[j].start_frame, y->bouts[j].end_frame, -1, tmp_features, false, false);
      for(i = 0; i < num_bout_features; i++) {
	bout_features[i].gamma += SQR(tmp_features[i]-bout_features[i].mu);
	assert(!isnan(bout_features[i].gamma));
      }
    }
#if KEEP_FEATURES_IN_MEMORY == 0
    x->Clear();
#endif
  }
  for(i = 0; i < num_bout_features; i++) {
    if(bout_features[i].gamma == 0)
      fprintf(stderr, "WARNING: feature %s seems to be the same for every example in the training set\n", bout_features[i].name);
    else
      bout_features[i].gamma = 1.0/sqrt(bout_features[i].gamma/num_bouts);
  }

  //saveBoutFeatures(dataset,"traindata_5.txt",true,false,0);
  //saveBoutFeatures(dataset,"traindata_5_unsphered.txt",false,false,0);
/*
	saveBoutFeatures(dataset, "train_bout_feat_4.txt",   false, true, 4);
	saveBoutFeatures(dataset, "train_bout_feat_8.txt",   false, true, 8);
	saveBoutFeatures(dataset, "train_bout_feat_16.txt",  false, true, 16);
	saveBoutFeatures(dataset, "train_bout_feat_32.txt",  false, true, 32);
	saveBoutFeatures(dataset, "train_bout_feat_64.txt",  false, true, 64);
	saveBoutFeatures(dataset, "train_bout_feat_128.txt", false, true, 128);
	//saveBoutFeatures(dataset, "train_bout_feat_unsphered.txt", false, true);
*/
  free(tmp_features);
}

/* 
* Compute bout-level features for a bout of class c on frames t_start to t_end.  This function
* supports several different kinds of expansions of frame-level features (e.g. tracker output or appearance
* features) into bout-level features (e.g. bout average, bout sum, histogram responses, temporal haar-features,
* harmonic motion, bout min and max, etc.)
* Currently, the features used are class-independent, and consist of the average frame features
* over the entire bout, and the change in feature response near the start and end of the bout.  
* If the parameter feat is non-NULL, these are stored into feat, otherwise the features
* are dynamically allocated.
*/
double *SVMBehaviorSequence::psi_bout(BehaviorBoutFeatures *b, int t_start, int t_end, int c, 
				      double *f, bool normalize, bool fast_update) {
  if(!f) 
    f = (double*)malloc(num_bout_features*sizeof(double));

  // Eyrun: temporary fix
  if(t_start == t_end)
    t_end += 1;

  int i = 0, j;
  double r;
  while(i < num_bout_features) {
    r = b->ComputeBoutFeature(&bout_features[i], t_start, t_end);
    if(!bout_features[i].num_thresholds)
      f[i++] = r;  // regular feature, raw feature response
    else {
      if(bout_features[i].num_thresholds==1) {
	// bout feature decision stump
	f[i] = r < bout_features[i].thresh;  
      } else {
	// bout feature histogram 
	for(j = i; j < i+bout_features[i].num_thresholds; j++) 
	  f[j] = 0;
	for(j = i; j < i+bout_features[i].num_thresholds; j++) {
	  if(r < bout_features[j].thresh) {
	    f[j] = 1;
	    break;
	  }
	}
      }
      i += bout_features[i].num_thresholds;
    } 
  }
  
  if(normalize)
    for(i = 0; i < num_bout_features; i++) 
      f[i] = (f[i]-bout_features[i].mu)*bout_features[i].gamma; 

  return f;
}


SparseVector SVMBehaviorSequence::Psi(StructuredData *x, StructuredLabel *yy) {
  BehaviorBoutSequence *y = (BehaviorBoutSequence*)yy, *y_partial = NULL;
  BehaviorBoutFeatures *b = (BehaviorBoutFeatures*)x;
  int i, j;
  double *tmp_features = (double*)malloc(sizeof(double)*(sizePsi+num_bout_features+10));
  double *all_features = tmp_features+num_bout_features+10;
  double *ptr = all_features, *class_features, *class_transitions, *class_counts; 
#if USE_DURATION_COST > 0
  double *class_durations;
  double duration, duration_diff;
#endif
  fill_unlabeled_gt_frames(y, y_partial);
  
  if(!b->memory_buffer)
    b->ComputeCaches(this);

  for(i = 0; i < sizePsi; i++)
    ptr[i] = 0;

  
  class_features = ptr; ptr += num_bout_features*behaviors->num_values;        // line 1
  class_transitions = ptr; ptr += behaviors->num_values*behaviors->num_values; // line 2
  class_counts = ptr; ptr += behaviors->num_values;   // unary feature
#if USE_DURATION_COST > 0
  class_durations = ptr; ptr += behaviors->num_values;   // duration feature
#endif

  for(i = 0; i < behaviors->num_values*num_bout_features; i++) 
    class_features[i] = 0;

  for(i = 0; i < behaviors->num_values*behaviors->num_values; i++) 
    class_transitions[i] = 0;

  for(i = 0; i < behaviors->num_values; i++) 
    class_counts[i] = 0;

#if USE_DURATION_COST > 0
  for(i = 0; i < behaviors->num_values; i++) 
    class_durations[i] = 0;
#endif

  for(i = 0; i < y->num_bouts; i++) {
    if(i) // count the number of times a transition from class c_p to c occurs
      class_transitions[y->bouts[i-1].behavior*behaviors->num_values + y->bouts[i].behavior]++;   
    if(!NONE_CLASS_HAS_NO_SCORE || y->bouts[i].behavior != 0)
      class_counts[y->bouts[i].behavior]++;     // count the number of times each class c occurs
#if USE_DURATION_COST > 0
    duration = y->bouts[i].end_frame-y->bouts[i].start_frame;
    duration_diff = 0;
#if USE_DURATION_COST > 1
    if (duration < min_frame_duration[y->bouts[i].behavior]) duration_diff = 10;
    else if (duration > max_frame_duration[y->bouts[i].behavior]) duration_diff = 10;
#else
    if (duration < min_frame_duration[y->bouts[i].behavior]) 
      duration_diff = min_frame_duration[y->bouts[i].behavior] - duration;
    else if (duration > max_frame_duration[y->bouts[i].behavior]) 
      duration_diff = duration - max_frame_duration[y->bouts[i].behavior];
    duration_diff *= duration_diff;
#endif
    if(!NONE_CLASS_HAS_NO_SCORE || y->bouts[i].behavior != 0)
      class_durations[y->bouts[i].behavior] += duration_diff;
#endif

    // The main appearance features used to compute the score of a bout for a given class.
    // Since the total score is the sum of the scores over bouts, we can simply add together
    // the features of the bouts with the same class label
    if(!NONE_CLASS_HAS_NO_SCORE || y->bouts[i].behavior != 0) {
      psi_bout(b, y->bouts[i].start_frame, y->bouts[i].end_frame, y->bouts[i].behavior, tmp_features, true, false);
      for(j = 0; j < num_bout_features; j++) 
	class_features[y->bouts[i].behavior*num_bout_features+j] += tmp_features[j];
    }
  }
        
  SparseVector retval = SparseVector(all_features, sizePsi);
  free(tmp_features);
  if(y != yy) delete y;
  if(y_partial) delete y_partial;

  return retval;
}

void SVMBehaviorSequence::print_features(FILE *fout, double *feat) {
  for(int i = 0; i < num_bout_features; i++)
    fprintf(fout, "    %lf %s\n", feat[i], bout_features[i].name);
}

void SVMBehaviorSequence::print_features(const char *fname, StructuredDataset *dataset, bool normalized) {
  double *tmp_features = (double*)malloc(sizeof(double)*num_bout_features);
  int i, j;
  BehaviorBoutSequence *y;
  BehaviorBoutFeatures *x;
  StructuredExample **ex = dataset->examples;
  FILE *fout = fopen(fname, "w");
  assert(fout);

  for(i = 0; i < dataset->num_examples; i++) {
    y = (BehaviorBoutSequence*)ex[i]->y;
    x = (BehaviorBoutFeatures*)ex[i]->x;
    
    for(j = 0; j < y->num_bouts; j++) {
      fprintf(fout, "  bout %d-%d %s\n", y->bouts[j].start_frame, y->bouts[j].end_frame, 
	      behaviors->values[y->bouts[j].behavior].name);
      psi_bout(x, y->bouts[j].start_frame, y->bouts[j].end_frame, y->bouts[j].behavior, tmp_features, normalized, false);
      print_features(fout, tmp_features);
    }
  }
  fclose(fout);
  free(tmp_features);
}


double *g_class_features, *g_class_transitions;
int g_num_features;
int cmp_feature_inds(const void * a, const void * b) { 
  int i1 = *(int*)a, i2 = *(int*)b;
  double f1 = i1 < g_num_features ? g_class_features[i1] : g_class_transitions[i1-g_num_features];
  double f2 = i2 < g_num_features ? g_class_features[i2] : g_class_transitions[i2-g_num_features];
  double d = my_abs(f1) - my_abs(f2);
  return (d < 0 ? 1 : (d > 0 ? -1 : 0)); 
}

void SVMBehaviorSequence::print_weights(FILE *fout, double *weights) {
  double *ptr = weights;
  double *class_features, *class_transitions;
  int i, j;
  int *inds = (int*)malloc(sizeof(int)*(num_bout_features+100000));

  
  class_features = ptr; ptr += num_bout_features*behaviors->num_values;
  class_transitions = ptr; ptr += behaviors->num_values*behaviors->num_values;
  g_num_features = num_bout_features;
  
  for(i = 0; i < behaviors->num_values; i++, class_features += num_bout_features, class_transitions += behaviors->num_values) {
    for(j = 0; j < num_bout_features+behaviors->num_values; j++)
      inds[j] = j;
    g_class_features = class_features;
    g_class_transitions = class_transitions;
    qsort(inds, num_bout_features+behaviors->num_values, sizeof(int), cmp_feature_inds);
    for(j = 0; j < num_bout_features+behaviors->num_values; j++) {
      if(inds[j] < num_bout_features)
	fprintf(fout, "%lf %s %s\n", class_features[inds[j]], behaviors->values[i].name, bout_features[inds[j]].name);
      else
	fprintf(fout, "%lf %s->%s\n", class_transitions[inds[j]-num_bout_features],behaviors->values[i].name, 
		behaviors->values[inds[j]-num_bout_features].name);
    }
  }
  free(inds);
}

void SVMBehaviorSequence::print_weights(const char *fname, double *weights) {
  FILE *fout = fopen(fname, "w");
  print_weights(fout, weights);
  fclose(fout);
}


bool SVMBehaviorSequence::SaveDataset(StructuredDataset *d, const char *fname, int start_from) {
  save_examples(fname, d);
  return true;
}

StructuredDataset *SVMBehaviorSequence::LoadDataset(const char *fname) {
  if(debugLevel > 0) fprintf(stderr, "Reading dataset %s...", fname);
  
  StructuredDataset *dataset = StructuredSVM::LoadDataset(fname);

  Lock();

  if(strlen(debugdir) && debug_features) 
    CreateDirectoryIfNecessary(debugdir);
	      

  int j;       
  //char **train_list = load_examples(fname, &num);
  bool computeClassTransitions = class_training_transitions ? false : true;
  //bool computeClassTransitions = true; // TEMP: EYRUN
  //StructuredDataset *dataset = new StructuredDataset();

  for(j = 0; j < dataset->num_examples; j++) {
    //StructuredExample *ex = read_struct_example(train_list[j], train_list[j], false);
    //dataset->AddExample(ex);
    BehaviorBoutSequence *y = (BehaviorBoutSequence*)dataset->examples[j]->y;
		
    // Temporary fix, ground truth segmentations should have bouts that span all the way to b->num_frames
    int T = ((BehaviorBoutFeatures*)dataset->examples[j]->x)->num_frames;
    if(y->num_bouts && y->bouts[y->num_bouts-1].end_frame==T-1)
      y->bouts[y->num_bouts-1].end_frame = T;
  }
  if(computeClassTransitions) {
    printf("Training examples %d behavior sequences\n", (int)dataset->num_examples);


    // Compute all feature cache data structures for all training examples.  This requires first computing some
    // statistics of the training set features for normalization purposes
    compute_feature_mean_variance_median_statistics(dataset);

    // Compute features for each training bout, normalized to have (0,1) mean and standard deviation
    for(j = 0; j < dataset->num_examples; j++) {
      BehaviorBoutFeatures *behavior_bout = (BehaviorBoutFeatures*)dataset->examples[j]->x;
      behavior_bout->fvec = Psi(dataset->examples[j]->x, dataset->examples[j]->y).ptr();
#if KEEP_FEATURES_IN_MEMORY == 0
      behavior_bout->Clear();
#endif
    }

    if(strlen(debugdir) && debug_features) {
      char fname[1000];
      CreateDirectoryIfNecessary(debugdir);
      sprintf(fname, "%s/features_unnormalized.txt", debugdir);
      print_features(fname, dataset, false);
      sprintf(fname, "%s/features_normalized.txt", debugdir);
      print_features(fname, dataset, true);
    }
  } /*else {
      for(j = 0; j < num; j++) {
      BehaviorBoutFeatures *x = (BehaviorBoutFeatures*)dataset->examples[j]->x;
      x->ComputeCaches(this);
      x->fvec = Psi(dataset->examples[j]->x, dataset->examples[j]->y).ptr();
      }
      }*/
  
  Unlock();
  
  return dataset;
}



// Keep track of the number of transitions between each pair of classes
void SVMBehaviorSequence::ComputeClassTrainsitionCounts() {
  int i, j, duration;
  bool allocate = !class_training_transitions;
  
  if(allocate) {
    class_training_transitions = (int**)malloc(behaviors->num_values*sizeof(int*));
    class_training_transitions_count = (int*)malloc(behaviors->num_values*sizeof(int));
    class_training_count = (int*)malloc(behaviors->num_values*sizeof(int));
  }
  for(i = 0; i < behaviors->num_values; i++) {
    class_training_transitions[i] = (int*)malloc(behaviors->num_values*sizeof(int));
    class_training_count[i] = class_training_transitions_count[i] = 0;
    for(j = 0; j < behaviors->num_values; j++) {
      class_training_transitions[i][j] = 0;
    }
  }
  for(j = 0; j < trainset->num_examples; j++) {
    BehaviorBoutFeatures *x = (BehaviorBoutFeatures*)trainset->examples[j]->x;
    BehaviorBoutSequence *y = (BehaviorBoutSequence*)trainset->examples[j]->y;
    for(i = 0; i < y->num_bouts; i++) {
      class_training_count[y->bouts[i].behavior]++;
#if USE_DURATION_COST > 0
      duration = y->bouts[i].end_frame - y->bouts[i].start_frame;
      min_frame_duration[y->bouts[i].behavior] = my_min(min_frame_duration[y->bouts[i].behavior], duration);
      max_frame_duration[y->bouts[i].behavior] = my_max(max_frame_duration[y->bouts[i].behavior], duration);
#endif
      if(i) {
	class_training_transitions[y->bouts[i].behavior][y->bouts[i-1].behavior]++;
	if(i==y->num_bouts-1 && y->bouts[i].behavior != NONE_BEHAVIOR && y->bouts[i].end_frame < x->num_frames)
	  class_training_transitions[NONE_BEHAVIOR][y->bouts[i].behavior]++;
      } else if(y->bouts[i].start_frame > 0 && y->bouts[i].behavior != NONE_BEHAVIOR)
	class_training_transitions[y->bouts[i].behavior][NONE_BEHAVIOR]++;
    }
  }

    // Compress the class transitions into a sparse array
  for(i = 0; i < behaviors->num_values; i++) {
    class_training_transitions_count[i] = 0;
    for(j = 0; j < behaviors->num_values; j++) {
      if(class_training_transitions[i][j])
	class_training_transitions[i][class_training_transitions_count[i]++] = j;
    }
  }
}

StructuredExample *SVMBehaviorSequence::read_struct_example(const char *label_fname, const char *features_fname, bool computeFeatures) {
  StructuredExample *ex = new StructuredExample;
  ex->x = NewStructuredData();
  ex->y = label_fname ? NewStructuredLabel(ex->x): NULL;
  BehaviorBoutFeatures *feature_cache = (BehaviorBoutFeatures*)ex->x;
  BehaviorBoutSequence *bouts = (BehaviorBoutSequence*)ex->y;

  strcpy(feature_cache->fname, features_fname);
  bool b = feature_cache->load(features_fname, this);
  assert(b);

  if(bouts) {
    strcpy(bouts->fname, label_fname);
    bool b = bouts->load(label_fname);
    assert(b);
  }

  if(computeFeatures) {
    feature_cache->ComputeCaches(this);
    feature_cache->fvec = label_fname ? Psi(ex->x, ex->y).ptr() : NULL;
  }
  return ex;
}





void SVMBehaviorSequence::OnFinishedIteration(StructuredData *x, StructuredLabel *y, SparseVector *w, StructuredLabel *ybar) {
  DebugExample(x, y, ybar, w);
  BehaviorBoutFeatures *b = (BehaviorBoutFeatures*)x;
  b->Clear();
}

void SVMBehaviorSequence::OnFinishedPassThroughTrainset() {
  Debug();
}

void SVMBehaviorSequence::Debug() {
  char fname[1000];
  int iter_num = this->t;
  if(strlen(debugdir) && debug_weights) {
    sprintf(fname, "%s/weights_%d.txt", debugdir, iter_num);
    SparseVector *w = GetCurrentWeights();
    double *ww = w->get_non_sparse<double>(sizePsi);
    print_weights(fname, ww);
    delete [] ww;
    delete w;
  }

  if(strlen(debugdir) && debug_model) {
    sprintf(fname, "%s/model_%d.txt", debugdir, iter_num);
    StructuredSVM::Save(fname, false);
  }
}

void SVMBehaviorSequence::DebugExample(StructuredData *x, StructuredLabel *y, StructuredLabel *ybar, SparseVector *w) {
  char fname[1000];
  int iter_num = this->t;
  double *ww = NULL;

  if(strlen(debugdir) && debug_predictions) {
    sprintf(fname, "%s/index.html", debugdir);
    FILE *fout = fopen(fname, "a");
    assert(fout);
    fprintf(fout, "<br><br><h2>Iteration %d</h2><a href=\"weights_%d.txt\">weights</a>|<a href=\"iter%d.html\">predictions</a>\n", iter_num, iter_num, iter_num);
    fclose(fout);
  
    // Save a visualization of all predicted bouts
    char *html=(char*)malloc(10000000), *html_gt=(char*)malloc(10000000), folder[1000], file[1000], fname[1000];
    ExtractFolderAndFileName(((BehaviorBoutFeatures*)x)->GetFileName(), folder, file);
    StripFileExtension(file);
    if(debug_features) {
      if(!ww) ww = w->get_non_sparse<double>(sizePsi);
      if(iter_num >= 0) sprintf(fname, "%s/%s_%d_features.html", debugdir, file, iter_num);
      else sprintf(fname, "%s/%s_features.html", debugdir, file);
      DebugFeatures(fname, (BehaviorBoutSequence*)ybar, ww);
    }
    if(iter_num >= 0) sprintf(fname, "%s/%s_%d", debugdir, file, iter_num);
    else sprintf(fname, "%s/%s", debugdir, file);
    ((BehaviorBoutSequence*)ybar)->Visualize(behaviors, fname, html);
    if(y) {
      BehaviorBoutSequence *yy = (BehaviorBoutSequence*)y;
      if(!yy->score) {
	if(ybar) yy->loss = Loss(y, ybar);
	
      }
      if(debug_features) {
	sprintf(fname, "%s/%s_gt_%d_features.html", debugdir, file, iter_num);
	DebugFeatures(fname, yy, ww);
      }
      sprintf(fname, "%s/%s_gt_%d", debugdir, file, iter_num);
      ((BehaviorBoutSequence*)y)->Visualize(behaviors, fname, html_gt);
    } 
    if(iter_num >= 0)
      sprintf(fname, "%s/iter%d.html", debugdir, iter_num);
    else
      sprintf(fname, "%s/index.html", debugdir);
    fout = fopen(fname, "a");
    assert(fout);
    BehaviorBoutSequence *yybar = (BehaviorBoutSequence*)ybar;
    if(y) {
      BehaviorBoutSequence *yy = (BehaviorBoutSequence*)y;
      fprintf(fout, "<br><br>%s: %s, score=%f, loss=%f<br>%s<br>%s: ground truth, score=%f, loss=%f<br>%s\n", file, iter_num >= 0 ? "most violated" : "best_score", (float)yybar->score, (float)yybar->loss, html, file, (float)yy->score, (float)yy->loss, html_gt);
    } else 
      fprintf(fout, "<br><br>%s: best score, score=%f, loss=%f<br>%s\n", file, (float)yybar->score, (float)yybar->loss, html);
    fclose(fout);
    free(html);
    free(html_gt);
  }
  if(ww) free(ww);
}

double **g_table = NULL;
BehaviorBout **g_states = NULL;
BehaviorBoutSequence *g_y = NULL;

char *getFilenameWithoutPath(char *filepath) {
	int lastChar = (int)strlen(filepath) - 1;

	for(int pos=lastChar; pos >=0; pos--) {
		if(filepath[pos]=='\\' || filepath[pos]=='/')
			return filepath+pos+1;
	}

	return filepath;
}




// Helper function for Inference().  Allocate memory for the label we are returning
void SVMBehaviorSequence::init_bout_label(BehaviorBoutSequence *ybar, BehaviorBoutSequence *y) {
  int sz = sizeof(BehaviorBout*);
  ybar->bouts = (BehaviorBout*)realloc(ybar->bouts, sz);
  memset(ybar->bouts, 0, sz);
  ybar->num_bouts = 0;
  ybar->slack = ybar->score = 0;
  ybar->loss = 0;
  if(y) {
    ybar->features = y->features;
    y->score = 0;
  }
}

// TODO: change this to be the hamming distance
// If a ground truth label y is given, update the componenent of loss(y,ybar) that is 
// attributable to the completed bout (c_prev,t_p,t)
double SVMBehaviorSequence::compute_updated_bout_loss(BehaviorBoutFeatures *b, BehaviorBoutSequence *y, 
						      int T, int t_p, int t, int c_prev, double *fn, 
						      int *gt_bout, double *dur_gt, double &loss_fp, double &loss_fn) {
  double loss_score = 0, fp;
  int TT = my_min(T,b->num_frames-1);  // int TT = T-1;

  if(y) { 
#if HAMMING_LOSS == 1
    loss_fp = loss_fn = 0;
    for(int j = gt_bout[t_p]; j <= gt_bout[t] && j < y->num_bouts; j++) {
      // Loop through all ground truth bouts intersecting with the completed bout
      double inter = (b->frame_times[my_min(my_min(y->bouts[j].end_frame,t),TT)] - 
                      b->frame_times[my_max(y->bouts[j].start_frame,t_p)]);
      if(inter <= 0) continue;
      if(c_prev != y->bouts[j].behavior) {
        loss_score += inter;
      }
      loss_fp = loss_fn = loss_score/2;
    }
#else
    double dur = b->frame_times[my_min(my_min(y->bouts[gt_bout[t]].end_frame,t),TT)]-b->frame_times[t_p];
    loss_fp = fp = match_false_positive_cost(dur, c_prev); // l^b_fp in (5) 
    loss_fn = 0;
    for(int j = gt_bout[t_p]; j <= gt_bout[t] && j < y->num_bouts; j++) {
      // Loop through all ground truth bouts intersecting with the completed bout
      double inter = (b->frame_times[my_min(my_min(y->bouts[j].end_frame,t),TT)] - 
                      b->frame_times[my_max(y->bouts[j].start_frame,t_p)]);
      if(inter <= 0) continue;
      if(c_prev == y->bouts[j].behavior) {
        // the absolute false negative cost is not used but rather the relative FN cost, 
        // where the "common shared" constant is left out, 
        // but since we are only interested in the maximum, that doesn't matter
        loss_fp -= fp*inter/dur; // ... and subtract agreeing frames from that maximum
      } else {
        loss_fn += fn[j]*inter/dur_gt[j];
      }
    }
    loss_score = loss_fn + loss_fp;
#endif
  } 

  return loss_score;
}


// Update the state portion of the dynamic programming cache tables
void SVMBehaviorSequence::store_solution(BehaviorBout &state, int t_p, int t, int c_prev, double bout_score, 
					 double transition_score, double unary_score, double loss_fn, 
					 double loss_fp, double duration_score) {
  // states[t][c_next] stores start/end/c_prev, such that we can backtrack to lookup the 
  // optimal solution corresponding to table[t][c_next]                                        
  state.start_frame =  t_p;
  state.end_frame =  t;
  state.behavior = c_prev;

  // The rest of this stuff is debug information
  state.bout_score = bout_score;                
  state.transition_score = transition_score; 
  state.unary_score = unary_score;
  state.loss_fn = loss_fn;
  state.loss_fp = loss_fp;
  state.duration_score = duration_score;
}


// Helper function for Inference() subject to partial labeling constraint
// Update transition counts (e.g., number of times a transition from a bout of behavior 
// c_p to a bout of behavior c_p) occurs in the training set using a partial label
void SVMBehaviorSequence::update_transition_counts_with_partial_label(BehaviorBoutSequence *y_partial, 
                 int** &old_class_transitions, int* &old_class_transition_counts, int* &old_class_training_counts) {
  int c, c_p, i, j;

  if(y_partial) {
    old_class_transition_counts = (int*)malloc(behaviors->num_values*(sizeof(int*)+sizeof(int)*(2+behaviors->num_values)));
    old_class_training_counts = old_class_transition_counts + behaviors->num_values;
    old_class_transitions = (int**)(old_class_training_counts + behaviors->num_values);
    int *ptr = (int*)(old_class_transitions + behaviors->num_values);
    for(c = 0; c < behaviors->num_values; c++) {
      old_class_transition_counts[c] = class_training_transitions_count[c];
      old_class_training_counts[c] = class_training_count[c];
      old_class_transitions[c] = ptr;
      for(int i = 0; i < class_training_transitions_count[c]; i++)
	old_class_transitions[c][i] = class_training_transitions[c][i];
      ptr += class_training_transitions_count[c];
    }

    for(i = 0; i < y_partial->num_bouts; i++) {
      c = y_partial->bouts[i].behavior;
      if(!class_training_count[c])
        class_training_count[c]++;
      if(i == 0 || y_partial->bouts[i].start_frame != y_partial->bouts[i-1].end_frame)
	continue;
      c_p = y_partial->bouts[i-1].behavior;
      
      if(c >= 0 && c_p >= 0) {
        for(j = 0; j < class_training_transitions_count[c_p]; j++)
          if(class_training_transitions[c][j] == c_p)
            break;
        if(j == class_training_transitions_count[c]) {
          class_training_transitions[c][class_training_transitions_count[c]++] = c_p;
	}
      }
    }
  }
}

int SVMBehaviorSequence::get_bout_start_time(int *durations, int &tt, int t_p, int t, int &next_duration, 
					     int &last_gt, int &last_partial, int *gt_bout, int *partial_label_bout, 
					     BehaviorBoutSequence *y, BehaviorBoutSequence *y_partial, 
					     int &restrict_c_prev, int &restrict_c_next, bool *allowable_time_frames) {
  bool isFirst = t_p == t;
  next_duration = 1;

  if(y) {
    // When given a groundtruth label y, stores the index of the bout corresponding to this timestep
    // in y->bouts
    if(isFirst) {
      last_gt = gt_bout[t] = gt_bout[t-1];
      while(y->bouts[gt_bout[t]].end_frame < t)
	last_gt = gt_bout[t] = gt_bout[t]+1;
    }
    while(y && (last_gt >= y->num_bouts || t_p <= y->bouts[last_gt].start_frame))
      last_gt--;
  }

  if(y_partial) {
    if(isFirst) {
      // Invoked the first time the dynamic programming algorithm gets to frame t.
      // Cache the index of the bout in the partial label that contains timestep t
      last_partial = partial_label_bout[t] = partial_label_bout[t-1] < y_partial->num_bouts && 
        y_partial->bouts[partial_label_bout[t-1]].end_frame <= t ? 
        partial_label_bout[t-1]+1 : partial_label_bout[t-1];

      // If the partial label contains a bout of a particular label at frame t, then it must be the case
      // the bout we are predicting also has a class c of that behavior.  So set restrict_c_next
      if(partial_label_bout[t] < y_partial->num_bouts && 
         y_partial->bouts[partial_label_bout[t]].start_frame <= t)
        restrict_c_next = y_partial->bouts[partial_label_bout[t]].behavior;
    }
  }

  int t_p_before = t_p;
  if(durations)
    t_p = t-durations[tt];
  else
    t_p = t-tt-1;

  // We may choose to add extra durations to make sure we explore the solutions
  // contained in y or b->partial_label
  if(y && last_gt >= 0 && (t_p < 0 ? -1 : gt_bout[t_p]) != last_gt && t_p != y->bouts[last_gt].start_frame && 
     (!allowable_time_frames || allowable_time_frames[y->bouts[last_gt].start_frame])) {
    t_p = y->bouts[last_gt].start_frame;
    next_duration = 0;
  }

  if(y_partial && y_partial->num_bouts) {
    if(last_partial < y_partial->num_bouts && t_p < y_partial->bouts[last_partial].start_frame && 
       t_p_before > y_partial->bouts[last_partial].start_frame) {
      t_p = y_partial->bouts[last_partial].start_frame;
      next_duration = 0;
    } else if(last_partial < y_partial->num_bouts && t_p < y_partial->bouts[last_partial].end_frame && 
       t_p_before > y_partial->bouts[last_partial].end_frame) {
      t_p = y_partial->bouts[last_partial].end_frame;
      next_duration = 0;
    }
    while(last_partial >= y_partial->num_bouts || (t_p < y_partial->bouts[last_partial].start_frame && last_partial > 0 && t_p < y_partial->bouts[last_partial-1].end_frame))
      last_partial--;
    if(last_partial < y_partial->num_bouts && t_p < y_partial->bouts[last_partial].start_frame && 
       t_p_before > y_partial->bouts[last_partial].start_frame) {
      t_p = y_partial->bouts[last_partial].start_frame;
      next_duration = 0;
    } else if(last_partial < y_partial->num_bouts && t_p < y_partial->bouts[last_partial].end_frame && 
       t_p_before > y_partial->bouts[last_partial].end_frame) {
      t_p = y_partial->bouts[last_partial].end_frame;
      next_duration = 0;
    }
  }

  tt += next_duration;
  if(t_p <= 0) {
    t_p = 0;
    next_duration = -1;
  }
  return t_p;
}

bool SVMBehaviorSequence::check_agreement_with_partial_label(BehaviorBoutSequence *y_partial, int t_p, int t,
							     int *partial_label_bout, int &restrict_c_prev) {
  if(y_partial) {
    // If the partial label contains a bout of a particular label at frame t_p, then it must be the case
    // the bout we are predicting also has a class c_prev of that behavior.  So set restrict_c_prev
    if(y_partial) {
      if(partial_label_bout[t_p] < y_partial->num_bouts && 
         y_partial->bouts[partial_label_bout[t_p]].start_frame <= t_p) {
        // If there are multiple bouts in  partial_label between t_p and t with different class labels, 
        // then it is definitely the case that the bout we are proposing doesn't agree with the partial label
        if(restrict_c_prev >= 0 && 
           y_partial->bouts[partial_label_bout[t_p]].behavior != restrict_c_prev)
          return false;
        restrict_c_prev = y_partial->bouts[partial_label_bout[t_p]].behavior;
      }
    }
  }

  return true;
}

// Helper function for Inference() subject to partial labeling constraint
void SVMBehaviorSequence::restore_transition_counts(BehaviorBoutSequence *y_partial, int** &old_class_training_transitions, 
						    int* &old_class_transition_counts, int* &old_class_training_counts) {
  if(y_partial && old_class_training_transitions) {
    for(int c = 0; c < behaviors->num_values; c++) {
      for(int i = 0; i < old_class_transition_counts[c]; i++)
	class_training_transitions[c][i] = old_class_training_transitions[c][i];
      class_training_transitions_count[c] = old_class_transition_counts[c];
      class_training_count[c] = old_class_training_counts[c];
    }
    free(old_class_transition_counts);
  }
}


// Helper function for Inference().  Extract the optimal solution after dynamic programming has run by backtracking
// through its cache tables
#if USE_DURATION_COST > 0
void SVMBehaviorSequence::backtrack_optimal_solution(BehaviorBoutSequence *ybar, double **table, 
						     BehaviorBout **states, double *duration_weights, int T, int c_prev, int t_p) {
#else
void SVMBehaviorSequence::backtrack_optimal_solution(BehaviorBoutSequence *ybar, double **table, 
						     BehaviorBout **states, int T, int c_prev, int t_p) {
#endif
  // First backtrack to count the number of bouts in the optimal solution  
  int t = T, c = 0; 
  ybar->num_bouts = 0;
#if USE_NEW_LOSS > 0
  if (c_prev > -1) {
    ybar->num_bouts++;
    t = t_p;
    c = c_prev;
  }
#endif
  while(t >= 0 && states[t][c].start_frame >= 0) { 
    int tt = states[t][c].start_frame;
    c = states[t][c].behavior;
    t = tt;
    ybar->num_bouts++; 
  }
  ybar->bouts = (BehaviorBout*)malloc(sizeof(BehaviorBout)*ybar->num_bouts);

  // Backtrack one more time to actually store that solution
  t = T; c = 0; 
  ybar->score = 0; 
  ybar->loss = 0;
  int i = ybar->num_bouts-1;
#if USE_NEW_LOSS > 0
  if (c_prev > -1) {
    ybar->bouts[i].start_frame = t_p;
    ybar->bouts[i].end_frame = T;
    ybar->bouts[i].behavior = c_prev;
    c = c_prev;
    t = t_p;
    i--;
  }
#endif
  while(t >= 0 && states[t][c].start_frame >= 0) {
#if USE_NEW_LOSS > 0 
   if (c_prev == -1){
#endif
    ybar->score += states[t][c].bout_score + states[t][c].transition_score + states[t][c].unary_score + states[t][c].duration_score;

    
    ybar->loss += states[t][c].loss_fn + states[t][c].loss_fp;
#if USE_NEW_LOSS > 0
   }
#endif
    ybar->bouts[i] = states[t][c];
    int tt = states[t][c].start_frame;
    c = states[t][c].behavior; 
    t = tt; 
    i--; 
  }
#if USE_NEW_LOSS > 0
 if (c_prev == -1)
#endif
}

void SVMBehaviorSequence::print_bout_sequence_scores(BehaviorBoutSequence *y) {
  double sum = 0;
  for(int i = 0; i < y->num_bouts; i++) {
    double score = y->bouts[i].bout_score + y->bouts[i].unary_score + y->bouts[i].transition_score + 
      y->bouts[i].loss_fn + y->bouts[i].loss_fp;
    sum += score;
    fprintf(stderr, "(%d %d), b=%d, sum=%lf, score=%lf, bout_score=%lf, unary_score=%lf, transition_score=%lf, loss_fn=%lf, loss_fp=%lf\n", 
	    y->bouts[i].start_frame, y->bouts[i].end_frame, y->bouts[i].behavior, sum, score,
	    y->bouts[i].bout_score, y->bouts[i].unary_score, y->bouts[i].transition_score, 
	    y->bouts[i].loss_fn, y->bouts[i].loss_fp);  
  }
}


 double SVMBehaviorSequence::get_duration_score(int duration, int c, double *duration_weights) {
   double duration_score = 0;
#if USE_DURATION_COST > 0
    double duration_diff = 0;
#if USE_DURATION_COST > 1
    if (duration < min_frame_duration[c]) duration_diff = 10;
    else if (duration > max_frame_duration[c]) duration_diff = 10;
#else
    if (duration < min_frame_duration[c]) 
	duration_diff = min_frame_duration[c] - duration;
    else if (duration > max_frame_duration[c]) 
	duration_diff = duration - max_frame_duration[c];
    duration_diff *= duration_diff; 
#endif
    duration_score = duration_weights[c] * duration_diff;
#endif
    if(NONE_CLASS_HAS_NO_SCORE && c == 0)
      duration_score = 0;
    return duration_score;
 }

// Helper function for Inference().  After dynamic programming runs, this does a bunch of sanity checks to test
// if the code anywhere has bugs
#if USE_DURATION_COST > 0
void SVMBehaviorSequence::sanity_check_dynamic_programming_solution(BehaviorBoutFeatures *b, BehaviorBoutSequence *ybar, 
				   BehaviorBoutSequence *y, SparseVector *w, double **class_weights, double **transition_weights, 
								    double *unary_weights, double *duration_weights, double **table, BehaviorBout **states, int T, BehaviorBoutSequence *y_partial) {
#else
void SVMBehaviorSequence::sanity_check_dynamic_programming_solution(BehaviorBoutFeatures *b, BehaviorBoutSequence *ybar, 
				   BehaviorBoutSequence *y, SparseVector *w, double **class_weights, double **transition_weights, 
				   double *unary_weights, double **table, BehaviorBout **states, int T, BehaviorBoutSequence *y_partial) {
#endif
  double *tmp_features = (double*)malloc((num_bout_features+1)*sizeof(double));

  double score = 0;
  double loss = 0;
  for(int i = 0; i < ybar->num_bouts; i++) {
    ybar->bouts[i].duration_score = get_duration_score(ybar->bouts[i].end_frame-ybar->bouts[i].start_frame, ybar->bouts[i].behavior, duration_weights);
    score += ybar->bouts[i].bout_score + ybar->bouts[i].transition_score + ybar->bouts[i].unary_score + ybar->bouts[i].duration_score;
    loss += ybar->bouts[i].loss_fn + ybar->bouts[i].loss_fp;
    
    // If this check fails, there is a problem with the basic dynamic programming algorithm
    assert(my_abs(score + loss - table[ybar->bouts[i].end_frame][i < ybar->num_bouts-1 ? ybar->bouts[i+1].behavior : 0]) <= .01); 
  }

  // Possibly redundant check: if this check fails, there is a problem with the basic dynamic programming algorithm
  assert(my_abs(ybar->score + ybar->loss - table[T][0]) <= .01); 


  // Make sure the score of ybar as accumulated over the dynamic programming algorithm is the same as the score
  // that is computed as the dot product between w and Psi(ybar).  If this fails, something is probably wrong with
  // the function Psi() or Inference(), such that they are inconsistent with each other
  SparseVector ybar_psi = Psi(b, ybar);
  double real_score = w->dot(ybar_psi);
  assert(my_abs(ybar->score - real_score) < .1);

  // y is the ground truth label when Inference() is invoked in learning mode
  if(y) {
    ybar->slack = ybar->score + ybar->loss;

    // Makes sure the computed components of the loss aggregated during dynamic programming are identical
    // to the loss when comparing the sequences y and ybar.  If this fails, something is probably wrong
    // with the computation of the loss during dynamic programming or with the function loss2()
    double l = loss2(y, ybar, 1);
    assert(my_abs(ybar->loss - l) < .01);


    // Compute the scores for the ground truth label y.  The cache tables in dynamic programming should always have
    // at least as high a score as the score yielded by the ground truth labeling.  If this fails, something
    // is probably wrong with the dynamic programming algorithm, such that it hasn't included the ground truth
    // label as a possible segmentation
    y->score = 0;
    for(int i = 0; i < y->num_bouts; i++) {
      psi_bout(b, y->bouts[i].start_frame, y->bouts[i].end_frame, -1, tmp_features, true, false);  
      y->bouts[i].bout_score = 0;
      for(int k = 0; k < num_bout_features; k++) 
        y->bouts[i].bout_score += class_weights[y->bouts[i].behavior][k]*tmp_features[k];

      if(i < y->num_bouts-1)
        y->bouts[i].transition_score = transition_weights[y->bouts[i].behavior][y->bouts[i+1].behavior];
      else
	y->bouts[i].transition_score = 0;
      y->bouts[i].unary_score = unary_weights[y->bouts[i].behavior];

      y->bouts[i].loss_fn = y->bouts[i].loss_fp = 0;
      y->bouts[i].duration_score = get_duration_score(y->bouts[i].end_frame-y->bouts[i].start_frame, y->bouts[i].behavior, duration_weights);
      y->score += y->bouts[i].bout_score + y->bouts[i].transition_score + y->bouts[i].unary_score + y->bouts[i].duration_score;

      // Making sure that ybar_score + ybar_loss >= y_score, so far
      if(!y->disable_checks && y->score > .01+(i < y->num_bouts-1 ? table[y->bouts[i+1].start_frame][y->bouts[i+1].behavior] : table[T][0]) && !y_partial ) {
        // Something went wrong, it might be informative for debugging (break here using gdb) to test the same dynamic 
        // programming problem but without using loss, then compare y_max to y
        g_table = table; g_states = states; g_y = y;
        assert(0);
      }
    }
    
    if(!y->disable_checks && !y_partial)
      assert(ybar->score+ybar->loss+.01 >= y->score);

    ybar->slack -= y->score;
    
    // Probably redundant with earlier checks, if this fails it means the dynamic programming algorithm hasn't 
    // correctly included the ground truth label as a possible segmentation
    if(!y->disable_checks && !y_partial)
      assert(ybar->slack >= -0.01);
  }

  free(tmp_features);
}


// Helper function for Inference().  To speedup Inference() when it is invoked during training, limit ourselves to 
// a restricted subset of the time frames.  That is, we only consider starting and ending bouts at particular
// time frames, which are chosen randomly.  As a consequence of calling this function, Inference() will take
// O(max_inference_learning_frames^2) time instead of O(T^2) time
//
// When Inference() is called at test time, it is assumed that we search all bout sizes exhaustively.  One could consider
// updating this function to also speedup test time (e.g. use a per-frame classifier to restrict certain frames from the
// search space).  Alternatively, one could speedup test time by breaking the entire video sequence into smaller segments
// (reducing T) 
bool *SVMBehaviorSequence::get_allowable_frame_times(BehaviorBoutSequence *y_gt, BehaviorBoutSequence *y_partial, int T) {
  if(!y_partial && (!y_gt || max_inference_learning_frames < 0 || T <= max_inference_learning_frames)) 
    return NULL;  // don't do approximate inference
  else {
    bool *allowable_time_frames = (bool*)malloc(sizeof(bool)*(T+1)*2+sizeof(int)*T);
    bool *allowable_tmp = allowable_time_frames+T+1;
    int *inds = (int*)(allowable_tmp+T+1);
    memset(allowable_time_frames, max_inference_learning_frames < 0 ? true : false, sizeof(bool)*(T+1));
    memset(allowable_tmp, true, sizeof(bool)*(T+1));
    allowable_time_frames[T] = true;
    
    // Ensure that the frames that bouts begin and end in y_gt and y_partial are included in the set of allowable frames
    // Ensure that the ground truth segmentation y_gt is included in the search space for Inference()
    if(y_gt) 
      for(int i = 0; i < y_gt->num_bouts; i++) 
	allowable_time_frames[y_gt->bouts[i].start_frame] = allowable_time_frames[y_gt->bouts[i].end_frame] = true;
    if(y_partial) {
      for(int i = 0; i < y_partial->num_bouts; i++) {
	// Ensure that the partial label is included in the search space for Inference()
	allowable_time_frames[y_partial->bouts[i].start_frame] = allowable_time_frames[y_partial->bouts[i].end_frame] = true;
	
      }
    }

    if(y_partial) {
      for(int i = 0; i < y_partial->num_bouts; i++) {
	// Never allow transitions within a bout in the partial label
	for(int j = y_partial->bouts[i].start_frame+1; j < y_partial->bouts[i].end_frame; j++) 
	  allowable_tmp[j] = allowable_time_frames[j] = false;
      }
      int num = 0;
      for(int i = 0; i < T; i++) {
	if(allowable_tmp[i])
	  inds[num++] = i;
      }
      
      // Now choose up to max_inference_learning_frames additional frames
      if(num) 
	for(int i = 0; i < max_inference_learning_frames; i++) 
	  allowable_time_frames[inds[rand()%num]] = true;
    } else {
      // Now choose up to max_inference_learning_frames additional frames
      for(int i = 0; i < max_inference_learning_frames; i++) 
	allowable_time_frames[rand()%T] = true;
    }

    return allowable_time_frames;
  }
}

void BehaviorBoutSequence::merge_bouts_with_same_behavior() {
  for(int i = 1; i < num_bouts; i++) {
    assert(bouts[i].start_frame >= bouts[i-1].end_frame);
  }

  int num = 1;
  for(int i = 1; i < num_bouts; i++) {
    if(bouts[i].behavior == bouts[num-1].behavior && 
       bouts[i].start_frame == bouts[num-1].end_frame)
      bouts[num-1].end_frame = bouts[i].end_frame;
    else
      bouts[num++] = bouts[i];
  }
  num_bouts = my_min(num_bouts, num);
  for(int i = 1; i < num_bouts; i++) {
    assert(bouts[i].start_frame >= bouts[i-1].end_frame);
  }
}

/*
void BehaviorBoutSequence::insert_bout(int ind, BehaviorBout b) {
  int i, j, k;
  bool inserted = false;
  for(i = 0; i < num_bouts[ind]; i++) {
    if(bouts[ind][i].end_frame > b.start_frame) {
      if(bouts[ind][i].end_frame <= b.end_frame) {
	assert(b.behavior == bouts[ind][i].behavior);
	bouts[ind][i].end_frame = b.end_frame;
	for(j = i; j < num_bouts[ind] && bouts[ind][j].start_frame < b.end_frame; j++) 
	  assert(b.behavior == bouts[ind][j].behavior);
	ptr = &bouts[ind][j];
	for(k = i; k < j; k++, ptr++) 
	  memcpy(&bouts[ind][k], ptr, sizeof(BehaviorBout));
	num_bouts[ind] = ;
      } else if(bouts[ind][i].start_frame < b.end_frame) {
	bouts[ind][i].start_frame = b.start_frame;
      } else {
	for(k = num_bouts[ind]; k > i; k--, ptr--) 
	  memcpy(&bouts[ind][k], &bouts[ind][k-1], sizeof(BehaviorBout));
	memcpy(&bouts[ind][i], &b, sizeof(BehaviorBout));
	num_bouts[ind]++;
      }
      inserted = true;
      break;
    }
  }
  if(!inserted) 
    memcpy(&bouts[ind][num_bouts[ind]++], &b, sizeof(BehaviorBout));
  merge_bouts_with_same_behavior();
}
*/

// If y_gt has unlabeled regions, the desired behavior is to ignore these regions during training.
// This is accomplished by setting the unlabeled regions set to the 'None' behavior, while at the same 
// time also setting the corresponding regions in y_partial to 'None'.  This forces all predictions to
// also have the 'None' behavior in these time frames, such that the bout score and features for both
// the ground truth label and predicted label are the same within these regions (thus the learning
// algorithm effectively ignores them).
bool SVMBehaviorSequence::fill_unlabeled_gt_frames(BehaviorBoutSequence *&y_gt, BehaviorBoutSequence *&y_partial) {
  if(y_gt) {
    // Check if there are any unlabeled frames in y_gt
    bool has_unlabeled = false;
    int lastframe = 0;
    for(int i = 0; i < y_gt->num_bouts; i++) {
      if(y_gt->bouts[i].start_frame > lastframe) {
	has_unlabeled = true;
	break;
      }
      lastframe = y_gt->bouts[i].end_frame;
    }
    if(!has_unlabeled)  // return if there are no unlabeled frames in y_gt
      return false;

    // Otherwise, create a copied version of y_gt with unlabeled regions set to the 'None' behavior.
    BehaviorBoutSequence *y_gt_new = (BehaviorBoutSequence*)NewStructuredLabel(y_gt->x);
    init_bout_label(y_gt_new, NULL);
    y_gt_new->filled = true;

    BehaviorBoutSequence *y_partial_new = NULL;
    if(!y_partial || !y_partial->filled) {
      y_partial_new = (BehaviorBoutSequence*)NewStructuredLabel(y_gt->x);
      init_bout_label(y_partial_new, NULL);
      y_partial_new->filled = true;
    }

    BehaviorBoutFeatures *x = (BehaviorBoutFeatures*)y_gt->x;
    
    int num_added = 0, partial_ind = 0, num_partial = 0;
    int s=0, e=0, b=0;
    lastframe = 0;
    for(int i = 0; i <= y_gt->num_bouts; i++) {
      if(i == y_gt->num_bouts) {
	if(e >= x->num_frames)
	  break;
	s = e = x->num_frames;
	b = 0;
      } else {
	s = y_gt->bouts[i].start_frame;
	e = y_gt->bouts[i].end_frame;
	b = y_gt->bouts[i].behavior;
      }

      if(y_partial && y_partial_new) {
	num_partial = y_partial->num_bouts;
	while(partial_ind < y_partial->num_bouts && y_partial->bouts[partial_ind].end_frame <= s) {
	  y_partial_new->bouts = (BehaviorBout*)realloc(y_partial_new->bouts, sizeof(BehaviorBout)*(num_partial+num_added+1));
	  y_partial_new->bouts[y_partial_new->num_bouts++] = y_partial->bouts[partial_ind++];
	}
      }
      if(s > lastframe || i == y_gt->num_bouts) {
	// Pad unlabelled bouts in y_gt as the "None" class
	y_gt_new->bouts = (BehaviorBout*)realloc(y_gt_new->bouts, sizeof(BehaviorBout)*(y_gt->num_bouts+num_added+1));
	y_gt_new->bouts[y_gt_new->num_bouts].start_frame = lastframe;
	y_gt_new->bouts[y_gt_new->num_bouts].end_frame = s;
	y_gt_new->bouts[y_gt_new->num_bouts++].behavior = 0;

	if(y_partial_new) {
	  // Set the corresponding regions in y_partial to also have the "None" class
	  y_partial_new->bouts = (BehaviorBout*)realloc(y_partial_new->bouts, sizeof(BehaviorBout)*(num_partial+num_added+1));
	  y_partial_new->bouts[y_partial_new->num_bouts].start_frame = lastframe;
	  y_partial_new->bouts[y_partial_new->num_bouts].end_frame = s;
	  y_partial_new->bouts[y_partial_new->num_bouts++].behavior = 0;
	}
	num_added++;
      }
      if(i != y_gt->num_bouts) {
	y_gt_new->bouts[y_gt_new->num_bouts].start_frame = s;
	y_gt_new->bouts[y_gt_new->num_bouts].end_frame = lastframe = e;
	y_gt_new->bouts[y_gt_new->num_bouts].behavior = b;
	y_gt_new->num_bouts++;
      }
    }
  
    if(y_gt_new) y_gt_new->merge_bouts_with_same_behavior();
    if(y_partial_new) y_partial_new->merge_bouts_with_same_behavior();

    y_gt = y_gt_new;
    if(y_partial_new) y_partial = y_partial_new;
    return true;
  }

  return false;
}

// Setup pointers to class_weights and transition weights.  It is assumed that
// the first class_weightsXnum_bout_features model parameters correspond to the
// class-feature weights, and the last num_classesXnum_classes weights correspond
// to the class transition weights
 void SVMBehaviorSequence::init_weight_pointers(double* &ptr, double** &class_weights, double** &transition_weights, double* &unary_weights, double* &duration_weights) {
   int i;
   class_weights = (double**)malloc(2*behaviors->num_values*sizeof(double*));
   transition_weights = class_weights+behaviors->num_values;
   for(i = 0; i < behaviors->num_values; i++, ptr += num_bout_features) 
     class_weights[i] = ptr;
   for(i = 0; i < behaviors->num_values; i++, ptr += behaviors->num_values) 
     transition_weights[i] = ptr; 
   unary_weights = ptr; 
   ptr += behaviors->num_values;
#if USE_DURATION_COST > 0
   duration_weights = ptr;
   ptr += behaviors->num_values;
#endif
}

/*
* Use dynamic programming to infer the label ybar that yields the highest score:
*   argmax_{ybar} w*psi(x,ybar)
* This is the predicted structured labelling of pattern x.
*
* When a groundtruth label y is specified (yy is non-NULL), it instead finds the label 
* ybar that has the highest score while also incorporating its loss with respect to
* the groundtruth label y
*   argmax_{ybar} loss(y,ybar)+w*psi(x,ybar)-w*psi(x,y)       // margin-scaling
*  or
*   argmax_{ybar} loss(y,ybar)*(1+w*psi(x,ybar)-w*psi(x,y))   // slack-scaling
* This is the most violated constraint with respect to example x for the current values of w
*
* It is also possible to supply a partial assignment to the desired output label, as specified in x->data.partial_label.
* When supplied a partial labelling, both ybar and y must agree entirely with the
* the partial labelling.  This is different from specifying a ground truth label y, where ybar is allowed to be
* different from y
*
*/
double SVMBehaviorSequence::Inference(StructuredData *x, StructuredLabel *y_bar, SparseVector *w, 
                                      StructuredLabel *yy_partial, StructuredLabel *y_gt, double w_scale) {
  BehaviorBoutSequence *ybar = (BehaviorBoutSequence*)(y_bar);
  BehaviorBoutSequence *y = (BehaviorBoutSequence*)(y_gt);
  BehaviorBoutSequence *y_partial = (BehaviorBoutSequence*)(yy_partial);
  BehaviorBoutFeatures *b = (BehaviorBoutFeatures*)x;
  double *tmp_features = (double*)malloc(2*(num_bout_features+1)*sizeof(double));
  double *tmp_features2 = tmp_features + (num_bout_features+1);
  double *ww = w->get_non_sparse<double>(sizePsi);
  double *ptr = ww;
  bool aborted = false;
  int T = b->num_frames;
  //if(y && y->num_bouts && y->num_bouts[0]) 
  //T = my_min(T, y->bouts[0][y->num_bouts[0]-1].end_frame);
  
  // If the ground truth label y has unlabeled regions, the desired behavior is to ignore these 
  // regions during training.  This is handled by altering y and y_partial
  fill_unlabeled_gt_frames(y, y_partial);
  ybar->filled = true;
  if(y) assert(y->bouts[y->num_bouts-1].end_frame == T);

  int *gt_bout = (int*)malloc(sizeof(int)*(T+1)*2);
  int *partial_label_bout = gt_bout + (T+1);
  bool *allowable_time_frames = get_allowable_frame_times(y, y_partial, T);
  double time_approx = time_approximation;

  if(!b->memory_buffer)
    b->ComputeCaches(this);

  // Initialize ybar
  if(y_gt) {
    sprintf(ybar->fname, "%s.pred.%d", b->fname, (int)this->t);
    int i = 1;
    while(FileExists(ybar->fname)) {
      sprintf(ybar->fname, "%s.pred.%d.%d", b->fname, (int)this->t, i++);
    }
  } else
    sprintf(ybar->fname, "%s.pred.model%d.iter%d", b->fname, modelId, iterId);
  init_bout_label(ybar, y);


  // During approximate inference, we can restrict searches over bout length to a subset of durations
  int *durations = (int*)malloc(sizeof(int)*(T+1)*4);
  int num_durations = 0;
  if(time_approximation) {
    for(int i = 1; i <= search_all_bout_durations_up_to; i++) // go only -time_approximation frames back
      durations[num_durations++] = i;
    if(time_approximation > 0) {
      // geometrically increasing series, where durations[k]= time_approximation^{k-1}
      double dur = search_all_bout_durations_up_to;
      while(dur <= T*(1+time_approximation)) {
	if(!num_durations || (int)(dur) != durations[num_durations-1])
	  durations[num_durations++] = (int)(dur);
	dur *= time_approximation;
      }
    }
  } else {
    for(int i = 1; i <= T; i++)
      durations[num_durations++] = i;
  }


  

  // Allocate buffers for dynamic programming
  // table[t][c] will store the maximum score for any sub-solution to frames 1...t in which 
  // a bout of class c begins at time t, and states[t][c] will store the corresponding bout labels
  double *bout_scores = (double*)malloc(behaviors->num_values*sizeof(double));
  double **class_weights, **transition_weights, *unary_weights;
#if USE_DURATION_COST > 0
  double *duration_weights;
#endif
  double *fn = y ? (double*)malloc(sizeof(double)*y->num_bouts*2) : NULL;
  double *dur_gt = fn ? fn + y->num_bouts : NULL;
  double **table = (double**)malloc((T+1)*(sizeof(double*)+behaviors->num_values*sizeof(double))), *ptr3;
  BehaviorBout **states = (BehaviorBout**)malloc((T+1)*(sizeof(BehaviorBout*)+behaviors->num_values*sizeof(BehaviorBout))), *ptr2;
  int **old_class_transitions = NULL, *old_class_transition_counts = NULL, *old_class_training_counts = NULL;
  int i;
  for(i = 0,  ptr3 = (double*)(table+T+1), ptr2 = (BehaviorBout*)(states+T+1); i <= T; 
      i++, ptr3 += behaviors->num_values, ptr2 += behaviors->num_values) {
    table[i] = ptr3;
    states[i] = ptr2;
  }

#if USE_DURATION_COST > 0
  init_weight_pointers(ptr, class_weights, transition_weights, unary_weights, duration_weights);
#else
  init_weight_pointers(ptr, class_weights, transition_weights, unary_weights, NULL);
#endif

  // When given a user supplied label, it may be the case that some class labels or label 
  // sequence in the partial labelling never appeared in the training set.  We add a couple of
  // checks here to protect against this case (otherwise our algorithm would find no solution)
  update_transition_counts_with_partial_label(y_partial, old_class_transitions, 
					      old_class_transition_counts, old_class_training_counts);
  
  // If evaluating loss with respect to a ground truth label y, compute the maximum false negative cost,
  // if every bout in y was missed entirely.  This loss will be subtracted off with each iteration of 
  // dynamic programming
  double max_fn_cost = 0;
  if(y) {
    int TT = my_min(T,b->num_frames-1);  
    for(i = 0; i < y->num_bouts; i++) {
      dur_gt[i] = (b->frame_times[my_min(y->bouts[i].end_frame,TT)] -
		   b->frame_times[y->bouts[i].start_frame]);
      fn[i] = match_false_negative_cost(dur_gt[i], y->bouts[i].behavior);
      max_fn_cost += fn[i];
    }
  }

  // Base case: initialize scores to 0
  for(int c = 0; c < behaviors->num_values; c++) {
    table[0][c] = 0;
    states[0][c].start_frame = states[0][c].end_frame = states[0][c].behavior = -1;
  }

  gt_bout[0] = 0;
  partial_label_bout[0] = 0;
  
  // Looping through all possible times t when an old bout could end and a new bout begins
  for(int t = 1; t <= T; t++) { 

    // Suppose a new bout of class c begins at time t.  Compute the optimal score for any labeling
    // through time 0...t that begins a bout of class c in time t.  We can prune the search space, 
    // because given the preceding completed bout (a bout beginning in some time step t_p and of some
    // class c_p) computation of the score function w.r.t. all frames at time t and onwards does not 
    // depend on any frames before t.  We can therefore enumerate all possible completed bouts and 
    // take the one with the highest score

    for(int c = 0; c < behaviors->num_values; c++) 
      table[t][c] = -INFINITY;
    
    // When given a manually supplied partial labelling, store the index of the bout corresponding to 
    // this timestep in b->partial_label->bouts
    int restrict_c_prev = -1, restrict_c_next = -1;

    // allowable_time_frames is a computational time saving trick, where we only allow bouts 
    // to start or end at a subset of allowable time frames
    if(allowable_time_frames && !allowable_time_frames[t]) {
      partial_label_bout[t] = partial_label_bout[t-1];
      gt_bout[t] = gt_bout[t-1];
      continue;   
    }

      
    // Looping through all possible times t_p when this bout begins (the bout we are considering
    // begins at t_p and ends at t
    bool is_first = true;
    int tt = 0, next_duration = 1, t_p = t;
    int last_partial, last_gt;
    while(next_duration >= 0 && tt < num_durations) {
	
      t_p = get_bout_start_time(durations, tt, t_p, t, next_duration, last_gt, last_partial, gt_bout, 
				partial_label_bout, y, y_partial, restrict_c_prev, restrict_c_next, allowable_time_frames);

      if(y && (relabelingExample||pauseWorkers) && t != T)  { // abort inference if training and someone wants to relabel an example
	aborted = true;
	continue;
      }

      // We can quickly discard all solutions where a candidate bout (t_p,t) overlaps a region in the partial
      // label in which there are multiple bouts of different classes.  Furthermore, we can compute
      // restrict_c_prev and restrict_c_next, which are restrictions on the possible labels of
      // c_prev and c_next, respectively
      if(!aborted && !check_agreement_with_partial_label(y_partial, t_p, t, partial_label_bout, restrict_c_prev))
	break;

      // allowable_time_frames is a computational time saving trick, where we only allow bouts 
      // to start or end at a subset of allowable time frames
      if(allowable_time_frames && !allowable_time_frames[t_p])
	continue; 
      
      // Compute bout-level features for the bout between t_p and t
      // For speed, CURRENTLY ASSUMING ALL CLASSES HAVE THE SAME BASIC FEATURE SPACE.  To get rid of this
      // assumption, move this below the for(c_prev...) loop and call psi_bout(b, t_p, t, c_prev, tmp_features).  
      // These loops were intentially ordered in a semi-funny way to make sure the computations
      // psi_bout(t_p,t) and w_c*psi() are computed as few times as possible
      // Compute bout scores for each class as the dot product between bout-level features and their 
      // corresponding weights
      psi_bout(b, t_p, t, -1, tmp_features, true, !is_first); 
      for(int c_prev = 0; c_prev < behaviors->num_values; c_prev++) {
	bout_scores[c_prev] = 0;
	if(!NONE_CLASS_HAS_NO_SCORE || c_prev != 0)
	  for(int k = 0; k < num_bout_features; k++) 
	    bout_scores[c_prev] += class_weights[c_prev][k]*tmp_features[k];
      }
        
      // Iterate through all classes c_next, where we are transitioning from a bout of class c_prev to a 
      // bout of class c_next, which is beginning at frame t.  
      for(int c_next = 0; c_next < (t == T ? 1 : behaviors->num_values); c_next++) {
	if(t < T && (!class_training_count[c_next] ||   // class c_next doesn't appear in the training set
		     (restrict_c_next >= 0 && c_next != restrict_c_next)))  // class c_next disagrees with the partial label
	  continue; 

	// Iterate through all classes c_prev, where we are transitioning from a bout of class c_prev to a
	// bout of class c_next, which is beginning at frame t.  Don't even consider class transition pairs from
	// c_prev to c_next that never occur in the training set
	for(int c_ind = 0; c_ind < (t == T ? behaviors->num_values : class_training_transitions_count[c_next]); c_ind++) { 
	  int c_prev = (t == T ? c_ind : class_training_transitions[c_next][c_ind]); 
	  if(!class_training_count[c_prev] ||     // class c_prev doesn't appear in the training set
	     (restrict_c_prev >= 0 && c_prev != restrict_c_prev))  // class c_prev disagrees with the partial label
	    continue; 
          
	  // Compute the score attributed to the proposed bout between (t_p,t) of label c_prev and transitioning to a 
	  // new bout c_next: 
	  //   bout_score: is the dot product between bout-level features and their corresponding weights
	  //   transition_score: is the score associated with transitioning from class c_prev to c_next
	  //   unary_score: unary prior score for each class
	  //   loss_score: if this function is invoked in traing mode, loss_score is the component of the 
	  //       segmentation loss Loss(y,ybar) that is attributable to the region from t_p to t
	  double loss_fn=0, loss_fp=0;
	  double bout_score = bout_scores[c_prev]; 
	  double transition_score = t != T ? transition_weights[c_prev][c_next] : 0;   
	  double unary_score = unary_weights[c_prev];
	  double loss_score = compute_updated_bout_loss(b, y, T, t_p, t, c_prev, fn, gt_bout, dur_gt, loss_fp, loss_fn);
#if USE_NEW_LOSS > 0
	  // compute loss differently
	  //BehaviorBoutSequence *ytemp = (BehaviorBoutSequence*)NewStructuredLabel(b);
	  backtrack_optimal_solution(ybar, table, states, duration_weights, t, c_prev, t_p); 
	  loss_score = Loss(y, ybar);
#endif
	  double duration_score = get_duration_score(t-t_p, c_prev, duration_weights);
	  double score = bout_score + transition_score + unary_score + loss_score + duration_score;       

	  // Check if the completed bout has a higher score than all previously examined solutions that 
	  // begin a bout of class c at time t
	  double f = table[t_p][c_prev] + score;
	  assert(!isnan(f));
	  if(f > table[t][c_next]) {
	    table[t][c_next] = f;
	    store_solution(states[t][c_next], t_p, t, c_prev, bout_score, transition_score, unary_score, loss_fn, loss_fp, duration_score);
	    //fprintf(stderr, "(t_p=%d, t=%d, c_prev=%d, c_next=%d), f=%f\n", t_p, t, c_prev, c_next, (float)f); 
	  } 
	  
	} // for(int c_ind = 0; ...),  c_next =...


	// If we don't search exhaustively through all possible bout durations, consider stretching the optimal solution
	// in which a behavior of class c_next begins at time t_p (it is preceded by a bout of some other behavior c_prev,
	// starting at some other time t_p_p),  such that the previous bout of behavior c_prev is stretched to end at time 
	// t instead of time t_p.  
	if(time_approx != 0 && t_p) {
	  int c_prev = states[t_p][c_next].behavior;
	  int t_p_p = states[t_p][c_next].start_frame;
	  if(c_prev >= 0 && table[t_p][c_next] > -INFINITY && (restrict_c_prev < 0 || (restrict_c_prev == c_prev && restrict_c_prev == c_next))) {
	    assert(states[t_p][c_next].end_frame == t_p);
	    
	    double bout_score = 0; 
	    if(!NONE_CLASS_HAS_NO_SCORE || c_prev != 0) {
	      psi_bout(b, t_p_p, t, -1, tmp_features2, true, !is_first);
	      for(int k = 0; k < num_bout_features; k++) 
		bout_score += class_weights[c_prev][k]*tmp_features2[k];
	    }
	    double loss_fn=0, loss_fp=0;
	    double transition_score = t != T ? transition_weights[c_prev][c_next] : 0;   
	    double unary_score = unary_weights[c_prev];
	    double loss_score = compute_updated_bout_loss(b, y, T, t_p_p, t, c_prev, fn, gt_bout, dur_gt, loss_fp, loss_fn);
#if USE_NEW_LOSS > 0
	    // compute loss differently
	    //BehaviorBoutSequence *ytemp = (BehaviorBoutSequence*)NewStructuredLabel(b);;
	    backtrack_optimal_solution(ybar, table, states, duration_weights, t, c_prev, t_p_p); 
	    loss_score = Loss(y, ybar);
#endif
	    double duration_score = get_duration_score(t-t_p_p, c_prev, duration_weights);
	    double score = bout_score + transition_score + unary_score + loss_score + duration_score;
	    double f = table[t_p_p][c_prev] + score;
	    if(f > table[t][c_next]) {
	      table[t][c_next] = f;
	      store_solution(states[t][c_next], t_p_p, t, c_prev, bout_score, transition_score, unary_score, loss_fn, loss_fp, duration_score);
	      //fprintf(stderr, "(t_prev=%d, t=%d, c_prev=%d, c_next=%d), f=%f\n", t_p, t, c_prev, c_next, (float)f); 
	    } 
	  }
	} // if(time_approx != 0 && t_p)
	
      } // for(int c_next = 0; ...)

      is_first = false;

    } // for(int tt; ...),  t_p = ...


    if(durations && NONE_CLASS_HAS_NO_SCORE) {
      // If the "none" behavior has no score, we can avoid doing an approximate inference, and exhaustively search through 
      // all possible durations for the "none" behavior.  
      // Currently, this isn't implemented in the fastest way possible (we simply replicated the code above, removed the 
      // part where we computed bout features, and restricted c_prev=0).  In the case where y=NULL and y_partial=NULL (e.g., 
      // regular inference, we could instead just cache in an array the highest scoring time to start a none behavior:
      //   t_opt[t]=argmax_{t'<=t}table[t_p][0]
      // and use it to update
      //   for all c_next: 
      //     f = (table[t_opt[t-1]][0]+transition_cost[c_next]["none"]);
      //     if(f > table[t][c_next]) store_solution(states[t][c_next], t_opt[t-1]][0], t, 0,...);
      //   t_opt[t] = max(t_opt[t-1],table[t][0])
      int tt = 0, next_duration = 1, t_p = t;
      int last_partial, last_gt;
      while(next_duration >= 0 && t_p >= 0) {
	// Set durations=NULL to exhaustively consider all durations
	t_p = get_bout_start_time(NULL, tt, t_p, t, next_duration, last_gt, last_partial, gt_bout, 
				  partial_label_bout, y, y_partial, restrict_c_prev, restrict_c_next, allowable_time_frames);
	
	// We can quickly discard all solutions where a candidate bout (t_p,t) overlaps a region in the partial
	// label in which there are multiple bouts of different classes.  Furthermore, we can compute
	// restrict_c_prev and restrict_c_next, which are restrictions on the possible labels of
	// c_prev and c_next, respectively
	if(!check_agreement_with_partial_label(y_partial, t_p, t, partial_label_bout, restrict_c_prev))
	  break;

	// allowable_time_frames is a computational time saving trick, where we only allow bouts 
	// to start or end at a subset of allowable time frames
	if(allowable_time_frames && !allowable_time_frames[t_p])
	  continue; 
       
        
	// Iterate through all classes c_next, where we are transitioning from a bout of class c_prev to a 
	// bout of class c_next, which is beginning at frame t.  
	for(int c_next = 0; c_next < (t == T ? 1 : behaviors->num_values); c_next++) {
	  if(t < T && (!class_training_count[c_next] ||   // class c_next doesn't appear in the training set
		       (restrict_c_next >= 0 && c_next != restrict_c_next)))  // class c_next disagrees with the partial label
	    continue; 
	  
	  // Select c_prev=0, the none behavior
	  for(int c_ind = 0; c_ind < (t == T ? behaviors->num_values : class_training_transitions_count[c_next]); c_ind++) {
	    int c_prev = (t == T ? c_ind : class_training_transitions[c_next][c_ind]); 
	    if(c_prev != 0 || !class_training_count[c_prev] ||     // class c_prev doesn't appear in the training set
	       (restrict_c_prev >= 0 && c_prev != restrict_c_prev))  // class c_prev disagrees with the partial label
	      continue; 
	    

	    double loss_fn=0, loss_fp=0;
	    double bout_score = 0;
	    double transition_score = t != T ? transition_weights[c_prev][c_next] : 0;   
	    double unary_score = unary_weights[c_prev];
	    double loss_score = compute_updated_bout_loss(b, y, T, t_p, t, c_prev, fn, gt_bout, dur_gt, loss_fp, loss_fn);
#if USE_NEW_LOSS > 0
	    // compute loss differently
	    //BehaviorBoutSequence *ytemp = (BehaviorBoutSequence*)NewStructuredLabel(b);
	    backtrack_optimal_solution(ybar, table, states, duration_weights, t, c_prev, t_p); 
	    loss_score = Loss(y, ybar);
#endif
	    double duration_score = get_duration_score(t-t_p, c_prev, duration_weights);
	    double score = bout_score + transition_score + unary_score + loss_score + duration_score;       
	    
	    // Check if the completed bout has a higher score than all previously examined solutions that 
	    // begin a bout of class c at time t
	    double f = table[t_p][c_prev] + score;
	    assert(!isnan(f));
	    if(f > table[t][c_next]) {
	      table[t][c_next] = f;
	      store_solution(states[t][c_next], t_p, t, c_prev, bout_score, transition_score, unary_score, loss_fn, loss_fp, duration_score);
	      //fprintf(stderr, "(t_p=%d, t=%d, c_prev=%d, c_next=%d), f=%f\n", t_p, t, c_prev, c_next, (float)f); 
	    } 
	  }
	} // for(int c_ind = 0; ...),  c_next =...
      }
    }
    

  } // for(int t = 0; ...)


#if USE_DURATION_COST > 0
  // Backtrack through table and states to extract the optimal solution
  backtrack_optimal_solution(ybar, table, states, duration_weights, T);
  if(!aborted) sanity_check_dynamic_programming_solution(b, ybar, y, w, class_weights, transition_weights, unary_weights, duration_weights, table, states, T, y_partial);
#else
  // Backtrack through table and states to extract the optimal solution
  backtrack_optimal_solution(ybar, table, states, T);
  if(!aborted) sanity_check_dynamic_programming_solution(b, ybar, y, w, class_weights, transition_weights, unary_weights, table, states, T, y_partial);
#endif
  
  // Restore modified transition tables, if necessary
  restore_transition_counts(y_partial, old_class_transitions, old_class_transition_counts, old_class_training_counts);
      
  // Cleanup
  free(table);
  free(states);
  free(class_weights);
  free(bout_scores);
  if(fn) free(fn);

  free(tmp_features);
  free(ww);
  free(gt_bout);
  free(durations);
  if(allowable_time_frames)
    free(allowable_time_frames);


  //if(y) { fprintf(stderr, "\ny:\n"); print_bout_sequence_scores(y,0); }
  //if(y_partial) { fprintf(stderr, "\ny_partial:\n"); print_bout_sequence_scores(y_partial,0); }
  //if(ybar) { fprintf(stderr, "\nybar:\n"); print_bout_sequence_scores(ybar,0); }

  if(y && y != y_gt) delete y;
  if(y_partial && y_partial != yy_partial) delete y_partial;

  return ybar->score + ybar->loss;
}
 
// Take a behavior sequence y_src and remove all bout labels between times t_start to t_end, such that labels
// in that section will be latent
BehaviorBoutSequence *SVMBehaviorSequence::bout_sequence_remove_section(BehaviorBoutSequence *y_src, int t_start, int t_end) {
  BehaviorBoutSequence *y_dst = (BehaviorBoutSequence*)NewStructuredLabel(y_src->x);
  bool changed = false;
  *y_dst = *y_src;
  y_dst->bouts = NULL;
  init_bout_label(y_dst, NULL);

  //Json::Value yy = y_src->save(this);
  //y_dst->load(yy, this);

  y_dst->num_bouts = 0;
  y_dst->bouts = (BehaviorBout*)realloc(y_dst->bouts, sizeof(BehaviorBout)*(y_src->num_bouts+2));
  for(int i = 0; i < y_src->num_bouts; i++) {
    if(y_src->bouts[i].start_frame < t_start || y_src->bouts[i].end_frame > t_end) {
      y_dst->bouts[y_dst->num_bouts] = y_src->bouts[i];
      if(y_dst->bouts[y_dst->num_bouts].start_frame < t_start && 
	 y_dst->bouts[y_dst->num_bouts].end_frame > t_start) {
	if(y_dst->bouts[y_dst->num_bouts].end_frame > t_end) {
	  y_dst->bouts[y_dst->num_bouts++].end_frame = t_start;
	  y_dst->bouts[y_dst->num_bouts] = y_src->bouts[i];
	  y_dst->bouts[y_dst->num_bouts].start_frame = t_end;
	} else
	  y_dst->bouts[y_dst->num_bouts].end_frame = t_start;
	changed = true;
      } else if(y_dst->bouts[y_dst->num_bouts].end_frame > t_end &&
		y_dst->bouts[y_dst->num_bouts].start_frame < t_end) {
	y_dst->bouts[y_dst->num_bouts].start_frame = t_end;
	changed = true;
      }
      y_dst->num_bouts++;
    } else
      changed = true;
  }
  
  if(!changed) {
    delete y_dst;
    return NULL;
  }
  return y_dst;
}



double SVMBehaviorSequence::ImportanceSample(StructuredData *x, SparseVector *w, StructuredLabel *y_gt, 
					     struct _SVM_cached_sample_set *set, double w_scale) {
  SparseVector *w_curr = GetCurrentWeights(true);
  double retval = 0;

  if(importance_sample_interval_size) {
    int t_start = 0;
    BehaviorBoutFeatures *b = (BehaviorBoutFeatures*)x;
    BehaviorBoutSequence *y = (BehaviorBoutSequence*)y_gt;
    int T = b->num_frames;
    //if(y && y->num_bouts && y->num_bouts[0]) 
    //  T = my_min(T, y->bouts[0][y->num_bouts[0]-1].end_frame);
    y->disable_checks = true;
    while(t_start < T && !relabelingExample && !pauseWorkers) {
      // Consider selecting the segmentation with the highest slack among all segmentations for this example that
      // have been considered so far.  We will see if there's some local change to this segmentation that will increase
      // the slack even more
      set->score_gt = w_curr->dot(*set->psi_gt);
      double best_score = set->score_gt;
      BehaviorBoutSequence *best_label = (BehaviorBoutSequence*)y_gt;
      for(int i = 0; i < set->num_samples; i++) {
	double score = w_curr->dot(*set->samples[i].psi) + set->samples[i].loss;
	set->samples[i].slack = score-set->score_gt;
	if(score > best_score) {
	  best_score = score;
	  best_label = (BehaviorBoutSequence*)set->samples[i].ybar;
	}
      }

      // Run inference on a subset of the sequence between t=t_start and t=t_start+importance_sample_interval_size
      BehaviorBoutSequence *ybar_partial = bout_sequence_remove_section(best_label, t_start, 
									my_min(t_start+importance_sample_interval_size,T));
      if(!ybar_partial) {
	t_start += importance_sample_interval_size;  // no labeled bouts in this particular region
	continue;
      }

      BehaviorBoutSequence *ybar = (BehaviorBoutSequence*)NewStructuredLabel(x);
      retval = Inference(x, ybar, w_curr, ybar_partial, y_gt, 1);
      sprintf(ybar->fname+strlen(ybar->fname), "_%d", t_start);
      fprintf(stderr, ".");
      if(retval > set->score_gt) { 
	Lock();
	SVM_cached_sample_set_add_sample(set, ybar);
	SVM_cached_sample_set_compute_features(set, trainset->examples[set->i]);
      
	// Use the current sample set to update the weights immediately.  Usually we would never do this, but
	// this allows us to update the model weights more frequently when bout sequences are really
	// long and inference is slow
	MultiSampleUpdate(set, trainset->examples[set->i], 1);
	delete w_curr;
	w_curr = GetCurrentWeights(false);
	Unlock();
      } else
	delete ybar;

      delete ybar_partial;
      t_start += importance_sample_interval_size;
    }
    y->disable_checks = false;
  }

  // Now run the full exhaustive Inference() procedure to find the most violated constraint
  StructuredLabel *ybar = NewStructuredLabel(x);
  retval = Inference(x, ybar, w_curr, NULL, y_gt, 1);
  Lock();
  set->score_gt = w_curr->dot(*set->psi_gt);
  SVM_cached_sample_set_add_sample(set, ybar);
  SVM_cached_sample_set_compute_features(set, trainset->examples[set->i]);
  if(set->num_samples) {
    SVM_cached_sample s = set->samples[0];
    set->samples[0] = set->samples[set->num_samples-1];
    set->samples[set->num_samples-1] = s;
  }
  Unlock();

  delete w_curr;

  return retval;
}



double      SVMBehaviorSequence::Loss(StructuredLabel *y_gt,  StructuredLabel *y_pred) {
  double l = 0;
  BehaviorBoutSequence *y = (BehaviorBoutSequence*)y_gt, *y_partial=NULL;
	
  fill_unlabeled_gt_frames(y, y_partial);

  l = loss2(y, y_pred, 0);
	
  if(y != y_gt) delete y;
  if(y_partial) delete y_partial;

  return l;
}

double      SVMBehaviorSequence::loss2(StructuredLabel *y_gt,  StructuredLabel *y_pred, int debug) {
  /* loss for correct label y and predicted label ybar. The loss for
     y==ybar has to be zero. sparm->loss_function is set with the -l option. */
  BehaviorBoutSequence *y = (BehaviorBoutSequence*)y_gt;
  BehaviorBoutSequence *ybar = (BehaviorBoutSequence*)y_pred;
  int curr_y=0, curr_ybar=0;  // current bout in y and ybar
  double l = 0;               // total loss so far
  double dur_y=0, dur_ybar=0; // duration in frames of the current bout
  double inter;               // intersection time between bouts
  double sum_y = 0, sum_ybar = 0;
  BehaviorBoutFeatures *b = y->features;
  int T = b->num_frames;
  if(y && y->num_bouts) 
    T = my_min(T, y->bouts[y->num_bouts-1].end_frame);
  int TT = my_min(T,b->num_frames-1);  
  double cl, l_fn = 0, l_fn2 = 0;

#if HAMMING_LOSS == 1
  while((curr_y < y->num_bouts || curr_ybar < ybar->num_bouts)) {
    if(curr_y < y->num_bouts && curr_ybar < ybar->num_bouts) {
      // Check if the current bout in y and ybar match.  
      dur_ybar = (b->frame_times[my_min(ybar->bouts[curr_ybar].end_frame,TT)] - 
		  b->frame_times[ybar->bouts[curr_ybar].start_frame]);
      dur_y = (b->frame_times[my_min(y->bouts[curr_y].end_frame,TT)] -
	       b->frame_times[y->bouts[curr_y].start_frame]);
      inter = (b->frame_times[my_min(my_min(y->bouts[curr_y].end_frame, ybar->bouts[curr_ybar].end_frame),TT)] - 
	       b->frame_times[my_max(y->bouts[curr_y].start_frame,ybar->bouts[curr_ybar].start_frame)]);
      assert(inter >= 0);
      if(y->bouts[curr_y].behavior != ybar->bouts[curr_ybar].behavior) 
	l_fn += inter;
    } else
      inter = 0;
    
    if(curr_ybar == ybar->num_bouts || (curr_y < y->num_bouts && 
					(y->bouts[curr_y].end_frame <= ybar->bouts[curr_ybar].end_frame))) {
      // Go to the next bout in y adding the appropriate loss based on whether or not
      // the bout was matched to a bout in ybar
      curr_y++;
    } else {
      // Go to the next bout in ybar, adding the appropriate loss based on whether or not
      // the bout was matched to a bout in y
      l += l_fn;
      ybar->bouts[curr_ybar].loss_fn = l_fn/2;
      ybar->bouts[curr_ybar].loss_fp = l_fn/2;
      curr_ybar++;
      l_fn = 0;
    }
  }
#else
// 	// TEMP EYRUN
// 	l = 0; 
// 	for(curr_ybar = 0; curr_ybar < ybar->num_bouts[0]; curr_ybar++)
// 		l += ybar->bouts[0][curr_ybar].loss_fp + ybar->bouts[0][curr_ybar].loss_fn;
// 	return l;
#if USE_NEW_LOSS > 0
  int beh_counts_gt[bevaiors->num_values], beh_counts_pred[behaviors->num_values];
  for (int i=0; i<behaviors->num_values; i++) {
    beh_counts_gt[i] = 0;
    beh_counts_pred[i] = 0;
  }
  for (int i=0; i<y->num_bouts; i++) {
    beh_counts_gt[y->bouts[i].behavior] ++;}
  for (int i=0; i<ybar->num_bouts; i++) {
    beh_counts_pred[ybar->bouts[i].behavior] ++;}
  for (int i=0; i<num_classes; i++) {
    beh_counts_gt[i] = my_max(beh_counts_gt[i],1);
    beh_counts_pred[i] = my_max(beh_counts_pred[i],1);
  }
#endif

  while((curr_y < y->num_bouts || curr_ybar < ybar->num_bouts)) {
    if(curr_y < y->num_bouts && curr_ybar < ybar->num_bouts) {
      // Check if the current bout in y and ybar match.  
      dur_ybar = (b->frame_times[my_min(ybar->bouts[curr_ybar].end_frame,TT)] - 
		  b->frame_times[ybar->bouts[curr_ybar].start_frame]);
      dur_y = (b->frame_times[my_min(y->bouts[curr_y].end_frame,TT)] -
	       b->frame_times[y->bouts[curr_y].start_frame]);
      inter = (b->frame_times[my_min(my_min(y->bouts[curr_y].end_frame, ybar->bouts[curr_ybar].end_frame),TT)] - 
	       b->frame_times[my_max(y->bouts[curr_y].start_frame,ybar->bouts[curr_ybar].start_frame)]);
      assert(inter >= 0);
      if(y->bouts[curr_y].behavior == ybar->bouts[curr_ybar].behavior) {
	sum_ybar += inter;
	sum_y += inter;
	//if(debug)
#if USE_NEW_LOSS > 0 
	l_fn -= match_false_negative_cost(dur_y, y->bouts[curr_y].behavior)*(dur_y ? (inter/dur_y) : 0) / beh_counts_gt[y->bouts[curr_y].behavior];
#else 
	l_fn -= match_false_negative_cost(dur_y, y->bouts[curr_y].behavior)*(dur_y ? (inter/dur_y) : 0);
# endif
      } else
#if USE_NEW_LOSS > 0
	l_fn2 += match_false_negative_cost(dur_y, y->bouts[curr_y].behavior)*(dur_y ? (inter/dur_y) : 0) / beh_counts_gt[y->bouts[curr_y].behavior];
#else
      l_fn2 += match_false_negative_cost(dur_y, y->bouts[curr_y].behavior)*(dur_y ? (inter/dur_y) : 0);
#endif
    } else
      inter = 0;

    if(curr_ybar == ybar->num_bouts || (curr_y < y->num_bouts && 
					(y->bouts[curr_y].end_frame <= ybar->bouts[curr_ybar].end_frame))) {
      // Go to the next bout in y adding the appropriate loss based on whether or not
      // the bout was matched to a bout in ybar
      cl = match_false_negative_cost(dur_y, y->bouts[curr_y].behavior)*(dur_y ? ((dur_y-sum_y)/dur_y) : 1);
#if USE_NEW_LOSS > 0
      cl /= beh_counts_gt[y->bouts[curr_y].behavior];
#endif
      l += cl; 
      if(debug == 2) fprintf(stderr, "y[%d]->%f %f\n", curr_y, cl, l);
      curr_y++;
      sum_y = 0;
    } else {
      // Go to the next bout in ybar, adding the appropriate loss based on whether or not
      // the bout was matched to a bout in y
      cl = match_false_positive_cost(dur_ybar, ybar->bouts[curr_ybar].behavior) *
	(dur_ybar ? ((dur_ybar-sum_ybar)/dur_ybar) : 1); 
#if USE_NEW_LOSS > 0
      cl /= beh_counts_pred[ybar->bouts[curr_ybar].behavior];
#endif
      l += cl;
      if(debug == 2) 
	fprintf(stderr, "ybar[%d]->%f %f\n", curr_ybar, cl, l);
      if(debug) {
	assert(my_abs(cl - ybar->bouts[curr_ybar].loss_fp) < .00001);
	assert(my_abs(l_fn2 - ybar->bouts[curr_ybar].loss_fn) < .00001);
      }
      ybar->bouts[curr_ybar].loss_fn = l_fn2;
      ybar->bouts[curr_ybar].loss_fp = cl;
      curr_ybar++;
      sum_ybar = 0;
      l_fn = 0; l_fn2 = 0;
    }
  }
#endif //HAMMING_LOSS == 1
  if(debug)
    assert(my_abs(l-ybar->loss) < .01);

  return l;
}




Json::Value SVMBehaviorSequence::Save() {
  Json::Value root, num_cl, trans;
  int i, j;
  root["version"] = INST_VERSION;

  if(behaviors)
    root["behaviors"] = SaveBehaviorDefinitions();

  Json::Value t;
  for(i = 0; i < behaviors->num_values; i++) {
    Json::Value o, tt;
    o["count"] = class_training_count[i];
    o["transitions_count"] = class_training_transitions_count[i];
    for(j = 0; j < class_training_transitions_count[i]; j++) {
      tt[j] = class_training_transitions[i][j];
    }
    o["transitions"] = tt;
#if USE_DURATION_COST > 0
    Json::Value limits;
    int idx = 0;
    limits[idx] = (int)min_frame_duration[i];
    idx = 1;
    limits[idx] = (int)max_frame_duration[i];
    o["limits"] = limits;
#endif
    t[i] = o;
  }
  root["transitions"] = t;

  root["smoothness_window"] = feature_sample_smoothness_window;

  if(bout_features) 
    root["bout_feature_params"] = SaveBoutFeatureParams();
	
  if(frame_features) 
    root["frame_feature_params"] = SaveFrameFeatureParams();
	
  Json::Value a(Json::arrayValue), b(Json::arrayValue);
  for(i=0; i<behaviors->num_values; i++){
    a[i] = false_negative_cost[i];
    b[i] = false_positive_cost[i];
  }
  root["false_negative_cost"] = a;
  root["false_positive_cost"] = b;
	
       
  return root;
}




bool SVMBehaviorSequence::Load(const Json::Value &root) {
  int i, j;

  if(root.isMember("version"))
    assert(!strcmp(root["version"].asString().c_str(), INST_VERSION));

  if(root.isMember("behaviors")) 
    if(!LoadBehaviorDefinitions(root["behaviors"]))
      return false;

  if(root.isMember("frame_feature_params")) 
    LoadFrameFeatureParams(root["frame_feature_params"]);
  if(root.isMember("bout_feature_params")) 
    LoadBoutFeatureParams(root["bout_feature_params"]);

  Init(behaviors, frame_features, num_frame_features, bout_features, num_bout_features);

  if(behaviors && root.isMember("transitions")) {
    class_training_transitions = (int**)malloc(behaviors->num_values*sizeof(int*));
    class_training_transitions_count = (int*)malloc(behaviors->num_values*sizeof(int));
    class_training_count = (int*)malloc(behaviors->num_values*sizeof(int));
    memset(class_training_transitions_count, 0, behaviors->num_values*sizeof(int));
    memset(class_training_count, 0, behaviors->num_values*sizeof(int));
    for(i = 0; i < behaviors->num_values; i++) {
      class_training_transitions[i] = (int*)malloc(behaviors->num_values*sizeof(int));
      memset(class_training_transitions[i], 0, behaviors->num_values*sizeof(int));
    }

    Json::Value t = root["transitions"];
    for(i = 0; i < behaviors->num_values; i++) {
      Json::Value o = t[i]; 
      Json::Value tt = o["transitions"];
      class_training_transitions_count[i] = o["transitions_count"].asInt();
      class_training_count[i] = o["count"].asInt();
      for(j = 0; j < class_training_transitions_count[i]; j++) {
	class_training_transitions[i][j] = tt[j].asInt();
      }
#if USE_DURATION_COST > 0
      Json::Value limits = o["limits"];
      int idx = 0;
      min_frame_duration[i] = limits[idx].asInt();
      idx = 1;
      max_frame_duration[i] = limits[idx].asInt();
#endif
    }
  }

  feature_sample_smoothness_window = root.get("smoothness_window",0).asInt();
  
  
  if(root.isMember("false_positive_cost") && root["false_positive_cost"].isArray()) {
    Json::Value a = root["false_negative_cost"], b=root["false_positive_cost"];
    for(i=0; i<behaviors->num_values; i++){
      false_negative_cost[i] = a[i].asDouble();
      false_positive_cost[i] = b[i].asDouble();
    }
  }

  sizePsi = getPsiSize(this->num_bout_features, behaviors->num_values);
  fprintf(stderr, "Psi is %d-dimensional, over %d frame features, %d bout features, %d behaviors\n", 
	  sizePsi, num_frame_features, num_bout_features, behaviors->num_values);

  return true;
}

bool SVMBehaviorSequence::LoadBehaviorDefinitions(const Json::Value &p) {
  FreeBehaviorDefinitions();

  int noneInd;
  Behaviors *behaviors = (Behaviors*)malloc(sizeof(Behaviors)); 
  behaviors->num_values = p.size();
  behaviors->values = (Behavior*)malloc(behaviors->num_values*sizeof(Behavior)); 
  memset(behaviors->values, 0, behaviors->num_values*sizeof(Behavior));
  for(int i = 0; i < (int)p.size(); i++) {
    strcpy(behaviors->values[i].name, p[i].get("name","").asString().c_str());
    strncpy(behaviors->values[i].abbreviation, behaviors->values[i].name, 3);
    behaviors->values[i].abbreviation[3] = '\0';
    behaviors->values[i].color = p[i].get("color",rand() & 0x00FFFFFF).asInt();
    if(!strcasecmp(behaviors->values[i].name, "none"))
      noneInd = i;
  }
  if(noneInd < 0) {
    fprintf(stderr, "Error: you must define a behavior 'none'\n");
    return false;
  } else if(noneInd > 0) {
    // make sure the "none" behavior is at index 0
    Behavior tmp = behaviors->values[0];
    behaviors->values[0] = behaviors->values[noneInd]; 
    behaviors->values[noneInd] = tmp;
    noneInd = 0;
  }
  this->behaviors = behaviors;

  return true;
}

Json::Value SVMBehaviorSequence::SaveBehaviorDefinitions() {
  Json::Value r;
  for(int i = 0; i < behaviors->num_values; i++) {
    Json::Value v(Json::arrayValue);
    v["name"] = behaviors->values[i].name;
    v["color"] = behaviors->values[i].color;
    r[i] = v;
  }
  return r;
}

void SVMBehaviorSequence::FreeBehaviorDefinitions() {
  if(behaviors) {
    free(behaviors->values);
    free(behaviors);
  }
  behaviors = NULL;
}

char *BoutOperatorToString(char *str, BoutOperator op) {
  switch(op) {
  case B_SUM: strcpy(str, "sum"); break;
  case B_AVE: strcpy(str, "ave"); break; 
  case B_VAR: strcpy(str, "var"); break; 
  case B_DEV: strcpy(str, "dev"); break; 
  case B_RAW: strcpy(str, "raw"); break; 
  case B_MIN: strcpy(str, "min"); break; 
  case B_MAX: strcpy(str, "max"); break; 
  case B_DUR: strcpy(str, "dur"); break; 
  }
  return str;
}

BoutOperator BoutOperatorFromString(char *str) {
  if(!strcmp(str, "sum")) return B_SUM;
  else if(!strcmp(str, "ave")) return B_AVE;
  else if(!strcmp(str, "var")) return B_VAR;
  else if(!strcmp(str, "dev")) return B_DEV;
  else if(!strcmp(str, "raw")) return B_RAW;
  else if(!strcmp(str, "min")) return B_MIN;
  else if(!strcmp(str, "max")) return B_MAX;
  else if(!strcmp(str, "dur")) return B_DUR;
  return (BoutOperator)-1;
}

Json::Value SVMBehaviorSequence::SaveBoutFeatureParams() {
   Json::Value fe;
   char op[1000];
   for(int i = 0; i < num_bout_features; i++) {
     Json::Value w, f, reg;
     for(int j = 0; j < bout_features[i].num_regions; j++) {
       Json::Value r;
       r["frame_feature"] = bout_features[i].regions[j].frame_feature;
       r["hist_bin"] = bout_features[i].regions[j].hist_bin;
       r["op"] = BoutOperatorToString(op,bout_features[i].regions[j].op);
       r["t_start"] = bout_features[i].regions[j].t_start;
       r["t_end"] = bout_features[i].regions[j].t_end;
       r["is_global"] = bout_features[i].regions[j].is_global;
       r["frame_coords"] = bout_features[i].regions[j].frame_coords;
       w[j] = bout_features[i].weights[j];
       reg[j] = r;
     }
     f["regions"] = reg;
     f["weights"] = w;
     f["thresh"] = bout_features[i].thresh;
     f["num_thresholds"] = bout_features[i].num_thresholds;
     f["mu"] = bout_features[i].mu;
     f["gamma"] = bout_features[i].gamma;
     f["absolute_value"] = bout_features[i].absolute_value;
     f["time_normalize"] = bout_features[i].time_normalize;
     fe[i] = f;
   }
   return fe;
 }

 void SVMBehaviorSequence::LoadBoutFeatureParams(const Json::Value &fe) {
   FreeBoutFeatureParams();

   num_bout_features = fe.size();
   bout_features = (BoutFeature*)realloc(bout_features, num_bout_features*sizeof(BoutFeature));

   char tmp[1000], op[1000];
   for(int i = 0; i < num_bout_features; i++) {
     Json::Value f = fe[i];
     Json::Value w = f["weights"];
     Json::Value reg = f["regions"];
     bout_features[i].num_regions = reg.isArray() ? reg.size() : 1;
     assert(bout_features[i].num_regions <= MAX_REGIONS);
     for(int j = 0; j < (int)bout_features[i].num_regions; j++) {
       Json::Value r = reg.isArray() ? reg[j] : reg;
       bout_features[i].regions[j].frame_feature = r["frame_feature"].asInt();
       strcpy(op, r["op"].asString().c_str());
       bout_features[i].regions[j].op = BoutOperatorFromString(op);
       bout_features[i].regions[j].t_start = r.get("t_start",0).asDouble();
       bout_features[i].regions[j].t_end = r.get("t_end",1).asDouble();
       bout_features[i].regions[j].hist_bin = r.get("hist_bin",-1).asInt();
       bout_features[i].regions[j].is_global = r.get("is_global",false).asBool();
       bout_features[i].regions[j].frame_coords = r.get("frame_coords",0).asInt();
       bout_features[i].weights[j] = w.isArray() ? w[j].asDouble() : w.asDouble();
     }
     bout_features[i].thresh = f.get("thresh",0).asDouble();
     bout_features[i].num_thresholds = f.get("num_thresholds",0).asInt();
     bout_features[i].mu = f.get("mu",0).asDouble();
     bout_features[i].gamma = f.get("gamma",0).asDouble();
     bout_features[i].absolute_value = f.get("absolute_value",false).asBool();
     bout_features[i].time_normalize = f.get("time_normalize",false).asBool();
     strcpy(tmp, f.get("name","").asString().c_str());
     bout_features[i].name = strlen(tmp) ? StringCopy(tmp) : NULL;
   }
 }

void SVMBehaviorSequence::FreeBoutFeatureParams() {
  if(bout_features) {
    for(int i = 0; i < num_bout_features; i++)
      if(bout_features[i].name)
	free(bout_features[i].name);
    free(bout_features);
  }
  bout_features = NULL;
}


 Json::Value SVMBehaviorSequence::SaveFrameFeatureParams() {
   Json::Value params;
   for(int i = 0; i < num_frame_features; i++) {
     Json::Value c, t;
     if(frame_features[i].name) c["name"] = frame_features[i].name;
     for(int j = 0; j < frame_features[i].num_histogram_bins-1; j++) 
       t[j] = frame_features[i].histogram_thresholds[j];
     c["histogram_thresholds"] = t;
     params[i] = c;
   }
   return params;
 }

void SVMBehaviorSequence::LoadFrameFeatureParams(const Json::Value &params) {
  FreeFrameFeatureParams();
  num_frame_features = params.size();
  frame_features = (FrameFeature*)realloc(frame_features, num_frame_features*sizeof(FrameFeature));
  memset(frame_features, 0, num_frame_features*sizeof(FrameFeature));
  
  char tmp[5000];
  for(int i = 0; i < num_frame_features; i++) {
    Json::Value p = params[i];
    strcpy(tmp, p.get("name","").asString().c_str());
    frame_features[i].name = strlen(tmp) ? StringCopy(tmp) : NULL;
    frame_features[i].num_histogram_bins = p.isMember("histogram_thresholds") ? 
      p["histogram_thresholds"].size()+1 : p.get("num_histogram_bins",0).asInt();
    frame_features[i].histogram_thresholds = (double*)malloc(sizeof(double)*frame_features[i].num_histogram_bins);
    for(int j = 0; j < frame_features[i].num_histogram_bins-1; j++) 
      frame_features[i].histogram_thresholds[j] = p["histogram_thresholds"][j].asDouble();
  }
}

void SVMBehaviorSequence::FreeFrameFeatureParams() {
  if(frame_features) {
    for(int i = 0; i < num_frame_features; i++) {
      if(frame_features[i].name)
	free(frame_features[i].name);
      if(frame_features[i].histogram_thresholds)
	free(frame_features[i].histogram_thresholds);
    }
    free(frame_features);
  }
  frame_features = NULL;
}
     

BehaviorBoutFeatures::BehaviorBoutFeatures() {
  partial_label = NULL;
  memory_buffer = NULL;
  features = NULL;
  frame_times = NULL;
  fvec = NULL;
  fps = 1;
}

BehaviorBoutFeatures::~BehaviorBoutFeatures() {
  if(features) free(features);
  if(memory_buffer) free(memory_buffer);
  if(fvec) delete fvec;
  if(partial_label) delete partial_label;
}

void BehaviorBoutFeatures::Clear() {
  if(memory_buffer) free(memory_buffer);
  memory_buffer = NULL;
}

 
/*
* Allocate space for feature caches used to compute bout-level features efficiently
*/
void BehaviorBoutFeatures::AllocateBuffers(SVMBehaviorSequence *svm, bool full) {
  int i, j;
  int T = num_frames;
  int log_T = (int)(LOG2(T)+1);
  int num_frame_features = svm->NumFrameFeatures();
  FrameFeature *frame_features = svm->GetFrameFeatureDefs();
  int num_bout_features = svm->NumBoutFeatures();
  BoutFeature *bout_features = svm->GetBoutFeatureDefs();
  
  // Compute a safe amount to pad feature buffers, such that all bout features can be computed
  // without checking to see if they go out out bounds
  pad1 = 0;
  pad2 = 0;
  for(i = 0; i < num_bout_features; i++) {
    for(j = 0; j < bout_features[i].num_regions; j++) {
      if(bout_features[i].regions[j].frame_coords == 0) {
	pad1 = my_max(pad1, -bout_features[i].regions[j].t_start*T);
	pad2 = my_max(pad2, (bout_features[i].regions[j].t_end-1)*T);
      } else {
	pad1 = my_max(pad1, -bout_features[i].regions[j].t_start);
	pad2 = my_max(pad2, bout_features[i].regions[j].t_end);
      } 
    }
  } 
  pad1++;
  pad2++;

  if(full) {
    long int cache_features_size = 
      num_frame_features*sizeof(double*) +        // integral_features
      num_frame_features*(T+pad1+pad2+1)*sizeof(double) + // integral_features[i]
      num_frame_features*sizeof(double*) +        // integral_sqr_features
      num_frame_features*(T+pad1+pad2+1)*sizeof(double) + // integral_sqr_features[i]
      num_frame_features*sizeof(double*) +        // smoothed_features
      num_frame_features*(T+pad1+pad2)*sizeof(double) + // smoothed_features[i]
      num_frame_features*sizeof(double*) +        // safe_features
      num_frame_features*(T+pad1+pad2)*sizeof(double) + // safe_features[i]
      num_frame_features*sizeof(int*) +           // histogram_bins
      num_frame_features*sizeof(int)*T +          // histogram_bins[i]
      num_frame_features*sizeof(double**) +       // integral_histogram_features
      2*num_frame_features*(sizeof(double**)+log_T*(sizeof(double*)+(T+pad1+pad2)*sizeof(double)));  // features_min, features_max
    for(i = 0; i < num_frame_features; i++) {
      cache_features_size +=
	sizeof(double*)*frame_features[i].num_histogram_bins + //integral_histogram_features[i]
	frame_features[i].num_histogram_bins*(T+pad1+pad2+1)*sizeof(double); //integral_histogram_features[i][j]
      //num_frame_features*sizeof(double*)*p->num_histogram_bins + //integral_histogram_features[i]
      //num_frame_features*p->num_histogram_bins*(T+1)*3*sizeof(double)+10000; //integral_histogram_features[i][j]
    }
    this->memory_buffer = (unsigned char*)malloc(cache_features_size);
  }

  if(!features) {
    // Extract regular raw features
    //this->features = (double**)malloc(num_frame_features*(sizeof(double*)+(num_frame_features+1)*T*sizeof(double)));
    this->features = (double**)malloc((num_frame_features+1)*(sizeof(double*)+(T+1)*sizeof(double)));
    double *ptr = (double*)(this->features+num_frame_features);
    for(i = 0; i < num_frame_features; i++) {
      this->features[i] = ptr;
      ptr += T;
    }
    this->frame_times = ptr;
  }

  this->partial_label = NULL;

  this->fvec = NULL;
}


/*
* Precompute certain data structures like integral images, such that bout features can
* be computed more efficiently
*/
void BehaviorBoutFeatures::ComputeCaches(SVMBehaviorSequence *svm) {
  //if(!svm->histogram_thresholds) return;
  if(!memory_buffer)
    AllocateBuffers(svm);

  int i, j, k, T = num_frames;
  double f;
  unsigned char *ptr = memory_buffer;
  int window;
  int num_frame_features = svm->NumFrameFeatures();
  FrameFeature *frame_features = svm->GetFrameFeatureDefs();

  // integral_features[i] store the running sum value of the i_th feature, such that the sum
  // value of that feature for an arbitrary interval [s,e) can be computed as 
  //   integral_features[e]-integral_features[s]
  // To allow the caller to index into integral_features[i] without worrying about bounds checking,
  // a buffer of size T+1 is added to the beginning and end, which behaves as if the feature in the
  // 0th timestep extends to frame -(T+1) and the (T-1)th pixel extends (T+1) frames into the future
  integral_features = (double**)ptr;
  ptr += num_frame_features*sizeof(double*);
  for(i = 0; i < num_frame_features; i++) {
    integral_features[i] = ((double*)ptr)+pad1;
    ptr += (T+pad1+pad2+1)*sizeof(double);

    integral_features[i][-pad1] = 0;
    for(j = -pad1; j < 0; j++) 
      integral_features[i][j+1] = integral_features[i][j] + features[i][0];
    for(j = 0; j < T; j++) 
      integral_features[i][j+1] = integral_features[i][j] + features[i][j];
    for(j = T; j < T+pad2; j++) 
      integral_features[i][j+1] = integral_features[i][j] + features[i][T-1];
  }

  integral_sqr_features = (double**)ptr;
  ptr += num_frame_features*sizeof(double*);
  for(i = 0; i < num_frame_features; i++) {
    integral_sqr_features[i] = ((double*)ptr)+pad1;
    ptr += (T+pad1+pad2+1)*sizeof(double);

    integral_sqr_features[i][-pad1] = 0;
    for(j = -pad1; j < 0; j++) 
      integral_sqr_features[i][j+1] = integral_sqr_features[i][j] + SQR(features[i][0]);
    for(j = 0; j < T; j++) 
      integral_sqr_features[i][j+1] = integral_sqr_features[i][j] + SQR(features[i][j]);
    for(j = T; j < T+pad2; j++) 
      integral_sqr_features[i][j+1] = integral_sqr_features[i][j] + SQR(features[i][T-1]);
  }

  // Compute smoothed versions of the raw features
  smoothed_features = (double**)ptr;
  ptr += num_frame_features*sizeof(double*);
  for(i = 0; i < num_frame_features; i++) {
    smoothed_features[i] = ((double*)ptr)+pad1;
    ptr += (T+pad1+pad2)*sizeof(double);
    window = svm->feature_sample_smoothness_window;
    for(j = 0; j < T; j++)
      smoothed_features[i][j] = (integral_features[i][j+window+1] - 
				       integral_features[i][j-window]) / (2*window+1);
    for(j = -pad1; j < 0; j++)
      smoothed_features[i][j] = smoothed_features[i][0];
    for(j = T; j < T+pad2; j++)
      smoothed_features[i][j] = smoothed_features[i][T-1];
  }

  // Compute safe (padded) versions of the raw features
  safe_features = (double**)ptr;
  ptr += num_frame_features*sizeof(double*);
  for(i = 0; i < num_frame_features; i++) {
    safe_features[i] = ((double*)ptr)+pad1;
    ptr += (T+pad1+pad2)*sizeof(double);
    for(j = 0; j < T; j++)
      safe_features[i][j] = features[i][j];
    for(j = -pad1; j < 0; j++) 
      safe_features[i][j] = safe_features[i][0];
    for(j = T; j < T+pad2; j++)
      safe_features[i][j] = safe_features[i][T-1];
  }

  // Compute histogram features.  We discretize the space of values for the i_th feature into num_histogram bins
  // of equal size between min_feature_responses[i] and max_feature_responses[i].  histogram_bins[i][j] stores the
  // index of the assigned bin for the i_th feature at the j_th time step
  histogram_bins = (int**)ptr;
  ptr += num_frame_features*sizeof(int*);
  for(i = 0; i < num_frame_features; i++, ptr += sizeof(int)*T) {
    histogram_bins[i] = (int*)ptr;
    for(j = 0; j < T; j++) {
      f = features[i][j];
      for(k = 0; k < frame_features[i].num_histogram_bins-1; k++) 
	if(f <= frame_features[i].histogram_thresholds[k]) 
	  break;
      histogram_bins[i][j] = k;
    }
  }

  // Compute integral features on the histogram features (e.g. using the same method as for integral_features),
  // such that the total count of the j_th histogram bin for the i_th feature can be computed using
  //   integral_histogram_features[i][j][e]-integral_histogram_features[i][j][s]
  integral_histogram_features = (double***)ptr;
  ptr += num_frame_features*sizeof(double**);
  for(i = 0; i < num_frame_features; i++) {
    integral_histogram_features[i] = (double**)ptr;
    ptr += sizeof(double*)*frame_features[i].num_histogram_bins;
    for(j = 0; j < frame_features[i].num_histogram_bins; j++) {
      integral_histogram_features[i][j] = ((double*)ptr)+pad1;
      ptr += (T+pad1+pad2+1)*sizeof(double);
      integral_histogram_features[i][j][-pad1] = 0;
      for(k = 0; k < T; k++) 
        integral_histogram_features[i][j][k+1] = integral_histogram_features[i][j][k] + (histogram_bins[i][k] == j ? 1 : 0);
      for(k = -pad1; k < 0; k++) 
        integral_histogram_features[i][j][k+1] = integral_histogram_features[i][j][k] + (histogram_bins[i][0] == j ? 1 : 0);
      for(k = T; k < T+pad2; k++) 
        integral_histogram_features[i][j][k+1] = integral_histogram_features[i][j][k] + (histogram_bins[i][T-1] == j ? 1 : 0);
    }
  }

  // Precompute the min and max in time intervals starting at each frame and of length 1,2,4,8,...
  int log_T = (int)(LOG2(T)+1);
  feature_mins = (double***)ptr;
  ptr += num_frame_features*sizeof(double**);
  feature_maxes = (double***)ptr;
  ptr += num_frame_features*sizeof(double**);
  for(i = 0; i < num_frame_features; i++) {
    feature_mins[i] = (double**)ptr;
    ptr += sizeof(double*)*log_T;
    feature_maxes[i] = (double**)ptr;
    ptr += sizeof(double*)*log_T;
    for(j = 0; j < log_T; j++) {
      feature_mins[i][j] = ((double*)ptr)+pad1;
      ptr += (T+pad1+pad2)*sizeof(double);
      feature_maxes[i][j] = ((double*)ptr)+pad1;
      ptr += (T+pad1+pad2)*sizeof(double);
    }
    for(k = -pad1; k < T+pad2; k++) 
      feature_mins[i][0][k] = feature_maxes[i][0][k] = safe_features[i][k];
    for(j = 1; j < log_T; j++) {
      int dt = j ? 1<<(j-1) : 0;
      for(k = -pad1; k < T+pad2-dt; k++) {
        feature_mins[i][j][k] = my_min(feature_mins[i][j-1][k],feature_mins[i][j-1][k+dt]);
        feature_maxes[i][j][k] = my_max(feature_maxes[i][j-1][k],feature_maxes[i][j-1][k+dt]);
      }
      int s = T+pad2-dt-1;
      for(k = T+pad2-dt; k < T+pad2; k++) {
        feature_mins[i][j][k] = feature_mins[i][j][s];
        feature_maxes[i][j][k] = feature_maxes[i][j][s];
      }
    }
  }
}

bool BehaviorBoutFeatures::load(const Json::Value &r, StructuredSVM *s) {
  strcpy(fname, r.get("fname", "").asString().c_str());
  return load(fname, (SVMBehaviorSequence*)s);
}

Json::Value BehaviorBoutFeatures::save(StructuredSVM *s) {
  Json::Value r;
  r["fname"] = fname;
  return r;
}




BehaviorBoutSequence::BehaviorBoutSequence(BehaviorBoutFeatures *x, SVMBehaviorSequence *svm) : StructuredLabel(x) {
  features = x;
  filled = false;
  num_bouts = 0;
  bouts = NULL;
  score = loss = slack = 0;
  disable_checks = false;
}

BehaviorBoutSequence::~BehaviorBoutSequence() {
  if(bouts) free(bouts);
}


/*
Json::Value BehaviorBoutSequence::save(StructuredSVM *s) {
  Json::Value r;
  r["fname"] = fname;
  save(fname);
  return r;
}

bool BehaviorBoutSequence::load(const Json::Value &r, StructuredSVM *s) {
  strcpy(fname, r.get("fname", "").asString().c_str());
  bool ret = load(fname);
  remove(fname);
  return ret;
  //return load(fname);
}
*/


Json::Value BehaviorBoutSequence::save(StructuredSVM *s) {
  Json::Value r;
  //r["num_bouts"] = num_bouts;

  if(bouts) {
    Json::Value b(Json::arrayValue);
    for(int j = 0; j < num_bouts; j++) {
      Json::Value c;
      c["start_frame"] = bouts[j].start_frame;
      c["end_frame"] = bouts[j].end_frame;
      c["behavior"] = bouts[j].behavior;
      b[j] = c;
    }
    r["bouts"] = b;
  }

  if(s->IsTesting()) {
    r["fname"] = fname;
    save(fname);
  }
  return r;
}


bool BehaviorBoutSequence::load(const Json::Value &r, StructuredSVM *s) {
  if(r.isMember("fname")) {
    char fname[1000];
    strcpy(fname, r.get("fname", "").asString().c_str());
    bool ret = load(fname);
    //remove(fname);
    return ret;
  } else {
    ((SVMBehaviorSequence*)s)->init_bout_label(this, NULL);
    num_bouts = r["bouts"].size();

    ((SVMBehaviorSequence*)s)->init_bout_label(this, NULL);
    this->bouts = (BehaviorBout*)malloc(sizeof(BehaviorBout)*num_bouts);
    for(int j = 0; j < num_bouts; j++) {
      Json::Value c = r["bouts"][j];
      bouts[j].start_frame = c["start_frame"].asInt();
      bouts[j].end_frame = c["end_frame"].asInt();
      bouts[j].behavior = c["behavior"].asInt();
    }
  }
  return true;
}

void SVMBehaviorSequence::DebugFeatures(const char *fname, BehaviorBoutSequence *y, double *ww) {
   BehaviorBoutFeatures *b = y->features;

  double *tmp_features = (double*)malloc((num_bout_features+3)*sizeof(double));
  double *class_features = (double*)malloc(sizeof(double)*(num_bout_features+behaviors->num_values+2));
  double *class_transitions = class_features + num_bout_features;
  double *unary = class_transitions + behaviors->num_values;
  double *duration = unary+1;
  int *inds = (int*)malloc(sizeof(int)*(num_bout_features+behaviors->num_values+2));
  g_num_features = num_bout_features;
  g_class_features = class_features;
  g_class_transitions = class_transitions;
  double **class_weights, **transition_weights, *unary_weights, *duration_weights;
  init_weight_pointers(ww, class_weights, transition_weights, unary_weights, duration_weights);
  FILE *fout = fopen(fname, "w");

  y->score = 0;
  y->loss = 0;
  for(int i = 0; i < y->num_bouts; i++) {
    int t_p = y->bouts[i].start_frame, t = y->bouts[i].end_frame, 
      c_prev = y->bouts[i].behavior, 
      c_next = i<y->num_bouts-1 ? y->bouts[i+1].behavior : -1;
    psi_bout(b, t_p, t, c_prev, tmp_features, true, false);  // compute bout features

    y->bouts[i].bout_score = 0;
    for(int k = 0; k < num_bout_features; k++) {
      class_features[k] = class_weights[c_prev][k]*tmp_features[k];  // update bout score
      y->bouts[i].bout_score += class_features[k];
    }
    for(int j = 0; j < behaviors->num_values; j++) 
      class_transitions[j] = 0; 
    class_transitions[c_next] = transition_weights[c_prev][c_next];
		      
    y->bouts[i].transition_score = c_next >=0 ? transition_weights[c_prev][c_next] : 0;   
    *unary = y->bouts[i].unary_score = unary_weights[c_prev];
    *duration = y->bouts[i].duration_score = get_duration_score(t-t_p, c_prev, duration_weights);
    double score = y->bouts[i].bout_score + y->bouts[i].transition_score + y->bouts[i].unary_score + y->bouts[i].duration_score;
    y->score += score;
    y->loss += y->bouts[i].loss_fn + y->bouts[i].loss_fp;
		    
    for(int j = 0; j < num_bout_features+behaviors->num_values+2; j++)
      inds[j] = j;
    qsort(inds, num_bout_features+behaviors->num_values+2, sizeof(int), cmp_feature_inds);
    fprintf(fout, "<br><br><a name=\"%d\">Behavior=%s score=%f t=(%d-%d)</a>: bout_score=%f transition_score=%f unary_score=%f duration_score=%f loss_fn=%f loss_fp=%f\n", i, behaviors->values[c_prev].name, (float)score, t_p, t,
	    (float)y->bouts[i].bout_score, (float)y->bouts[i].transition_score, (float)y->bouts[i].unary_score, (float)y->bouts[i].duration_score, (float)y->bouts[i].loss_fn,  (float)y->bouts[i].loss_fp);
    for(int j = 0; j < num_bout_features+behaviors->num_values; j++) {
      if(inds[j] < num_bout_features)
	fprintf(fout, "<br>%lf=%f*%f, %s %s\n", class_features[inds[j]], class_weights[c_prev][inds[j]], tmp_features[inds[j]], behaviors->values[c_prev].name, bout_features[inds[j]].name);
      else if(inds[j] < num_bout_features+behaviors->num_values)
	fprintf(fout, "<br>%lf, %s->%s\n", class_transitions[inds[j]-num_bout_features], behaviors->values[c_prev].name, behaviors->values[inds[j]-num_bout_features].name);
      else if(inds[j] < num_bout_features+behaviors->num_values+1)
	fprintf(fout, "<br>%lf, %s unary\n", *unary, behaviors->values[c_prev].name);
      else 
	fprintf(fout, "<br>%lf, %s duration\n", *duration, behaviors->values[c_prev].name);
    }
  }
  free(class_weights);
  free(class_features);
  free(tmp_features);
  free(inds);
  fclose(fout);
}

#ifdef HAVE_OPENCV
#include "cv.h"
#include "highgui.h"
#endif

#define LABEL_BOUTS 0
 void BehaviorBoutSequence::Visualize(Behaviors *group, const char *fname, char *html) { 
   if(!this->num_bouts) 
     return;

   int h = LABEL_BOUTS ? 50 : 10;
   int T = this->bouts[this->num_bouts-1].end_frame;
   int i;

#ifdef HAVE_OPENCV
   CvFont font;
   cvInitFont(&font, CV_FONT_VECTOR0, 0.5f, 0.4f, 0, 2);

   IplImage *img = cvCreateImage(cvSize(T, h), IPL_DEPTH_8U, 3);
   cvZero(img);

   // Draw bouts as colored rectangles
   for(i = 0; i < this->num_bouts; i++) {
     int c = this->bouts[i].behavior;
     cvRectangle(img, cvPoint(this->bouts[i].start_frame,0), cvPoint(this->bouts[i].end_frame,h), 
		 CV_RGB((group->values[c].color & 0xff0000)>>16,
			(group->values[c].color & 0x00ff00)>>8,
			(group->values[c].color & 0xff)), CV_FILLED);
   }

#if LABEL_BOUTS
   // Draw labels for bouts
   int prev_prev_max_x = -100000, prev_max_x = -100000;
   int last_pos = 2, last_last_pos = 1;
   int y[3], ymin;
   for(i = 0; i < this->num_bouts; i++) {
     int c = this->bouts[i].behavior;
     int m = (this->bouts[i].start_frame+this->bouts[i].end_frame)/2;

     cvGetTextSize(group->values[c].abbreviation, &font, &sz, &ymin);

     if(i == 0) {
       y[0] = my_min(h-2, (h+sz.height)/2);
       y[1] = h-6;
       y[2] = sz.height;
     }

     int pos = 0;
     if(m-sz.width/2 < prev_prev_max_x+4) 
       pos = 3 - last_last_pos - last_pos;
     else if(m-sz.width/2 < prev_max_x+4) 
       pos = last_pos > 0 ? 0 : 1;
     else {
       pos = 0;
       last_pos = 2;
     }
     last_last_pos = last_pos;
     last_pos = pos;
     prev_prev_max_x = prev_max_x;
     prev_max_x = m+sz.width/2;
     cvPutText(img, group->values[c].abbreviation, cvPoint(m-sz.width/2, y[pos]), 
	       &font, CV_RGB(255,255,255));
   }
#endif
#endif

   if(fname) {
     char *html_tmp=(char*)malloc(10000000), folder[1000], fname2[1000], fname3[1000];
     if(!strstr(fname, ".png")) sprintf(fname2, "%s.png", fname);
     else strcpy(fname2, fname);

#ifdef HAVE_OPENCV
     cvSaveImage(fname2, img);
     cvReleaseImage(&img);
#endif

     if(html) {
       ExtractFolderAndFileName(fname, folder, fname2);    

		  

       sprintf(html_tmp, "<img src=\"%s.png\" usemap=\"#%s\" height=50 />\n<map name=\"%s\">\n", fname2, fname2, fname2);

       char str[10000], alt[1000];
       char *ptr = html_tmp+strlen(html_tmp);
       for(i = 0; i < this->num_bouts; i++) {
	 int c = this->bouts[i].behavior;
	 float z = 50.0f/h;
	 double score = this->bouts[i].bout_score + this->bouts[i].transition_score + this->bouts[i].unary_score;
#if USE_DURATION_COST > 0
	 score += this->bouts[i].duration_score;
#endif
	 sprintf(alt, "behavior=%s score=%f bout_score=%f transition_score=%f unary_score=%f duration_score=%f loss_fn=%f loss_fp=%f", group->values[c].name, (float)score,
		 (float)this->bouts[i].bout_score, (float)this->bouts[i].transition_score, (float)this->bouts[i].unary_score, (float)this->bouts[i].duration_score, (float)this->bouts[i].loss_fn,  (float)this->bouts[i].loss_fp);
	 sprintf(str, "  <area shape=\"rect\" coords=\"%d,%d,%d,%d\" href=\"%s_features.html#%d\" title=\"%s\" />\n", (int)(z*this->bouts[i].start_frame), 0, 
		 (int)(z*this->bouts[i].end_frame), (int)(z*h), fname2, i, alt);


	 strcpy(ptr, str);
	 ptr += strlen(str);
       }
       strcpy(ptr, "</map>");
       strcpy(html, html_tmp);
       free(html_tmp);
     }
   }
}
