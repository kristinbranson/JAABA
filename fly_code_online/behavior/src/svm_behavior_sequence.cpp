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

#if DEBUG > 0
const char *bout_feature_names[] = {"feature_sample_smoothness_window", "num_temporal_levels", "num_bout_max_thresholds", 
			 "num_bout_min_thresholds", "num_bout_change_points", "num_histogram_bins", 
		     "num_histogram_temporal_levels", "num_difference_temporal_levels", "num_harmonic_features",
		     "use_bout_sum_features", "use_bout_ave_features", "use_bout_sum_absolute_features", "use_bout_ave_absolute_features", "use_standard_deviation", "use_sum_variance",
		     "use_bout_max_feature", "use_bout_min_feature", 
		     "use_global_difference_max_ave_features", "use_global_difference_min_ave_features", "use_global_difference_ave_ave_features", 
		     "use_global_difference_max_sum_features", "use_global_difference_min_sum_features", "use_global_difference_ave_sum_features",
		     "use_bout_change", "use_bout_absolute_change", "use_histogram_sum_features",
		     "use_histogram_ave_features", "use_sum_harmonic_features", "use_ave_harmonic_features", "use_sum_absolute_harmonic_features",
		     "use_ave_absolute_harmonic_features", "use_start_sum_absolute_diff_haar_features",
		     "use_end_sum_absolute_diff_haar_features", "use_start_sum_diff_haar_features", "use_end_sum_diff_haar_features",
		     "use_start_ave_absolute_diff_haar_features", "use_end_ave_absolute_diff_haar_features", "use_start_ave_diff_haar_features",
			 "use_end_ave_diff_haar_features"};
char *g_feature_names[MAX_FEATURES];
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




SVMBehaviorSequence::SVMBehaviorSequence(int num_feat, struct _BehaviorGroups *behaviors, int beh, SVMFeatureParams *sparams) : StructuredSVM() {
	for(int i = 0; i < behaviors->num; i++) {
		if(beh == -1 || i == beh) {
			restrict_behavior_features[i] = (bool**)malloc(behaviors->behaviors[i].num_values*sizeof(bool*));
			memset(restrict_behavior_features[i], 0, sizeof(bool*)*behaviors->behaviors[i].num_values);
			/*
			CSC & SB 20110324: changed beh to i
			*/
		}
	}

	Init(num_feat, behaviors, beh, sparams);
}

SVMBehaviorSequence::SVMBehaviorSequence(struct _BehaviorGroups *behaviors, int beh)  : StructuredSVM(){
	Init(0, behaviors, beh);
	for(int i = 0; i < behaviors->num; i++) {
		if(beh == -1 || i == beh) {
			restrict_behavior_features[i] = (bool**)malloc(behaviors->behaviors[i].num_values*sizeof(bool*));
			memset(restrict_behavior_features[i], 0, sizeof(bool*)*behaviors->behaviors[i].num_values);
			/*
			CSC & SB 20110324: changed beh to i
			restrict_behavior_features[beh] = (bool**)malloc(behaviors->behaviors[beh].num_values*sizeof(bool*));
			memset(restrict_behavior_features[beh], 0, sizeof(bool*)*behaviors->behaviors[beh].num_values);
			*/
		}
	}
}

void SVMBehaviorSequence::Init(int num_feat, struct _BehaviorGroups *behaviors, int beh, SVMFeatureParams *sparams) {
        eps = .015;
	C = 10.0;
	method = SPO_DUAL_MULTI_SAMPLE_UPDATE_WITH_CACHE;// SPO_DUAL_UPDATE_WITH_CACHE;
	debugLevel = 3;

	max_inference_learning_frames = -1;  // don't disregard any frames at train time to speedup inference
	//max_inference_learning_frames = 200;

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
	updateFromCacheThread = false;//true;

	strcpy(debugdir, "");
	debug_predictions = debug_weights = debug_features = debug_model = false;

	this->behaviors = behaviors;
	this->behavior = beh;
	this->feature_diff_frames = 10;

	class_training_transitions = NULL;
	class_training_transitions_count = class_training_count = NULL;
	features_mu = features_gamma = NULL;
	histogram_thresholds = min_thresholds = max_thresholds = NULL;
	min_bout_duration = 1;
	feature_names = NULL;

	int i;
	num_base_features = num_feat;
	for(i = 0; i < num_base_features; i++) {
		feature_params[i] = sparams ? sparams[i] : DefaultParams();
		assert(feature_params[i].num_temporal_levels <= MAX_TEMPORAL_LEVELS && 
			feature_params[i].num_histogram_temporal_levels <= MAX_TEMPORAL_LEVELS &&
			feature_params[i].num_difference_temporal_levels <= MAX_TEMPORAL_LEVELS &&
			feature_params[i].num_harmonic_features <= MAX_HARMONIC_LEVELS);
	}


	
	sizePsi = 0;
	compute_feature_space_size();
	for(int bh = 0; bh < behaviors->num; bh++) {
		num_classes[bh] = behaviors->behaviors[bh].num_values;

		false_negative_cost[bh] = (double*)malloc(2*sizeof(double)*num_classes[bh]);
		false_positive_cost[bh] = false_negative_cost[bh]+num_classes[bh];
#if USE_DURATION_COST > 0
		min_frame_duration[bh]  = (double*)malloc(2*sizeof(double)*num_classes[bh]);
		max_frame_duration[bh]  = min_frame_duration[bh]+num_classes[bh];
#endif
		// It is assumed here that the label 0 is the unknown label.  The user could define custom
		// class-specific values for these constants .
		// In the future, we could implement this with a per class-pair confusion cost
#if SCALED_LOSS > 0 
		false_negative_cost[bh][0] = 0;  // 0 when beh 0 represents "other"
		false_positive_cost[bh][0] = 0;  // 0 when beh 0 represents "other"
		for(i = 1; i < num_classes[bh]; i++) { 
			false_negative_cost[bh][i] = 10; // 100 for eric's videos
			false_positive_cost[bh][i] = 10; // 100 for eric's videos
		}
#else
		false_negative_cost[bh][0] = 0;
		false_positive_cost[bh][0] = 0;  //
		for(i = 1; i < num_classes[bh]; i++) {  
			false_negative_cost[bh][i] = 1; 
			false_positive_cost[bh][i] = 1;
		}
#endif
		//for(i = 0; i < num_classes[bh]; i++)
		//	num_bouts_per_class[bh][i] = 0;

#if USE_DURATION_COST > 0
		for(i = 0; i < num_classes[bh]; i++) {
			min_frame_duration[bh][i] = INFINITY;
			max_frame_duration[bh][i] = -INFINITY;
		}
#endif
		sizePsi += getPsiSize(this->num_features, this->num_classes[bh]);
	}
}

SVMBehaviorSequence::~SVMBehaviorSequence() {
	int i;

	if(feature_names) {
		for(i = 0; i < num_features; i++)
			if(feature_names[i])
				free(feature_names[i]);
		free(feature_names);
	}
	if(class_training_transitions) free(class_training_transitions);
	if(class_training_transitions_count) free(class_training_transitions_count);
	if(class_training_count) free(class_training_count);
}




/*
* Default parameters for expanding per-frame features into bout-level features
*/
SVMFeatureParams SVMBehaviorSequence::DefaultParams() {
	SVMFeatureParams p;
	memset(&p, 0, sizeof(SVMFeatureParams));

	// When computing the min and max feature responses, average over a window of 2*5+1 frames
	p.feature_sample_smoothness_window = 5;

	// Compute the average feature response over the entire bout, the 1st half of the bout, 
	// and the 2nd half of the bout.  Add features for both the bout average and raw sum response
	p.num_temporal_levels = 2;   
	p.use_bout_sum_features = p.use_bout_ave_features = true;
	p.use_bout_sum_absolute_features = p.use_bout_ave_absolute_features = false;

	// Use the standard deviation of the feature over the entire bout
	p.use_standard_deviation = true;
	p.use_sum_variance = false;

	// Don't use the bout-wise minimum and maximum or thresholded versions of them
	p.use_bout_max_feature = p.use_bout_min_feature = false;
	p.num_bout_max_thresholds = p.num_bout_min_thresholds = 0;

	// Compare the bout average features to the global minimum, maximum, and average
	p.use_global_difference_max_ave_features = p.use_global_difference_min_ave_features = p.use_global_difference_ave_ave_features = false;
	p.use_global_difference_max_sum_features = p.use_global_difference_min_sum_features = p.use_global_difference_ave_sum_features = true;

	// Don't use the difference at the start of the bout and end of the bout as a feature
	p.num_bout_change_points = 0;
	p.use_bout_change = p.use_bout_absolute_change = false;

	// Quantize the feature space into 6 bins, then histogram over the duration of the bout, computing both the sum
	// and average histogram response
	p.num_histogram_bins = 6;  
	p.num_histogram_temporal_levels = 1;
	p.use_histogram_sum_features = p.use_histogram_ave_features = true;

	// Compute haar-like filter responses comparing the average feature response of this bout to
	// frames preceding (and following) this bout of size 1/2, and 1/4 the size of this bout 
	p.num_difference_temporal_levels = 2;
	p.use_start_ave_absolute_diff_haar_features = p.use_end_ave_absolute_diff_haar_features = 
		p.use_start_ave_diff_haar_features = p.use_end_ave_diff_haar_features = false;
	p.use_start_sum_absolute_diff_haar_features = p.use_end_sum_absolute_diff_haar_features = 
		p.use_start_sum_diff_haar_features = p.use_end_sum_diff_haar_features = true;

	// Don't use harmonic features in this bout
	p.num_harmonic_features = 0;
	p.use_sum_harmonic_features = true;
	p.use_ave_harmonic_features = p.use_sum_absolute_harmonic_features = p.use_ave_absolute_harmonic_features = false;


	return p;
}


int cmp_double(const void * a, const void * b) { double d = *(double*)a - *(double*)b;  return (d < 0 ? -1 : (d > 0 ? 1 : 0)); }

void SVMBehaviorSequence::saveBoutFeatures(StructuredDataset *dataset, const char *filename, bool sphered, bool addRandBouts, int window_size) {
	// EYRUN: Store normalized bout features of training set in a txt file (want to verify that bout features are good enough to separate the classes)
	// Count the total number of bouts in the training set
	int num_bouts = 0;
	int num_rand_factor = 5;
	//double window_size = 4; // 8, 16, 32, 64, 128
	StructuredExample **ex = dataset->examples;
	BehaviorBoutSequence *y;
	BehaviorBoutFeatures *x;
	double *tmp_features = (double*)malloc(sizeof(double)*num_features);

	int num_bouts_extra = 0;
	for(int n = 0; n < dataset->num_examples; n++) {
		for(int beh = 0; beh < behaviors->num; beh++) {
			if(behavior < 0 || beh == behavior) {
				num_bouts += ((BehaviorBoutSequence*)ex[n]->y)->num_bouts[beh];
			}
		}
		num_bouts_extra += ((BehaviorBoutFeatures*)ex[n]->x)->num_frames - window_size;
	}

	if(addRandBouts)
		num_bouts = num_bouts + num_bouts_extra; //num_bouts = num_bouts * (num_rand_factor+1); //

	FILE *featureFile, *featureFileUnsphered;
	featureFile = fopen(filename,"w");
	fprintf(featureFile, "num_features = %d\nnum_bouts = %d", num_features, num_bouts);
	for(int n = 0; n < dataset->num_examples; n++) {
		y = ((BehaviorBoutSequence*)ex[n]->y);
		x = ((BehaviorBoutFeatures*)ex[n]->x);
		int frame_labels[x->num_frames];
		if(!x->memory_buffer)
		  x->ComputeCaches(this);
		for(int beh = 0; beh < behaviors->num; beh++) {
			if(behavior < 0 || beh == behavior) {
				for(int j = 0; j < y->num_bouts[beh]; j++) {
					psi_bout(x, y->bouts[beh][j].start_frame, y->bouts[beh][j].end_frame, beh, y->bouts[beh][j].behavior, tmp_features, sphered, false);
					fprintf(featureFile, "\n%d %d %d %d ", n, y->bouts[beh][j].start_frame, y->bouts[beh][j].end_frame, y->bouts[beh][j].behavior);
					for(int i = 0; i < num_features; i++) {
						fprintf(featureFile, "%f ", tmp_features[i]);
					}
					for(int i = y->bouts[beh][j].start_frame; i < y->bouts[beh][j].end_frame; i++) 
						frame_labels[i] = y->bouts[beh][j].behavior;
				}
			}
		}
		if(addRandBouts) {
			// Generate random intervals and calculate the bout-features for those
			//int num_rand = y->num_bouts[0] * num_rand_factor;
			int num_rand = x->num_frames - window_size;
			for(int r = window_size/2; r<x->num_frames-window_size/2; r++) {
				int f_start, f_dur, f_end;
				// generate random start frame between 0 and x->num_frames - 5
			//	f_start = rand() % (x->num_frames - 5);
				// generate random duration between 4 and 50 (the number of frames we go back)
			//	f_dur = rand() % 47 + 4;
			//	f_end = my_min(f_start + f_dur - 1, x->num_frames - 1);

				f_start = r-window_size/2;
				f_end = f_start+window_size;

				// calculate bout score for that interval
				psi_bout(x, f_start, f_end, 0, -1, tmp_features, sphered, false);
				fprintf(featureFile, "\n%d %d %d %d ", n, f_start, f_end, frame_labels[r]);
				for(int i = 0; i < num_features; i++) {
					fprintf(featureFile, "%f ", tmp_features[i]);
				}
			}	
		}
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
	int i, j, k, n, ind, beh, num_bouts = 0, curr = 0;
	double **feat, *ptr, mean, num_histogram_bins = 0, num_max_thresholds = 0, num_min_thresholds = 0, w, target, m;
	BehaviorBoutSequence *y;
	double *tmp_features = (double*)malloc(sizeof(double)*num_features);
	SVMFeatureParams *p;
	BehaviorBoutFeatures *x;
	StructuredExample **ex = dataset->examples;

      
	// Count the total number of bouts in the training set
	for(n = 0; n < dataset->num_examples; n++) {
		for(beh = 0; beh < behaviors->num; beh++) {
			if(behavior < 0 || beh == behavior) {
				y = ((BehaviorBoutSequence*)ex[n]->y);
				num_bouts += y->num_bouts[beh];
				//for(i=0; i<y->num_bouts[beh]; i++)
				//	num_bouts_per_class[beh][y->bouts[beh][i].behavior]++;
			}
		}
	}
	

#if SCALED_LOSS > 0
	for (beh = 0; beh < behaviors->num; beh++) {
		for(i=0; i<num_classes[beh]; i++){
		  //false_negative_cost[beh][i] /= num_bouts_per_class[beh][i];
		  //false_positive_cost[beh][i] /= num_bouts_per_class[beh][i];
			false_negative_cost[beh][i] /= class_training_count[beh][i];
			false_positive_cost[beh][i] /= class_training_count[beh][i];
		}
	}
#endif

	feat = (double**)malloc(sizeof(double*)*num_base_features+num_bouts*num_base_features*sizeof(double));
	ptr = (double*)(feat+num_base_features);

	// Compute the bout-average feature response for every bout in the training set, and store them into feat[][]
	for(i = 0; i < num_base_features; i++, ptr += num_bouts) {
		p = &feature_params[i];
		feat[i] = ptr;
		num_histogram_bins += p->num_histogram_bins;
		num_max_thresholds += p->num_bout_max_thresholds;
		num_min_thresholds += p->num_bout_min_thresholds;
		curr = 0;
		for(n = 0; n < dataset->num_examples; n++) {
			y = ((BehaviorBoutSequence*)ex[n]->y);
			x = ((BehaviorBoutFeatures*)ex[n]->x);
			for(beh = 0; beh < behaviors->num; beh++) {
				if(behavior < 0 || beh == behavior) {
					for(j = 0; j < y->num_bouts[beh]; j++) {
						mean = 0;
						for(k = y->bouts[beh][j].start_frame; k < y->bouts[beh][j].end_frame; k++) {
							assert(k >= 0);  // CSC 20110324
							mean += x->features[i][k];
						}
						if (y->bouts[beh][j].end_frame == y->bouts[beh][j].start_frame)
							printf("Warning: bout has 0 frames!");
						mean /= (y->bouts[beh][j].end_frame-y->bouts[beh][j].start_frame);
						feat[i][curr++] = mean;
					}
				}
			}
		}
		assert(curr == num_bouts);
	}

	// Compute thresholds for assigning histogram bins.  The thresholds are chosen such that an equal number
	// of training examples (training example bouts) would lie in each histogram bin (with respect to the mean of each bout)
	histogram_thresholds = (double**)realloc(histogram_thresholds, 3*sizeof(double*)*num_base_features + 
		(num_histogram_bins+num_max_thresholds+num_min_thresholds)*sizeof(double));
	min_thresholds = histogram_thresholds + num_base_features;
	max_thresholds = min_thresholds + num_base_features;
	ptr = (double*)(max_thresholds + num_base_features);
	for(i = 0; i < num_base_features; i++) {
		p = &feature_params[i];
		qsort(feat[i], num_bouts, sizeof(double), cmp_double);

		histogram_thresholds[i] = p->num_histogram_bins ? ptr : NULL; 
		ptr += p->num_histogram_bins;
		min_thresholds[i] = p->num_bout_min_thresholds ? ptr : NULL; 
		ptr += p->num_bout_min_thresholds;
		max_thresholds[i] = p->num_bout_max_thresholds ? ptr : NULL; 
		ptr += p->num_bout_max_thresholds;

		for(j = 0; j < p->num_histogram_bins; j++) {
			target = (j+1)*num_bouts / (double)p->num_histogram_bins;
			ind = (int)target;
			w = 1-(target-ind);
			while(j && ind < num_bouts && feat[i][ind] == histogram_thresholds[i][j-1]) ind++;
			histogram_thresholds[i][j] = feat[i][my_min(ind,num_bouts-1)]*w + feat[i][my_min(ind+1,num_bouts-1)]*(1-w);
			assert(!isnan(histogram_thresholds[i][j]));
		}
	}

	for(i = 0; i < num_base_features; i++) {
		// Compute the bout-level min over all training examples
		curr = 0;
		for(n = 0; n < dataset->num_examples; n++) {
			y = ((BehaviorBoutSequence*)ex[n]->y);
			x = ((BehaviorBoutFeatures*)ex[n]->x);
			for(beh = 0; beh < behaviors->num; beh++) {
				if(behavior < 0 || beh == behavior) {
					for(j = 0; j < y->num_bouts[beh]; j++) {
						m = 1000000;
						for(k = y->bouts[beh][j].start_frame; k < y->bouts[beh][j].end_frame; k++) {
							m = my_min(x->features[i][k],m);
						}
						feat[i][curr++] = m;
					}
				}
			}
		}


		p = &feature_params[i];
		qsort(feat[i], num_bouts, sizeof(double), cmp_double);

		for(j = 0; j < p->num_bout_min_thresholds; j++) {
			target = (j+1)*num_bouts / (double)(p->num_bout_min_thresholds+1);
			ind = (int)target;
			w = 1-(target-ind);
			while(j && ind < num_bouts && feat[i][ind] == min_thresholds[i][j-1]) ind++;
			min_thresholds[i][j] = feat[i][my_min(ind,num_bouts-1)]*w + feat[i][my_min(ind+1,num_bouts-1)]*(1-w);
		}
	}

	for(i = 0; i < num_base_features; i++) {
		p = &feature_params[i];

		// Compute the bout-level max over all training examples
		curr = 0;
		for(n = 0; n < dataset->num_examples; n++) {
			y = ((BehaviorBoutSequence*)ex[n]->y);
			x = ((BehaviorBoutFeatures*)ex[n]->x);
			for(beh = 0; beh < behaviors->num; beh++) {
				if(behavior < 0 || beh == behavior) {
					for(j = 0; j < y->num_bouts[beh]; j++) {
						m = -100000000;
						for(k = y->bouts[beh][j].start_frame; k < y->bouts[beh][j].end_frame; k++) {
							m = my_max(x->features[i][k],m);
						}
						feat[i][curr++] = m;
					}
				}
			}
		}

		qsort(feat[i], num_bouts, sizeof(double), cmp_double);

		for(j = 0; j < p->num_bout_max_thresholds; j++) {
			target = (j+1)*num_bouts / (double)(p->num_bout_max_thresholds+1);
			ind = (int)target;
			w = 1-(target-ind);
			while(j && ind < num_bouts && feat[i][ind] == max_thresholds[i][j-1]) ind++;
			max_thresholds[i][j] = feat[i][my_min(ind,num_bouts-1)]*w + feat[i][my_min(ind+1,num_bouts-1)]*(1-w);
		}
	}

	feature_names = (char**)realloc(feature_names, sizeof(char*)*num_features);
	compute_feature_space_size();

	features_mu = (double*)realloc(features_mu, num_features*sizeof(double)*2);
	features_gamma = features_mu + num_features;
	for(i = 0; i < num_features; i++) 
		features_mu[i] = features_gamma[i] = 0;
	for(n = 0; n < dataset->num_examples; n++) {
		y = ((BehaviorBoutSequence*)ex[n]->y);
		x = ((BehaviorBoutFeatures*)ex[n]->x);

		// Now that histogram thresholds are set, we can compute bout-wise features for the training set
		if(!x->memory_buffer)
		  x->ComputeCaches(this);
		for(beh = 0; beh < behaviors->num; beh++) {
			if(behavior < 0 || beh == behavior) {
				for(j = 0; j < y->num_bouts[beh]; j++) {
					psi_bout(x, y->bouts[beh][j].start_frame, y->bouts[beh][j].end_frame, beh, y->bouts[beh][j].behavior, tmp_features, false, false);
					for(i = 0; i < num_features; i++) {
					  // Now that we have bout features, we can compute the mean of each bout-wise feature over the training set
						assert(!isnan(tmp_features[i]));
						features_mu[i] += tmp_features[i];
						assert(!isnan(features_mu[i]));
					}
				}
			}
		}
#if KEEP_FEATURES_IN_MEMORY == 0
		x->Clear();
#endif
	}
	for(i = 0; i < num_features; i++) 
		features_mu[i] /= num_bouts;

	// Compute the (inverse of the) standard deviation bout-average feature response
	for(n = 0; n < dataset->num_examples; n++) {
		y = ((BehaviorBoutSequence*)ex[n]->y);
		x = ((BehaviorBoutFeatures*)ex[n]->x);
		if(!x->memory_buffer)
		  x->ComputeCaches(this);
		for(beh = 0; beh < behaviors->num; beh++) {
			if(behavior < 0 || beh == behavior) {
				for(j = 0; j < y->num_bouts[beh]; j++) {
					psi_bout(x, y->bouts[beh][j].start_frame, y->bouts[beh][j].end_frame, beh, -1, tmp_features, false, false);
					for(i = 0; i < num_features; i++) {
						features_gamma[i] += SQR(tmp_features[i]-features_mu[i]);
						assert(!isnan(features_gamma[i]));
					}
				}
			}
		}
#if KEEP_FEATURES_IN_MEMORY == 0
		x->Clear();
#endif
	}
	for(i = 0; i < num_features; i++) {
		if(features_gamma[i] == 0)
			fprintf(stderr, "WARNING: feature %s seems to be the same for every example in the training set\n", feature_names[i]);
		features_gamma[i] = features_gamma[i] ? 1.0/sqrt(features_gamma[i]/num_bouts) : 0;
	}

	saveBoutFeatures(dataset,"traindata_5.txt",true,false,0);
	saveBoutFeatures(dataset,"traindata_5_unsphered.txt",false,false,0);
/*
	saveBoutFeatures(dataset, "train_bout_feat_4.txt",   false, true, 4);
	saveBoutFeatures(dataset, "train_bout_feat_8.txt",   false, true, 8);
	saveBoutFeatures(dataset, "train_bout_feat_16.txt",  false, true, 16);
	saveBoutFeatures(dataset, "train_bout_feat_32.txt",  false, true, 32);
	saveBoutFeatures(dataset, "train_bout_feat_64.txt",  false, true, 64);
	saveBoutFeatures(dataset, "train_bout_feat_128.txt", false, true, 128);
	//saveBoutFeatures(dataset, "train_bout_feat_unsphered.txt", false, true);
*/
	free(feat);
	free(tmp_features);
}

#if DEBUG > 1
#define MAX_FEATURES 10000
int g_feature_map[MAX_FEATURES] = {0}; // if first element is initialized the rest will be zero-filled
int g_feature_map2[MAX_FEATURES] = {0}; // if first element is initialized the rest will be zero-filled
int g_feature_map3[MAX_FEATURES] = {0}; // if first element is initialized the rest will be zero-filled
int g_csc_ind;

#define ADD_FEATURE(feat, ind, val, feature_id, feature2_id, feature3_id) {\
	g_csc_ind = ind; \
	(feat)[(ind)] = (val); \
	ind++; \
	assert(abs(f) < 100000000); \
	g_feature_map[(g_csc_ind)] = feature_id; \
	g_feature_map2[(g_csc_ind)] = feature2_id;\
	g_feature_map3[(g_csc_ind)] = feature3_id; \
}
#else
#define ADD_FEATURE(feat, ind, val, feature_id, feature2_id, feature3_id)	{ (feat)[(ind)] = (val); ind++; }
#endif

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
double *SVMBehaviorSequence::psi_bout(BehaviorBoutFeatures *b, int t_start, int t_end, int beh, int c, 
				      double *feat, bool normalize, bool fast_update, 
				      double get_extreme_vals[2][NUMFEAT], double use_extreme_vals[2][NUMFEAT]) {
  if(!feat) feat = (double*)malloc(num_features*sizeof(double));

  // Eyrun: temporary fix
  if(t_start == t_end)
    t_end += 1;

  int i, j, k, l, ind = 0, td, beg;
  double start, temporal_grid_size[MAX_TEMPORAL_LEVELS], harmonic_grid_size[MAX_HARMONIC_LEVELS+1], sum, f, t, ave, sum_sqr;
  SVMFeatureParams *p;
  double dur = t_end-t_start;
  double inv = dur ? 1/dur : 0;

  double max_vals[NUMFEAT], min_vals[NUMFEAT];

  //p = &feature_params[0];
  //if(p->use_bout_max_feature | p->use_bout_min_feature) { // Eyrun: only update cache if needed (note this assumes that all feature stats are the same as for 0)
    if(!fast_update && !use_extreme_vals)
      b->bout_start = b->bout_end = -1;
    else
      assert(b->bout_end == t_end);

    if(use_extreme_vals) {
      // The min/max values for this bout were passed into this function, such that we don't need to compute them
      for(i=0; i<num_base_features; i++) {
	max_vals[i] = use_extreme_vals[0][i];
	min_vals[i] = use_extreme_vals[1][i];
      } 
    } else {
      // Update the min/max values for this bout
      b->UpdateCaches(t_start, t_end, c);
      for(i=0; i<num_base_features; i++) {
	max_vals[i] = b->bout_max_feature_responses[i];
	min_vals[i] = b->bout_min_feature_responses[i];
      }
    }
    if(get_extreme_vals) {
      // Return the min/max values for this bout back to the caller
      for(i=0; i<num_base_features; i++) {
	get_extreme_vals[0][i] = max_vals[i];
	get_extreme_vals[1][i] = min_vals[i];
      }
    }
    // }

  assert(dur > 0);
  //double log_dur = log10(dur);
  //ADD_FEATURE(feat, ind, log_dur, true, -1, 0);
  for(i = 0; i < MAX_TEMPORAL_LEVELS; i++) 
    temporal_grid_size[i] = (t_end-t_start) / ((double)(1<<i));
  for(i = 0; i < MAX_HARMONIC_LEVELS+1; i++) 
    harmonic_grid_size[i] = (t_end-t_start) / ((double)(i+1));

  for(i = 0; i < num_base_features; i++) {
    p = &feature_params[i];
    beg = ind;

    sum = (b->integral_features[i][t_end]-b->integral_features[i][t_start]);
    ave = sum * inv;

    // Sum feature response in the interval (t_start,t_end).  Can also divide into multiple temporal
    // levels, and compute the average response in grids of size (t_end-t_start)/2, (t_end-t_start)/4, 
    // (t_end-t_start)/8, ... This is an average feature response in temporal intervals of different size,
    // where the outputted feature vector is the concatenation of features over each temporal interval
    // The output feature space has size 
    //   F_NUM_GOOD_GLOBAL_FEATURES*(sum_{l=0}^num_temporal_levels(2^i))
    for(j = 0; j < p->num_temporal_levels; j++) {
      for(k = 0, start = t_start; k < (1<<j); k++, start += temporal_grid_size[j]) {
	f = (b->integral_features[i][(int)(start+temporal_grid_size[j]+.5)] -
	     b->integral_features[i][(int)(start+.5)]);
	if(p->use_bout_sum_features) ADD_FEATURE(feat, ind, f, USE_BOUT_SUM_FEATURES, i, j);
	if(p->use_bout_ave_features) ADD_FEATURE(feat, ind, f*inv, USE_BOUT_AVE_FEATURES, i, j);
	if(p->use_bout_sum_absolute_features) ADD_FEATURE(feat, ind, my_abs(f), USE_BOUT_SUM_ABSOLUTE_FEATURES, i, j);
	if(p->use_bout_ave_absolute_features) ADD_FEATURE(feat, ind, my_abs(f)*inv, USE_BOUT_AVE_ABSOLUTE_FEATURES, i, j);

// 				// Add max and min feature responses over the temporal region
// 				double min = INFINITY;
// 				double max = -INFINITY;
// 				int f_end = my_min((int)(start+temporal_grid_size[j]),t_end);
// 				int f_start = my_max((int)(start),t_start);
// 				for (int ff=f_start; ff<=f_end; ff++) {
// 					min = my_min(min,b->smoothed_features[i][ff]);
// 					max = my_max(max,b->smoothed_features[i][ff]);
// 				}
// 				if(p->use_bout_max_feature) ADD_FEATURE(feat, ind, max, USE_BOUT_MAX_FEATURES, i, j);
// 				if(p->use_bout_min_feature) ADD_FEATURE(feat, ind, min, USE_BOUT_MIN_FEATURES, i, j);
      }
    }

    // The standard deviation of the feature in the interval (t_start,t_end),  Can also use the total sum variance
    if(p->use_standard_deviation || p->use_sum_variance) {
      sum_sqr = (b->integral_sqr_features[i][t_end]-b->integral_sqr_features[i][t_start]);
      f = sum_sqr - 2*ave*sum + SQR(ave)*dur;
      if(f < 0) f = 0;  // avoid precision-related errors
      if(p->use_sum_variance) ADD_FEATURE(feat, ind, f, USE_SUM_VARIANCE, i, 0);
      if(p->use_standard_deviation) ADD_FEATURE(feat, ind, sqrt(f*inv), USE_STANDARD_DEVIATION, i, 0);
    }

    // Add features for the min/max feature response in the current bout, or thresholded versions of the min/max values
    if(p->use_bout_max_feature) 
      ADD_FEATURE(feat, ind, max_vals[i], USE_BOUT_MAX_FEATURE, i, j);
    if(p->use_bout_min_feature) 
      ADD_FEATURE(feat, ind, min_vals[i], USE_BOUT_MIN_FEATURE, i, j);
    for(l = 0; l < p->num_bout_max_thresholds; l++) 
      ADD_FEATURE(feat, ind, max_vals[i] > max_thresholds[i][l] ? 1 : 0, NUM_BOUT_MAX_THRESHOLDS, i, j);
    for(l = 0; l < p->num_bout_min_thresholds; l++) 
      ADD_FEATURE(feat, ind, min_vals[i] < min_thresholds[i][l] ? 1 : 0, NUM_BOUT_MIN_THRESHOLDS, i, j);


    // These are differences in the bout average or sum response for the bout as compared to the global max, min,
    // and average over the entire behavioral sequence
    if(p->use_global_difference_max_sum_features) 
      ADD_FEATURE(feat, ind, sum - max_vals[i]*(t_end-t_start), USE_GLOBAL_DIFFERENCE_MAX_SUM_FEATURES, i, 0);
    if(p->use_global_difference_max_ave_features) 
      ADD_FEATURE(feat, ind, ave - max_vals[i], USE_GLOBAL_DIFFERENCE_MAX_AVE_FEATURES, i, 0);
    if(p->use_global_difference_min_sum_features) 
      ADD_FEATURE(feat, ind, sum - min_vals[i]*(t_end-t_start), USE_GLOBAL_DIFFERENCE_MIN_SUM_FEATURES, i, 0);
    if(p->use_global_difference_min_ave_features) 
      ADD_FEATURE(feat, ind, ave - min_vals[i], USE_GLOBAL_DIFFERENCE_MIN_AVE_FEATURES, i, 0);
    if(p->use_global_difference_ave_sum_features) 
      ADD_FEATURE(feat, ind, sum - b->ave_feature_responses[i]*(t_end-t_start), USE_GLOBAL_DIFFERENCE_AVE_SUM_FEATURES, i, 0);
    if(p->use_global_difference_ave_ave_features) 
      ADD_FEATURE(feat, ind, ave - b->ave_feature_responses[i], USE_GLOBAL_DIFFERENCE_AVE_AVE_FEATURES, i, 0);


    // The total change from the beginning of the bout to the end of the bout
    if(p->num_bout_change_points) {
      for(j = 0, t=t_start; j < p->num_bout_change_points; j++, t+=harmonic_grid_size[p->num_bout_change_points]) {
	f = b->features[i][my_max((int)t,t_end-1)] - b->features[i][t_start];
	if(p->use_bout_change) 
	  ADD_FEATURE(feat, ind, f, USE_BOUT_CHANGE, i, j);
	if(p->use_bout_absolute_change) 
	  ADD_FEATURE(feat, ind, my_abs(f), USE_BOUT_ABSOLUTE_CHANGE, i, j);
      }
    }


    // Discretize the feature space into different bins, and then compute the count of the number
    // of frames in each bin.  This is a histogram of feature values.  The histograms can also be
    // computed over different temporal intervals, yielding an output feature space of size
    //    F_NUM_GOOD_GLOBAL_FEATURES*num_histogram_bins*(sum_{l=0}^num_histogram_bins_temporal_levels(2^i))
    for(j = 0; j < p->num_histogram_temporal_levels; j++) {
      for(l = 0; l < p->num_histogram_bins; l++) {
	for(k = 0, start = t_start; k < (1<<j); k++, start += temporal_grid_size[j]) {
	  f = (b->integral_histogram_features[i][l][(int)(start+temporal_grid_size[j]+.5)] -
	       b->integral_histogram_features[i][l][(int)(start+.5)]);
	  if(p->use_histogram_sum_features) ADD_FEATURE(feat, ind, f, USE_HISTOGRAM_SUM_FEATURES, i, j);
	  if(p->use_histogram_ave_features) ADD_FEATURE(feat, ind, f*inv, USE_HISTOGRAM_AVE_FEATURES, i, j);
	}
      }
    }

    // Kind of like 1D Haar-like features, these are the difference in sum feature response
    // in some temporal interval minus the difference in sum feature response of some other interval.
    // Can also use features of the absolute value of these feature responses
    // The total num of features (if absolute features are also used) is
    //   F_NUM_GOOD_GLOBAL_FEATURES*num_difference_temporal_levels*2*2
    for(j = 0; j <  p->num_difference_temporal_levels; j++) {
      // Difference in sum feature response in the region inside the bout (t_start,t_end) and 
      // the region of duration (t_end-t_start)/(2^j) immediately before t_start

      td = (int)(temporal_grid_size[j+1]+.5);
      f = b->integral_features[i][t_start+td] - 2*b->integral_features[i][t_start] + b->integral_features[i][t_start-td];
      if(p->use_start_sum_diff_haar_features)
	ADD_FEATURE(feat, ind, f, USE_START_SUM_DIFF_HAAR_FEATURES, i, 0);
      if(p->use_start_ave_diff_haar_features)
	ADD_FEATURE(feat, ind, f*inv, USE_START_AVE_DIFF_HAAR_FEATURES, i, 0);
      if(p->use_start_sum_absolute_diff_haar_features) 
	ADD_FEATURE(feat, ind, my_abs(f), USE_START_SUM_ABSOLUTE_DIFF_HAAR_FEATURES, i, 0);
      if(p->use_start_ave_absolute_diff_haar_features) 
	ADD_FEATURE(feat, ind, my_abs(f*inv), USE_START_AVE_ABSOLUTE_DIFF_HAAR_FEATURES, i, 0);


      // Difference in average feature response in the region inside the bout (t_start,t_end) and 
      // the region of duration (t_end-t_start)/(2^j) immediately after t_end
      f = -b->integral_features[i][t_end+td] + 2*b->integral_features[i][t_end] - b->integral_features[i][t_end-td];
      if(p->use_end_sum_diff_haar_features)
	ADD_FEATURE(feat, ind, f, USE_END_SUM_DIFF_HAAR_FEATURES, i, 0);
      if(p->use_end_ave_diff_haar_features)
	ADD_FEATURE(feat, ind, f*inv, USE_END_AVE_DIFF_HAAR_FEATURES, i, 0);
      //ADD_FEATURE(feat, ind, f*inv, USE_END_AVE_DIFF_HAAR_FEATURES);
      if(p->use_end_sum_absolute_diff_haar_features) 
	ADD_FEATURE(feat, ind, my_abs(f), USE_END_SUM_ABSOLUTE_DIFF_HAAR_FEATURES, i, 0);
      if(p->use_end_ave_absolute_diff_haar_features) 
	ADD_FEATURE(feat, ind, my_abs(f*inv), USE_END_AVE_ABSOLUTE_DIFF_HAAR_FEATURES, i, 0);
    }

    // An extension of the 1D Haar-like features, these capture changes or harmonic motion within the
    // bout.  When j=0, the response is the 1st half of the bout minus the end half of the bout.  When
    // j=1, the response is the 1st 3rd minus the 2nd 3rd plus the 3rd 3rd. And so on.
    // The total num of features (if absolute features are also used) is
    //   F_NUM_GOOD_GLOBAL_FEATURES*num_harmonic_features*2
    for(j = 1; j <= p->num_harmonic_features; j++) {
      f = 0;
      for(k = 0, start = t_start; k < j+1; k++, start += harmonic_grid_size[j]) {
	f += (1-(k&1)*2)*(b->integral_features[i][(int)(start+harmonic_grid_size[j]+.5)] -
			  b->integral_features[i][(int)(start+.5)]);
      }
      if(p->use_sum_harmonic_features)
	ADD_FEATURE(feat, ind, f, USE_SUM_HARMONIC_FEATURES, i, j);
      if(p->use_ave_harmonic_features)
	ADD_FEATURE(feat, ind, f*inv, USE_AVE_HARMONIC_FEATURES, i, j);
      if(p->use_sum_absolute_harmonic_features)
	ADD_FEATURE(feat, ind, my_abs(f), USE_SUM_ABSOLUTE_HARMONIC_FEATURES, i, j);
      if(p->use_ave_absolute_harmonic_features)
	ADD_FEATURE(feat, ind, my_abs(f*inv), USE_AVE_ABSOLUTE_HARMONIC_FEATURES, i, j);
    }
    assert(ind-beg == p->num_features);
  }

  // TODO: add duration features, taking the duration (dividing 
  // log(frame_times[t_end]-frame_times[t_start]) into p->num_duration_features bins based on median statistics of train set)


  assert(ind == num_features);

  // Normalize each feature to have mean 0 and standard deviation 1 (as determined by the training set)
  //normalize = false;
  if(normalize) {
    for(i = 0; i < ind; i++)  {
      //assert(!isnan(feat[i]));
      //assert(abs(feat[i]) < 100000000);  // EYRUN: inrcreased to 100,000,000 from 10,000,000
      feat[i] = (feat[i]-features_mu[i])*features_gamma[i];
      //double div = my_max(b->max_feature_responses[i],-b->min_feature_responses[i]);
      //feat[i] = feat[i]/(div < 0.001 ? 1 : div);
      //assert(!isnan(feat[i]));
      //assert(abs(feat[i]) < 100000000);
#ifdef CAP_BOUT_FEATURES
      if(feat[i] < -CAP_BOUT_FEATURES) feat[i] = -CAP_BOUT_FEATURES;
      if(feat[i] > CAP_BOUT_FEATURES) feat[i] = CAP_BOUT_FEATURES;
#endif
    }
  }

  // For a given behavior class, optionally use only a subset of the available bout features
  if(beh >= 0 && c >= 0 && restrict_behavior_features[beh] && restrict_behavior_features[beh][c]) {
    for(i = 0; i < ind; i++) 
      if(!restrict_behavior_features[beh][c][i])
	feat[i] = 0;
  }

// 	//EYRUN: test (weigh the bout feature with the duration of the interval)
// 	for (i=0; i<num_features; i++)
// 		feat[i] = feat[i] * dur / 20000.0;

	return feat;
}


/*
* This is called to update the bout feature cache structures (like the bout-wise min and max)
* in a computationally efficient way, e.g. when doing dynamic programming we can update the
* min/max from the bout explored in the previous loop iteration in constant time (as opposed
* to time proportional to the number of frames in the bout)
*/
void BehaviorBoutFeatures::UpdateCaches(int t_start, int t_end, int c) {
	/*
	  BehaviorBoutFeatures *b: contains b-> start and b->and, which contain

	  Two usage scenarios:
	  (a) ev
	  (b) within 
	*/
	bool reset = false;
	int i, j;

	if(t_end == bout_end) {
		if(t_start > bout_start) {
			printf("Upcoming Error...\n");
			//reset = true;
			//bout_start = t_end;
		}
		assert(t_start <= bout_start);
	} else {
		reset = true;
		bout_start = bout_end = t_end;
	}
	for(i = 0; i < num_base_features; i++) {
		if(reset) {
			bout_max_feature_responses[i] = -INFINITY;
			bout_min_feature_responses[i] = INFINITY;
		}    
		for(j = t_start; j < bout_start; j++) {
			bout_max_feature_responses[i] = my_max(smoothed_features[i][j],  bout_max_feature_responses[i]);
			bout_min_feature_responses[i] = my_min(smoothed_features[i][j],  bout_min_feature_responses[i]);
		}
	}
	bout_start = t_start;
}

SparseVector SVMBehaviorSequence::Psi(StructuredData *x, StructuredLabel *yy) {
  BehaviorBoutSequence *y = (BehaviorBoutSequence*)yy, *y_partial = NULL;
  BehaviorBoutFeatures *b = (BehaviorBoutFeatures*)x;
  int i, j, beh;
  double *tmp_features = (double*)malloc(sizeof(double)*(sizePsi+num_features+10));
  double *all_features = tmp_features+num_features+10;
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

  for(beh = 0; beh < behaviors->num; beh++) {
    if(behavior >= 0 && beh != behavior) {
      ptr += getPsiSize(num_features, num_classes[beh]);
      continue;
    }
    class_features = ptr; ptr += num_features*num_classes[beh];        // line 1
    class_transitions = ptr; ptr += num_classes[beh]*num_classes[beh]; // line 2
    class_counts = ptr; ptr += num_classes[beh];   // unary feature
#if USE_DURATION_COST > 0
    class_durations = ptr; ptr += num_classes[beh];   // duration feature
#endif

    for(i = 0; i < num_classes[beh]*num_features; i++) 
      class_features[i] = 0;

    for(i = 0; i < num_classes[beh]*num_classes[beh]; i++) 
      class_transitions[i] = 0;

    for(i = 0; i < num_classes[beh]; i++) 
      class_counts[i] = 0;

#if USE_DURATION_COST > 0
    for(i = 0; i < num_classes[beh]; i++) 
      class_durations[i] = 0;
#endif

    for(i = 0; i < y->num_bouts[beh]; i++) {
      if (i) // count the number of times a transition from class c_p to c occurs
        class_transitions[y->bouts[beh][i-1].behavior*num_classes[beh] + y->bouts[beh][i].behavior]++;   
      if(!NONE_CLASS_HAS_NO_SCORE || y->bouts[beh][i].behavior != 0)
	class_counts[y->bouts[beh][i].behavior]++;     // count the number of times each class c occurs
#if USE_DURATION_COST > 0
      duration = y->bouts[beh][i].end_frame-y->bouts[beh][i].start_frame;
      duration_diff = 0;
#if USE_DURATION_COST > 1
      if (duration < min_frame_duration[beh][y->bouts[beh][i].behavior]) duration_diff = 10;
      else if (duration > max_frame_duration[beh][y->bouts[beh][i].behavior]) duration_diff = 10;
#else
      if (duration < min_frame_duration[beh][y->bouts[beh][i].behavior]) 
	duration_diff = min_frame_duration[beh][y->bouts[beh][i].behavior] - duration;
      else if (duration > max_frame_duration[beh][y->bouts[beh][i].behavior]) 
	duration_diff = duration - max_frame_duration[beh][y->bouts[beh][i].behavior];
      duration_diff *= duration_diff;
#endif
      if(!NONE_CLASS_HAS_NO_SCORE || y->bouts[beh][i].behavior != 0)
	class_durations[y->bouts[beh][i].behavior] += duration_diff;
#endif
      // The main appearance features used to compute the score of a bout for a given class.
      // Since the total score is the sum of the scores over bouts, we can simply add together
      // the features of the bouts with the same class label
      if(!NONE_CLASS_HAS_NO_SCORE || y->bouts[beh][i].behavior != 0) {
	psi_bout(b, y->bouts[beh][i].start_frame, y->bouts[beh][i].end_frame, beh, y->bouts[beh][i].behavior, tmp_features, true, false);
	for(j = 0; j < num_features; j++) 
	  class_features[y->bouts[beh][i].behavior*num_features+j] += tmp_features[j];
      }
    }  
  }
        
  SparseVector retval = SparseVector(all_features, sizePsi);
  free(tmp_features);
  if(y != yy) delete y;
  if(y_partial) delete y_partial;

  return retval;
}


void SVMBehaviorSequence::set_feature_name(int feature_ind, int base_feature_ind, const char *name) {
	char str[1000];
	if(feature_names) {
		sprintf(str, "%s %s", get_base_feature_name(base_feature_ind), name);
		feature_names[feature_ind] = StringCopy(str);
	}
}

void SVMBehaviorSequence::print_features(FILE *fout, double *feat) {
	for(int i = 0; i < num_features; i++)
		fprintf(fout, "    %lf %s\n", feat[i], feature_names[i]);
}

void SVMBehaviorSequence::print_features(const char *fname, StructuredDataset *dataset, bool normalized) {
	double *tmp_features = (double*)malloc(sizeof(double)*num_features);
	int i, j, beh;
	BehaviorBoutSequence *y;
	BehaviorBoutFeatures *x;
	StructuredExample **ex = dataset->examples;
	FILE *fout = fopen(fname, "w");
	assert(fout);

	for(i = 0; i < dataset->num_examples; i++) {
		y = (BehaviorBoutSequence*)ex[i]->y;
		x = (BehaviorBoutFeatures*)ex[i]->x;
		for(beh = 0; beh < behaviors->num; beh++) {
			if(behavior >= 0 && beh != behavior) 
				continue;
			for(j = 0; j < y->num_bouts[beh]; j++) {
				fprintf(fout, "  bout %d-%d %s:%s\n", y->bouts[beh][j].start_frame, y->bouts[beh][j].end_frame, behaviors->behaviors[beh].name, 
					behaviors->behaviors[beh].values[y->bouts[beh][j].behavior].name);
				psi_bout(x, y->bouts[beh][j].start_frame, y->bouts[beh][j].end_frame, beh, y->bouts[beh][j].behavior, tmp_features, normalized, false);
				print_features(fout, tmp_features);
			}
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
	int beh, i, j;
	int *inds = (int*)malloc(sizeof(int)*(num_features+100000));

	for(beh = 0; beh < behaviors->num; beh++) {
		if(behavior >= 0 && beh != behavior) {
			ptr += num_classes[beh]*(num_classes[beh]+num_features);
			continue;
		}
		class_features = ptr; ptr += num_features*num_classes[beh];
		class_transitions = ptr; ptr += num_classes[beh]*num_classes[beh];
		g_num_features = num_features;

		for(i = 0; i < behaviors->behaviors[beh].num_values; i++, class_features += num_features, class_transitions += behaviors->behaviors[beh].num_values) {
			for(j = 0; j < num_features+behaviors->behaviors[beh].num_values; j++)
				inds[j] = j;
			g_class_features = class_features;
			g_class_transitions = class_transitions;
			qsort(inds, num_features+behaviors->behaviors[beh].num_values, sizeof(int), cmp_feature_inds);
			for(j = 0; j < num_features+behaviors->behaviors[beh].num_values; j++) {
				if(inds[j] < num_features)
					fprintf(fout, "%lf %s:%s %s\n", class_features[inds[j]], behaviors->behaviors[beh].name, behaviors->behaviors[beh].values[i].name, feature_names[inds[j]]);
				else
					fprintf(fout, "%lf %s:%s->%s\n", class_transitions[inds[j]-num_features], behaviors->behaviors[beh].name, behaviors->behaviors[beh].values[i].name, behaviors->behaviors[beh].values[inds[j]-num_features].name);
			}
		}
	}
	free(inds);
}

void SVMBehaviorSequence::print_weights(const char *fname, double *weights) {
	FILE *fout = fopen(fname, "w");
	print_weights(fout, weights);
	fclose(fout);
}

/*
* Compute the size of the bout-level feature space.  This should be identical to psi_bout(), but without computing
* any features
*/
int SVMBehaviorSequence::compute_feature_space_size() {
	int i, j, k, l, curr;
	SVMFeatureParams *p;
	char str[1000];

	// Compute the total number of features used
	num_features = 0;
	//sprintf(str,"log duration");
	//set_feature_name(num_features++,-1,str);
	for(i = 0; i < num_base_features; i++) {
		p = &feature_params[i];
		curr = num_features;

		for(j = 0; j < p->num_temporal_levels; j++) {
			for(k = 0; k < (1<<j); k++) {
				if(p->use_bout_sum_features) { sprintf(str, "bout_sum_%.3f<=t<%.3f", k/(float)(1<<j), (k+1)/(float)(1<<j)); set_feature_name(num_features++, i, str); }
				if(p->use_bout_ave_features) { sprintf(str, "bout_average_%.3f<=t<%.3f", k/(float)(1<<j), (k+1)/(float)(1<<j)); set_feature_name(num_features++, i, str); }
				if(p->use_bout_sum_absolute_features) { sprintf(str, "bout_sum_abs_%.3f<=t<%.3f", k/(float)(1<<j), (k+1)/(float)(1<<j)); set_feature_name(num_features++, i, str); }
				if(p->use_bout_ave_absolute_features) { sprintf(str, "bout_average_abs_%.3f<=t<%.3f", k/(float)(1<<j), (k+1)/(float)(1<<j)); set_feature_name(num_features++, i, str); }
// 				if(p->use_bout_max_feature) { sprintf(str, "bout_max_%.3f<=t<%.3f", k/(float)(1<<j), (k+1)/(float)(1<<j)); set_feature_name(num_features++, i, str); }
// 				if(p->use_bout_min_feature) { sprintf(str, "bout_min_%.3f<=t<%.3f", k/(float)(1<<j), (k+1)/(float)(1<<j)); set_feature_name(num_features++, i, str); }
			}
		}
		if(p->use_sum_variance) set_feature_name(num_features++, i, "sum_variance"); 
		if(p->use_standard_deviation) set_feature_name(num_features++, i, "standard_deviation"); 
		if(p->use_bout_max_feature) set_feature_name(num_features++, i, "bout_max"); 
		if(p->use_bout_min_feature) set_feature_name(num_features++, i, "bout_min"); 
		for(l = 0; l < p->num_bout_max_thresholds; l++) { sprintf(str, "bout_max>%.5f", max_thresholds ? max_thresholds[i][l] : 0); set_feature_name(num_features++, i, str); }
		for(l = 0; l < p->num_bout_min_thresholds; l++) { sprintf(str, "bout_min<%.5f", min_thresholds ? min_thresholds[i][l] : 0); set_feature_name(num_features++, i, str); }
		if(p->use_global_difference_max_sum_features) set_feature_name(num_features++, i, "sum_bout_difference_from_global_max"); 
		if(p->use_global_difference_max_ave_features) set_feature_name(num_features++, i, "ave_bout_difference_from_global_max"); 
		if(p->use_global_difference_min_sum_features) set_feature_name(num_features++, i, "sum_bout_difference_from_global_min"); 
		if(p->use_global_difference_min_ave_features) set_feature_name(num_features++, i, "ave_bout_difference_from_global_min"); 
		if(p->use_global_difference_ave_sum_features) set_feature_name(num_features++, i, "sum_bout_difference_from_global_ave"); 
		if(p->use_global_difference_ave_ave_features) set_feature_name(num_features++, i, "ave_bout_difference_from_global_ave");

		for(j = 0; j < p->num_bout_change_points; j++) {
			if(p->use_bout_change) { sprintf(str, "bout_change_through_t=%.5f", (j+1)/(float)p->num_bout_change_points); set_feature_name(num_features++, i, str); }
			if(p->use_bout_absolute_change) { sprintf(str, "bout_absolute_change_through_t=%.5f", (j+1)/(float)p->num_bout_change_points); set_feature_name(num_features++, i, str); }
		}

		for(j = 0; j < p->num_histogram_temporal_levels; j++) {
			for(l = 0; l < p->num_histogram_bins; l++) {
				for(k = 0; k < (1<<j); k++) {
					if(p->use_histogram_sum_features) { 
						if(histogram_thresholds) sprintf(str, "histogram_sum_%.4f<=f<%.4f_%.3f<=t<%.3f", 
							(float)(l ? histogram_thresholds[i][l-1] : -INFINITY), (float)histogram_thresholds[i][l], k/(float)(1<<j), (k+1)/(float)(1<<j)); 
						set_feature_name(num_features++, i, str);
					}
					if(p->use_histogram_ave_features) {
						if(histogram_thresholds) sprintf(str, "histogram_normalized_%.4f<=f<%.4f_%.3f<=t<%.3f", 
							(float)(l ? histogram_thresholds[i][l-1] : -INFINITY), (float)histogram_thresholds[i][l], k/(float)(1<<j), (k+1)/(float)(1<<j)); 
						set_feature_name(num_features++, i, str);
					}
				}
			}
		}

		for(j = 0; j <  p->num_difference_temporal_levels; j++) {
			if(p->use_start_sum_diff_haar_features) { sprintf(str, "haar_sum_(0<=t<%.3f)-(-%.3f<=t<0)", 1.0f/(1<<(j+1)), 1.0f/(1<<(j+1))); set_feature_name(num_features++, i, str); }
			if(p->use_start_ave_diff_haar_features) { sprintf(str, "haar_normalized_(0<=t<%.3f)-(-%.3f<=t<0)", 1.0f/(1<<(j+1)), 1.0f/(1<<(j+1))); set_feature_name(num_features++, i, str); }
			if(p->use_start_sum_absolute_diff_haar_features) { sprintf(str, "haar_sum_|(0<=t<%.3f)-(-%.3f<=t<0)|", 1.0f/(1<<(j+1)), 1.0f/(1<<(j+1))); set_feature_name(num_features++, i, str); }
			if(p->use_start_ave_absolute_diff_haar_features) { sprintf(str, "haar_normalized_|(0<=t<%.3f)-(-%.3f<=t<0)|", 1.0f/(1<<(j+1)), 1.0f/(1<<(j+1))); set_feature_name(num_features++, i, str); }

			if(p->use_end_sum_diff_haar_features) { sprintf(str, "haar_sum_(%.3f<=t<1)-(1<=t<%.3f)", 1.0f-1.0f/(1<<(j+1)), 1.0f+1.0f/(1<<(j+1))); set_feature_name(num_features++, i, str); }
			if(p->use_end_ave_diff_haar_features) { sprintf(str, "haar_normalized_|(%.3f<=t<1)-(1<=t<%.3f)|", 1.0f-1.0f/(1<<(j+1)), 1.0f+1.0f/(1<<(j+1))); set_feature_name(num_features++, i, str); }
			if(p->use_end_sum_absolute_diff_haar_features) { sprintf(str, "haar_sum_(%.3f<=t<1)-(1<=t<%.3f)", 1.0f-1.0f/(1<<(j+1)), 1.0f+1.0f/(1<<(j+1))); set_feature_name(num_features++, i, str); }
			if(p->use_end_ave_absolute_diff_haar_features) { sprintf(str, "haar_normalized_|(%.3f<=t<1)-(1<=t<%.3f)|", 1.0f-1.0f/(1<<(j+1)), 1.0f+1.0f/(1<<(j+1))); set_feature_name(num_features++, i, str); }
		}

		for(j = 1; j <= p->num_harmonic_features; j++) {
			if(p->use_sum_harmonic_features) { sprintf(str, "harmonic_sum_%d", j+1); set_feature_name(num_features++, i, str); }
			if(p->use_ave_harmonic_features) { sprintf(str, "harmonic_normalized_%d", j+1); set_feature_name(num_features++, i, str); }
			if(p->use_sum_absolute_harmonic_features) { sprintf(str, "harmonic_sum_absolute_%d", j+1); set_feature_name(num_features++, i, str); }
			if(p->use_ave_absolute_harmonic_features) { sprintf(str, "harmonic_normalized_absolute_%d", j+1); set_feature_name(num_features++, i, str); }
		}

		p->num_features = num_features - curr;
	} 

	return num_features;
}


bool SVMBehaviorSequence::SaveDataset(StructuredDataset *d, const char *fname, int start_from) {
  save_examples(fname, d);
  return true;
}

StructuredDataset *SVMBehaviorSequence::LoadDataset(const char *fname) {
	if(debugLevel > 0) fprintf(stderr, "Reading dataset %s...", fname);
  
	Lock();

	if(strlen(debugdir) && debug_features) 
	  CreateDirectoryIfNecessary(debugdir);
	      

	int num, i, j, beh;       
	char **train_list = load_examples(fname, &num);
	bool computeClassTransitions = class_training_transitions ? false : true;
	//bool computeClassTransitions = true; // TEMP: EYRUN
	StructuredDataset *dataset = new StructuredDataset();
	int duration = 0;

	// Keep track of the number of transitions between each pair of classes
	if(computeClassTransitions) {
		class_training_transitions = (int***)malloc(behaviors->num*sizeof(int**));
		class_training_transitions_count = (int**)malloc(behaviors->num*sizeof(int*));
		class_training_count = (int**)malloc(behaviors->num*sizeof(int*));
		for(beh = 0; beh < behaviors->num; beh++) {
			class_training_transitions[beh] = (int**)malloc(num_classes[beh]*sizeof(int*));
			class_training_transitions_count[beh] = (int*)malloc(num_classes[beh]*sizeof(int));
			class_training_count[beh] = (int*)malloc(num_classes[beh]*sizeof(int));
			for(i = 0; i < num_classes[beh]; i++) {
				class_training_transitions[beh][i] = (int*)malloc(num_classes[beh]*sizeof(int));
				class_training_count[beh][i] = class_training_transitions_count[beh][i] = 0;
				for(j = 0; j < num_classes[beh]; j++) {
					class_training_transitions[beh][i][j] = 0;
				}
			}
		}
	}

	for(j = 0; j < num; j++) {
		StructuredExample *ex = read_struct_example(train_list[j], train_list[j], false);
		dataset->AddExample(ex);
		BehaviorBoutSequence *y = (BehaviorBoutSequence*)ex->y;
		if(computeClassTransitions) {
			for(beh = 0; beh < behaviors->num; beh++) {
				for(i = 0; i < y->num_bouts[beh]; i++) {
					class_training_count[beh][y->bouts[beh][i].behavior]++;
#if USE_DURATION_COST > 0
					duration = y->bouts[beh][i].end_frame - y->bouts[beh][i].start_frame;
					min_frame_duration[beh][y->bouts[beh][i].behavior] = my_min(min_frame_duration[beh][y->bouts[beh][i].behavior], duration);
					max_frame_duration[beh][y->bouts[beh][i].behavior] = my_max(max_frame_duration[beh][y->bouts[beh][i].behavior], duration);
#endif
					if(i)
						class_training_transitions[beh][y->bouts[beh][i].behavior][y->bouts[beh][i-1].behavior]++;
				}
			}	 
		} 
		// Temporary fix, ground truth segmentations should have bouts that span all the way to b->num_frames
		int T = ((BehaviorBoutFeatures*)ex->x)->num_frames;
		if(y->num_bouts && y->num_bouts[0] && y->bouts[0][y->num_bouts[0]-1].end_frame==T-1)
		  y->bouts[0][y->num_bouts[0]-1].end_frame = T;
	}
	if(computeClassTransitions) {
		printf("Training examples %d behavior sequences\n", (int)num);


		// Compute all feature cache data structures for all training examples.  This requires first computing some
		// statistics of the training set features for normalization purposes
		compute_feature_mean_variance_median_statistics(dataset);

		// Compute features for each training bout, normalized to have (0,1) mean and standard deviation
		for(j = 0; j < num; j++) {
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


		// Compress the class transitions into a sparse array
		for(beh = 0; beh < behaviors->num; beh++) {
			for(i = 0; i < num_classes[beh]; i++) {
				class_training_transitions_count[beh][i] = 0;
				for(j = 0; j < num_classes[beh]; j++) {
					if(class_training_transitions[beh][i][j])
						class_training_transitions[beh][i][class_training_transitions_count[beh][i]++] = j;
				}
			}
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
    int beh = this->behavior >= 0 ? this->behavior : 0;
    if(debug_features) {
      if(!ww) ww = w->get_non_sparse<double>(sizePsi);
      if(iter_num >= 0) sprintf(fname, "%s/%s_%d_features.html", debugdir, file, iter_num);
      else sprintf(fname, "%s/%s_features.html", debugdir, file);
      DebugFeatures(fname, (BehaviorBoutSequence*)ybar, ww, beh);
    }
    if(iter_num >= 0) sprintf(fname, "%s/%s_%d", debugdir, file, iter_num);
    else sprintf(fname, "%s/%s", debugdir, file);
    ((BehaviorBoutSequence*)ybar)->Visualize(behaviors, beh, fname, html);
    if(y) {
      BehaviorBoutSequence *yy = (BehaviorBoutSequence*)y;
      if(!yy->score) {
	if(ybar) yy->loss = Loss(y, ybar);
	
      }
      if(debug_features) {
	sprintf(fname, "%s/%s_gt_%d_features.html", debugdir, file, iter_num);
	DebugFeatures(fname, yy, ww, beh);
      }
      sprintf(fname, "%s/%s_gt_%d", debugdir, file, iter_num);
      ((BehaviorBoutSequence*)y)->Visualize(behaviors, beh, fname, html_gt);
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
	int lastChar = strlen(filepath) - 1;

	for(int pos=lastChar; pos >=0; pos--) {
		if(filepath[pos]=='\\' || filepath[pos]=='/')
			return filepath+pos+1;
	}

	return filepath;
}




// Helper function for Inference().  Allocate memory for the label we are returning
void SVMBehaviorSequence::init_bout_label(BehaviorBoutSequence *ybar, BehaviorBoutSequence *y) {
  int sz = (sizeof(BehaviorBout*)+sizeof(int*)+2*sizeof(double))*behaviors->num;
  ybar->bouts = (BehaviorBout**)realloc(ybar->bouts, sz);
  memset(ybar->bouts, 0, sz);
  ybar->num_bouts = (int*)(ybar->bouts+behaviors->num);
  ybar->scores = (double*)(ybar->num_bouts+behaviors->num);
  ybar->losses = (double*)(ybar->scores+behaviors->num);
  ybar->slack = ybar->score = 0;
  ybar->loss = 0;
  ybar->behaviors = behaviors;
  if(y) {
    ybar->features = y->features;
    ybar->behaviors = y->behaviors;
    y->score = 0;
  }
}

// TODO: change this to be the hamming distance
// If a ground truth label y is given, update the componenent of loss(y,ybar) that is 
// attributable to the completed bout (c_prev,t_p,t)
double SVMBehaviorSequence::compute_updated_bout_loss(BehaviorBoutFeatures *b, BehaviorBoutSequence *y, 
						      int beh, int T, int t_p, int t, int c_prev, double *fn, 
						      int *gt_bout, double *dur_gt, double &loss_fp, double &loss_fn) {
  double loss_score = 0, fp;
  int TT = my_min(T,b->num_frames-1);  // int TT = T-1;

  if(y) { 
#if HAMMING_LOSS == 1
    loss_fp = loss_fn = 0;
    for(int j = gt_bout[t_p]; j <= gt_bout[t] && j < y->num_bouts[beh]; j++) {
      // Loop through all ground truth bouts intersecting with the completed bout
      double inter = (b->frame_times[my_min(my_min(y->bouts[beh][j].end_frame,t),TT)] - 
                      b->frame_times[my_max(y->bouts[beh][j].start_frame,t_p)]);
      if(inter <= 0) continue;
      if(c_prev != y->bouts[beh][j].behavior) {
        loss_score += inter;
      }
      loss_fp = loss_fn = loss_score/2;
    }
#else
    double dur = b->frame_times[my_min(my_min(y->bouts[beh][gt_bout[t]].end_frame,t),TT)]-b->frame_times[t_p];
    loss_fp = fp = match_false_positive_cost(dur, beh, c_prev); // l^b_fp in (5) 
    loss_fn = 0;
    for(int j = gt_bout[t_p]; j <= gt_bout[t] && j < y->num_bouts[beh]; j++) {
      // Loop through all ground truth bouts intersecting with the completed bout
      double inter = (b->frame_times[my_min(my_min(y->bouts[beh][j].end_frame,t),TT)] - 
                      b->frame_times[my_max(y->bouts[beh][j].start_frame,t_p)]);
      if(inter <= 0) continue;
      if(c_prev == y->bouts[beh][j].behavior) {
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
					 double loss_fp, double extreme_vals[2][NUMFEAT], double duration_score) {
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
  if(extreme_vals) {
    memcpy(state.extreme_vals[0], extreme_vals[0], sizeof(double)*NUMFEAT);
    memcpy(state.extreme_vals[1], extreme_vals[1], sizeof(double)*NUMFEAT);
  }
  state.duration_score = duration_score;
}


// Helper function for Inference() subject to partial labeling constraint
// Update transition counts (e.g., number of times a transition from a bout of behavior 
// c_p to a bout of behavior c_p) occurs in the training set using a partial label
void SVMBehaviorSequence::update_transition_counts_with_partial_label(int beh, BehaviorBoutSequence *y_partial, 
                 int** &old_class_transitions, int* &old_class_transition_counts, int* &old_class_training_counts) {
  int c, c_p, i, j;

  if(y_partial) {
    old_class_transition_counts = (int*)malloc(num_classes[beh]*(sizeof(int*)+sizeof(int)*(2+num_classes[beh])));
    old_class_training_counts = old_class_transition_counts + num_classes[beh];
    old_class_transitions = (int**)(old_class_training_counts + num_classes[beh]);
    int *ptr = (int*)(old_class_transitions + num_classes[beh]);
    for(c = 0; c < num_classes[beh]; c++) {
      old_class_transition_counts[c] = class_training_transitions_count[beh][c];
      old_class_training_counts[c] = class_training_count[beh][c];
      old_class_transitions[c] = ptr;
      for(int i = 0; i < class_training_transitions_count[beh][c]; i++)
	old_class_transitions[c][i] = class_training_transitions[beh][c][i];
      ptr += class_training_transitions_count[beh][c];
    }

    for(i = 0; i < y_partial->num_bouts[beh]; i++) {
      c = y_partial->bouts[beh][i].behavior;
      if(!class_training_count[beh][c])
        class_training_count[beh][c]++;
      if(i == 0 || y_partial->bouts[beh][i].start_frame != y_partial->bouts[beh][i-1].end_frame)
	continue;
      c_p = y_partial->bouts[beh][i-1].behavior;
      
      if(c >= 0 && c_p >= 0) {
        for(j = 0; j < class_training_transitions_count[beh][c_p]; j++)
          if(class_training_transitions[beh][c][j] == c_p)
            break;
        if(j == class_training_transitions_count[beh][c]) {
          class_training_transitions[beh][c][class_training_transitions_count[beh][c]++] = c_p;
	}
      }
    }
  }
}

int SVMBehaviorSequence::get_bout_start_time(int beh, int *durations, int &tt, int t_p, int t, int &next_duration, 
					     int &last_gt, int &last_partial, int *gt_bout, int *partial_label_bout, 
					     BehaviorBoutSequence *y, BehaviorBoutSequence *y_partial, 
					     int &restrict_c_prev, int &restrict_c_next, bool *allowable_time_frames) {
  bool isFirst = t_p == t;
  next_duration = 1;

  if(y) {
    // When given a groundtruth label y, stores the index of the bout corresponding to this timestep
    // in y->bouts[beh]
    if(isFirst) {
      last_gt = gt_bout[t] = gt_bout[t-1];
      while(y->bouts[beh][gt_bout[t]].end_frame < t)
	last_gt = gt_bout[t] = gt_bout[t]+1;
    }
    while(y && (last_gt >= y->num_bouts[beh] || t_p <= y->bouts[beh][last_gt].start_frame))
      last_gt--;
  }

  if(y_partial) {
    if(isFirst) {
      // Invoked the first time the dynamic programming algorithm gets to frame t.
      // Cache the index of the bout in the partial label that contains timestep t
      last_partial = partial_label_bout[t] = partial_label_bout[t-1] < y_partial->num_bouts[beh] && 
        y_partial->bouts[beh][partial_label_bout[t-1]].end_frame <= t ? 
        partial_label_bout[t-1]+1 : partial_label_bout[t-1];

      // If the partial label contains a bout of a particular label at frame t, then it must be the case
      // the bout we are predicting also has a class c of that behavior.  So set restrict_c_next
      if(partial_label_bout[t] < y_partial->num_bouts[beh] && 
         y_partial->bouts[beh][partial_label_bout[t]].start_frame <= t)
        restrict_c_next = y_partial->bouts[beh][partial_label_bout[t]].behavior;
    }
  }

  int t_p_before = t_p;
  if(durations)
    t_p = t-durations[tt];
  else
    t_p = t-tt-1;

  // We may choose to add extra durations to make sure we explore the solutions
  // contained in y or b->partial_label
  if(y && last_gt >= 0 && (t_p < 0 ? -1 : gt_bout[t_p]) != last_gt && t_p != y->bouts[beh][last_gt].start_frame && 
     (!allowable_time_frames || allowable_time_frames[y->bouts[beh][last_gt].start_frame])) {
    t_p = y->bouts[beh][last_gt].start_frame;
    next_duration = 0;
  }

  if(y_partial && y_partial->num_bouts[beh]) {
    if(last_partial < y_partial->num_bouts[beh] && t_p < y_partial->bouts[beh][last_partial].start_frame && 
       t_p_before > y_partial->bouts[beh][last_partial].start_frame) {
      t_p = y_partial->bouts[beh][last_partial].start_frame;
      next_duration = 0;
    } else if(last_partial < y_partial->num_bouts[beh] && t_p < y_partial->bouts[beh][last_partial].end_frame && 
       t_p_before > y_partial->bouts[beh][last_partial].end_frame) {
      t_p = y_partial->bouts[beh][last_partial].end_frame;
      next_duration = 0;
    }
    while(last_partial >= y_partial->num_bouts[beh] || (t_p < y_partial->bouts[beh][last_partial].start_frame && last_partial > 0 && t_p < y_partial->bouts[beh][last_partial-1].end_frame))
      last_partial--;
    if(last_partial < y_partial->num_bouts[beh] && t_p < y_partial->bouts[beh][last_partial].start_frame && 
       t_p_before > y_partial->bouts[beh][last_partial].start_frame) {
      t_p = y_partial->bouts[beh][last_partial].start_frame;
      next_duration = 0;
    } else if(last_partial < y_partial->num_bouts[beh] && t_p < y_partial->bouts[beh][last_partial].end_frame && 
       t_p_before > y_partial->bouts[beh][last_partial].end_frame) {
      t_p = y_partial->bouts[beh][last_partial].end_frame;
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

bool SVMBehaviorSequence::check_agreement_with_partial_label(BehaviorBoutSequence *y_partial, int beh, int t_p, int t,
							     int *partial_label_bout, int &restrict_c_prev) {
  if(y_partial) {
    // If the partial label contains a bout of a particular label at frame t_p, then it must be the case
    // the bout we are predicting also has a class c_prev of that behavior.  So set restrict_c_prev
    if(y_partial) {
      if(partial_label_bout[t_p] < y_partial->num_bouts[beh] && 
         y_partial->bouts[beh][partial_label_bout[t_p]].start_frame <= t_p) {
        // If there are multiple bouts in  partial_label between t_p and t with different class labels, 
        // then it is definitely the case that the bout we are proposing doesn't agree with the partial label
        if(restrict_c_prev >= 0 && 
           y_partial->bouts[beh][partial_label_bout[t_p]].behavior != restrict_c_prev)
          return false;
        restrict_c_prev = y_partial->bouts[beh][partial_label_bout[t_p]].behavior;
      }
    }
  }

  return true;
}

// Helper function for Inference() subject to partial labeling constraint
void SVMBehaviorSequence::restore_transition_counts(int beh, BehaviorBoutSequence *y_partial, int** &old_class_training_transitions, 
						    int* &old_class_transition_counts, int* &old_class_training_counts) {
  if(y_partial && old_class_training_transitions) {
    for(int c = 0; c < num_classes[beh]; c++) {
      for(int i = 0; i < old_class_transition_counts[c]; i++)
	class_training_transitions[beh][c][i] = old_class_training_transitions[c][i];
      class_training_transitions_count[beh][c] = old_class_transition_counts[c];
      class_training_count[beh][c] = old_class_training_counts[c];
    }
    free(old_class_transition_counts);
  }
}


// Helper function for Inference().  Extract the optimal solution after dynamic programming has run by backtracking
// through its cache tables
#if USE_DURATION_COST > 0
void SVMBehaviorSequence::backtrack_optimal_solution(BehaviorBoutSequence *ybar, int beh, double **table, 
						     BehaviorBout **states, double *duration_weights, int T, int c_prev, int t_p) {
#else
void SVMBehaviorSequence::backtrack_optimal_solution(BehaviorBoutSequence *ybar, int beh, double **table, 
						     BehaviorBout **states, int T, int c_prev, int t_p) {
#endif
  // First backtrack to count the number of bouts in the optimal solution  
  int t = T, c = 0; 
  ybar->num_bouts[beh] = 0;
#if USE_NEW_LOSS > 0
  if (c_prev > -1) {
    ybar->num_bouts[beh]++;
    t = t_p;
    c = c_prev;
  }
#endif
  while(t >= 0 && states[t][c].start_frame >= 0) { 
    int tt = states[t][c].start_frame;
    c = states[t][c].behavior;
    t = tt;
    ybar->num_bouts[beh]++; 
  }
  ybar->bouts[beh] = (BehaviorBout*)malloc(sizeof(BehaviorBout)*ybar->num_bouts[beh]);

  // Backtrack one more time to actually store that solution
  t = T; c = 0; 
  ybar->scores[beh] = 0; 
  ybar->losses[beh] = 0;
  int i = ybar->num_bouts[beh]-1;
#if USE_NEW_LOSS > 0
  if (c_prev > -1) {
    ybar->bouts[beh][i].start_frame = t_p;
    ybar->bouts[beh][i].end_frame = T;
    ybar->bouts[beh][i].behavior = c_prev;
    c = c_prev;
    t = t_p;
    i--;
  }
#endif
  while(t >= 0 && states[t][c].start_frame >= 0) {
#if USE_NEW_LOSS > 0 
   if (c_prev == -1){
#endif
    ybar->scores[beh] += states[t][c].bout_score + states[t][c].transition_score + states[t][c].unary_score + states[t][c].duration_score;

    
    ybar->losses[beh] += states[t][c].loss_fn + states[t][c].loss_fp;
#if USE_NEW_LOSS > 0
   }
#endif
    ybar->bouts[beh][i] = states[t][c];
    int tt = states[t][c].start_frame;
    c = states[t][c].behavior; 
    t = tt; 
    i--; 
  }
#if USE_NEW_LOSS > 0
 if (c_prev == -1)
#endif
  ybar->score += ybar->scores[beh];
}

void SVMBehaviorSequence::print_bout_sequence_scores(BehaviorBoutSequence *y, int beh) {
  double sum = 0;
  for(int i = 0; i < y->num_bouts[beh]; i++) {
    double score = y->bouts[beh][i].bout_score + y->bouts[beh][i].unary_score + y->bouts[beh][i].transition_score + y->bouts[beh][i].loss_fn +
      y->bouts[beh][i].loss_fp;
    sum += score;
    fprintf(stderr, "(%d %d), b=%d, sum=%lf, score=%lf, bout_score=%lf, unary_score=%lf, transition_score=%lf, loss_fn=%lf, loss_fp=%lf\n", 
	    y->bouts[beh][i].start_frame, y->bouts[beh][i].end_frame, y->bouts[beh][i].behavior, sum, score,
	    y->bouts[beh][i].bout_score, y->bouts[beh][i].unary_score, y->bouts[beh][i].transition_score, 
	    y->bouts[beh][i].loss_fn, y->bouts[beh][i].loss_fp);  
  }
}


 double SVMBehaviorSequence::get_duration_score(int duration, int c, int beh, double *duration_weights) {
   double duration_score = 0;
#if USE_DURATION_COST > 0
    double duration_diff = 0;
#if USE_DURATION_COST > 1
    if (duration < min_frame_duration[beh][c]) duration_diff = 10;
    else if (duration > max_frame_duration[beh][c]) duration_diff = 10;
#else
    if (duration < min_frame_duration[beh][c]) 
	duration_diff = min_frame_duration[beh][c] - duration;
    else if (duration > max_frame_duration[beh][c]) 
	duration_diff = duration - max_frame_duration[beh][c];
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
void SVMBehaviorSequence::sanity_check_dynamic_programming_solution(int beh, BehaviorBoutFeatures *b, BehaviorBoutSequence *ybar, 
				   BehaviorBoutSequence *y, SparseVector *w, double **class_weights, double **transition_weights, 
								    double *unary_weights, double *duration_weights, double **table, BehaviorBout **states, int T, BehaviorBoutSequence *y_partial) {
#else
void SVMBehaviorSequence::sanity_check_dynamic_programming_solution(int beh, BehaviorBoutFeatures *b, BehaviorBoutSequence *ybar, 
				   BehaviorBoutSequence *y, SparseVector *w, double **class_weights, double **transition_weights, 
				   double *unary_weights, double **table, BehaviorBout **states, int T, BehaviorBoutSequence *y_partial) {
#endif
  double *tmp_features = (double*)malloc((num_features+1)*sizeof(double));

  double score = 0;
  double loss = 0;
  for(int i = 0; i < ybar->num_bouts[beh]; i++) {
    ybar->bouts[beh][i].duration_score = get_duration_score(ybar->bouts[beh][i].end_frame-ybar->bouts[beh][i].start_frame, ybar->bouts[beh][i].behavior, beh, duration_weights);
    score += ybar->bouts[beh][i].bout_score + ybar->bouts[beh][i].transition_score + ybar->bouts[beh][i].unary_score + ybar->bouts[beh][i].duration_score;
    loss += ybar->bouts[beh][i].loss_fn + ybar->bouts[beh][i].loss_fp;
    
    // If this check fails, there is a problem with the basic dynamic programming algorithm
    assert(my_abs(score + loss - table[ybar->bouts[beh][i].end_frame][i < ybar->num_bouts[beh]-1 ? ybar->bouts[beh][i+1].behavior : 0]) <= .01); 
  }

  // Possibly redundant check: if this check fails, there is a problem with the basic dynamic programming algorithm
  if (my_abs(ybar->scores[beh] + ybar->losses[beh] - table[T][0]) > .01)
	int a = 1;
  assert(my_abs(ybar->scores[beh] + ybar->losses[beh] - table[T][0]) <= .01); 


  // Make sure the score of ybar as accumulated over the dynamic programming algorithm is the same as the score
  // that is computed as the dot product between w and Psi(ybar).  If this fails, something is probably wrong with
  // the function Psi() or Inference(), such that they are inconsistent with each other
  SparseVector ybar_psi = Psi(b, ybar);
  double real_score = w->dot(ybar_psi);
  assert(my_abs(ybar->score - real_score) < .1);

  // y is the ground truth label when Inference() is invoked in learning mode
  if(y) {
    ybar->slack += ybar->scores[beh] + ybar->losses[beh];
    ybar->loss += ybar->losses[beh];

    // Makes sure the computed components of the loss aggregated during dynamic programming are identical
    // to the loss when comparing the sequences y and ybar.  If this fails, something is probably wrong
    // with the computation of the loss during dynamic programming or with the function loss2()
    double l = loss2(y, ybar, beh, 1);
    assert(my_abs(ybar->losses[beh] - l) < .01);


    // Compute the scores for the ground truth label y.  The cache tables in dynamic programming should always have
    // at least as high a score as the score yielded by the ground truth labeling.  If this fails, something
    // is probably wrong with the dynamic programming algorithm, such that it hasn't included the ground truth
    // label as a possible segmentation
    y->scores[beh] = 0;
    for(int i = 0; i < y->num_bouts[beh]; i++) {
      psi_bout(b, y->bouts[beh][i].start_frame, y->bouts[beh][i].end_frame, beh, -1, tmp_features, true, false);  
      y->bouts[beh][i].bout_score = 0;
      for(int k = 0; k < num_features; k++) 
        y->bouts[beh][i].bout_score += class_weights[y->bouts[beh][i].behavior][k]*tmp_features[k];

      if(i < y->num_bouts[beh]-1)
        y->bouts[beh][i].transition_score = transition_weights[y->bouts[beh][i].behavior][y->bouts[beh][i+1].behavior];
      else
	y->bouts[beh][i].transition_score = 0;
      y->bouts[beh][i].unary_score = unary_weights[y->bouts[beh][i].behavior];

      y->bouts[beh][i].loss_fn = y->bouts[beh][i].loss_fp = 0;
      y->bouts[beh][i].duration_score = get_duration_score(y->bouts[beh][i].end_frame-y->bouts[beh][i].start_frame, y->bouts[beh][i].behavior, beh, duration_weights);
      y->scores[beh] += y->bouts[beh][i].bout_score + y->bouts[beh][i].transition_score + y->bouts[beh][i].unary_score + y->bouts[beh][i].duration_score;

      // Making sure that ybar_score + ybar_loss >= y_score, so far
      if(!y->disable_checks && y->scores[beh] > .01+(i < y->num_bouts[beh]-1 ? table[y->bouts[beh][i+1].start_frame][y->bouts[beh][i+1].behavior] : table[T][0]) && !y_partial ) {
        // Something went wrong, it might be informative for debugging (break here using gdb) to test the same dynamic 
        // programming problem but without using loss, then compare y_max to y
        g_table = table; g_states = states; g_y = y;
        assert(0);
      }
    }
    y->score += y->scores[beh];
    if(!y->disable_checks)
      assert(ybar->score+ybar->loss+.01 >= y->score);

    ybar->slack -= y->score;
    
    // Probably redundant with earlier checks, if this fails it means the dynamic programming algorithm hasn't 
    // correctly included the ground truth label as a possible segmentation
    if(!y->disable_checks)
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
    for(int beh = 0; beh < behaviors->num; beh++) {
      // Ensure that the ground truth segmentation y_gt is included in the search space for Inference()
      if(y_gt) 
	for(int i = 0; i < y_gt->num_bouts[beh]; i++) 
	  allowable_time_frames[y_gt->bouts[beh][i].start_frame] = allowable_time_frames[y_gt->bouts[beh][i].end_frame] = true;
      if(y_partial) {
	for(int i = 0; i < y_partial->num_bouts[beh]; i++) {
	  // Ensure that the partial label is included in the search space for Inference()
	  allowable_time_frames[y_partial->bouts[beh][i].start_frame] = allowable_time_frames[y_partial->bouts[beh][i].end_frame] = true;

	}
      }
    }

    if(y_partial) {
      for(int beh = 0; beh < behaviors->num; beh++) {
	for(int i = 0; i < y_partial->num_bouts[beh]; i++) {
	  // Never allow transitions within a bout in the partial label
	  for(int j = y_partial->bouts[beh][i].start_frame+1; j < y_partial->bouts[beh][i].end_frame; j++) 
	    allowable_tmp[j] = allowable_time_frames[j] = false;
	}
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
  for(int ind = 0; ind < behaviors->num; ind++) {
    for(int i = 1; i < num_bouts[ind]; i++) {
      assert(bouts[ind][i].start_frame >= bouts[ind][i-1].end_frame);
    }

    int num = 1;
    for(int i = 1; i < num_bouts[ind]; i++) {
      if(bouts[ind][i].behavior == bouts[ind][num-1].behavior && 
	 bouts[ind][i].start_frame == bouts[ind][num-1].end_frame)
	bouts[ind][num-1].end_frame = bouts[ind][i].end_frame;
      else
	bouts[ind][num++] = bouts[ind][i];
    }
    num_bouts[ind] = my_min(num_bouts[ind], num);
    for(int i = 1; i < num_bouts[ind]; i++) {
      assert(bouts[ind][i].start_frame >= bouts[ind][i-1].end_frame);
    }
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
    for(int ind = 0; ind < behaviors->num && !has_unlabeled; ind++) {
      int lastframe = 0;
      for(int i = 0; i < y_gt->num_bouts[ind]; i++) {
	if(y_gt->bouts[ind][i].start_frame > lastframe) {
	  has_unlabeled = true;
	  break;
	}
	lastframe = y_gt->bouts[ind][i].end_frame;
      }
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
    for(int ind = 0; ind < behaviors->num; ind++) {
      int num_added = 0, lastframe = 0, partial_ind = 0, num_partial = 0;
      int s=0, e=0, b=0;
      for(int i = 0; i <= y_gt->num_bouts[ind]; i++) {
	if(i == y_gt->num_bouts[ind]) {
	  if(e >= x->num_frames)
	    break;
	  s = e = x->num_frames;
	  b = 0;
	} else {
	  s = y_gt->bouts[ind][i].start_frame;
	  e = y_gt->bouts[ind][i].end_frame;
	  b = y_gt->bouts[ind][i].behavior;
	}

	if(y_partial && y_partial_new) {
	  num_partial = y_partial->num_bouts[ind];
	  while(partial_ind < y_partial->num_bouts[ind] && y_partial->bouts[ind][partial_ind].end_frame <= s) {
	    y_partial_new->bouts[ind] = (BehaviorBout*)realloc(y_partial_new->bouts[ind], sizeof(BehaviorBout)*(num_partial+num_added+1));
	    y_partial_new->bouts[ind][y_partial_new->num_bouts[ind]++] = y_partial->bouts[ind][partial_ind++];
	  }
	}
	if(s > lastframe || i == y_gt->num_bouts[ind]) {
	  // Pad unlabelled bouts in y_gt as the "None" class
	  y_gt_new->bouts[ind] = (BehaviorBout*)realloc(y_gt_new->bouts[ind], sizeof(BehaviorBout)*(y_gt->num_bouts[ind]+num_added+1));
	  y_gt_new->bouts[ind][y_gt_new->num_bouts[ind]].start_frame = lastframe;
	  y_gt_new->bouts[ind][y_gt_new->num_bouts[ind]].end_frame = s;
	  y_gt_new->bouts[ind][y_gt_new->num_bouts[ind]++].behavior = 0;

	  if(y_partial_new) {
	    // Set the corresponding regions in y_partial to also have the "None" class
	    y_partial_new->bouts[ind] = (BehaviorBout*)realloc(y_partial_new->bouts[ind], sizeof(BehaviorBout)*(num_partial+num_added+1));
	    y_partial_new->bouts[ind][y_partial_new->num_bouts[ind]].start_frame = lastframe;
	    y_partial_new->bouts[ind][y_partial_new->num_bouts[ind]].end_frame = s;
	    y_partial_new->bouts[ind][y_partial_new->num_bouts[ind]++].behavior = 0;
	  }
	  num_added++;
	}
	if(i != y_gt->num_bouts[ind]) {
	  y_gt_new->bouts[ind][y_gt_new->num_bouts[ind]].start_frame = s;
	  y_gt_new->bouts[ind][y_gt_new->num_bouts[ind]].end_frame = lastframe = e;
	  y_gt_new->bouts[ind][y_gt_new->num_bouts[ind]].behavior = b;
	  y_gt_new->num_bouts[ind]++;
	}
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
// the first class_weightsXnum_features model parameters correspond to the
// class-feature weights, and the last num_classesXnum_classes weights correspond
// to the class transition weights
 void SVMBehaviorSequence::init_weight_pointers(int beh, double* &ptr, double** &class_weights, double** &transition_weights, double* &unary_weights, double* &duration_weights) {
   int i;
   class_weights = (double**)malloc(2*num_classes[beh]*sizeof(double*));
   transition_weights = class_weights+num_classes[beh];
   for(i = 0; i < num_classes[beh]; i++, ptr += num_features) 
     class_weights[i] = ptr;
   for(i = 0; i < num_classes[beh]; i++, ptr += num_classes[beh]) 
     transition_weights[i] = ptr; 
   unary_weights = ptr; 
   ptr += num_classes[beh];
#if USE_DURATION_COST > 0
   duration_weights = ptr;
   ptr += num_classes[beh];
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
  double *tmp_features = (double*)malloc(2*(num_features+1)*sizeof(double));
  double *tmp_features2 = tmp_features + (num_features+1);
  double *ww = w->get_non_sparse<double>(sizePsi);
  double *ptr = ww;
  int T = b->num_frames;
  //if(y && y->num_bouts && y->num_bouts[0]) 
  //T = my_min(T, y->bouts[0][y->num_bouts[0]-1].end_frame);
  
  // If the ground truth label y has unlabeled regions, the desired behavior is to ignore these 
  // regions during training.  This is handled by altering y and y_partial
  fill_unlabeled_gt_frames(y, y_partial);
  ybar->filled = true;
  if(y) assert(y->bouts[0][y->num_bouts[0]-1].end_frame == T);

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


  // Currently optimizes each behavior group independently
  for(int beh = 0; beh < behaviors->num; beh++) {
    if(behavior >= 0 && beh != behavior) {
      // If behavior != -1, it means we are only optimizing a single behavior group
      // and can skip all others
      ptr += getPsiSize(num_features, num_classes[beh]);
      continue;
    }

    // Allocate buffers for dynamic programming
    // table[t][c] will store the maximum score for any sub-solution to frames 1...t in which 
    // a bout of class c begins at time t, and states[t][c] will store the corresponding bout labels
    double *bout_scores = (double*)malloc(num_classes[beh]*sizeof(double));
    double **class_weights, **transition_weights, *unary_weights;
#if USE_DURATION_COST > 0
    double *duration_weights;
#endif
    double *fn = y ? (double*)malloc(sizeof(double)*y->num_bouts[beh]*2) : NULL;
    double *dur_gt = fn ? fn + y->num_bouts[beh] : NULL;
    double **table = (double**)malloc((T+1)*(sizeof(double*)+num_classes[beh]*sizeof(double))), *ptr3;
    BehaviorBout **states = (BehaviorBout**)malloc((T+1)*(sizeof(BehaviorBout*)+num_classes[beh]*sizeof(BehaviorBout))), *ptr2;
    int **old_class_transitions = NULL, *old_class_transition_counts = NULL, *old_class_training_counts = NULL;
    double extreme_vals[2][NUMFEAT], extreme_vals2[2][NUMFEAT];
    int i;
    for(i = 0,  ptr3 = (double*)(table+T+1), ptr2 = (BehaviorBout*)(states+T+1); i <= T; 
        i++, ptr3 += num_classes[beh], ptr2 += num_classes[beh]) {
      table[i] = ptr3;
      states[i] = ptr2;
    }

#if USE_DURATION_COST > 0
    init_weight_pointers(beh, ptr, class_weights, transition_weights, unary_weights, duration_weights);
#else
    init_weight_pointers(beh, ptr, class_weights, transition_weights, unary_weights, NULL);
#endif

    // When given a user supplied label, it may be the case that some class labels or label 
    // sequence in the partial labelling never appeared in the training set.  We add a couple of
    // checks here to protect against this case (otherwise our algorithm would find no solution)
    update_transition_counts_with_partial_label(beh, y_partial, old_class_transitions, 
						old_class_transition_counts, old_class_training_counts);

    // If evaluating loss with respect to a ground truth label y, compute the maximum false negative cost,
    // if every bout in y was missed entirely.  This loss will be subtracted off with each iteration of 
    // dynamic programming
    double max_fn_cost = 0;
    if(y) {
      int TT = my_min(T,b->num_frames-1);  
      for(i = 0; i < y->num_bouts[beh]; i++) {
        dur_gt[i] = (b->frame_times[my_min(y->bouts[beh][i].end_frame,TT)] -
                     b->frame_times[y->bouts[beh][i].start_frame]);
        fn[i] = match_false_negative_cost(dur_gt[i], beh, y->bouts[beh][i].behavior);
        max_fn_cost += fn[i];
      }
    }

    // Base case: initialize scores to 0
    for(int c = 0; c < num_classes[beh]; c++) {
      table[0][c] = 0;
      states[0][c].start_frame = states[0][c].end_frame = states[0][c].behavior = -1;
      for(int i=0; i<NUMFEAT; i++){
	states[0][c].extreme_vals[0][i] = -INFINITY;
	states[0][c].extreme_vals[1][i] = INFINITY;
      }
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

      for(int c = 0; c < num_classes[beh]; c++) 
        table[t][c] = -INFINITY;
    
      // When given a manually supplied partial labelling, store the index of the bout corresponding to 
      // this timestep in b->partial_label->bouts[beh]
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
	
	t_p = get_bout_start_time(beh, durations, tt, t_p, t, next_duration, last_gt, last_partial, gt_bout, 
				  partial_label_bout, y, y_partial, restrict_c_prev, restrict_c_next, allowable_time_frames);

        // We can quickly discard all solutions where a candidate bout (t_p,t) overlaps a region in the partial
        // label in which there are multiple bouts of different classes.  Furthermore, we can compute
        // restrict_c_prev and restrict_c_next, which are restrictions on the possible labels of
        // c_prev and c_next, respectively
        if(!check_agreement_with_partial_label(y_partial, beh, t_p, t, partial_label_bout, restrict_c_prev))
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
	psi_bout(b, t_p, t, beh, -1, tmp_features, true, !is_first, extreme_vals); 
        for(int c_prev = 0; c_prev < num_classes[beh]; c_prev++) {
          bout_scores[c_prev] = 0;
	  if(!NONE_CLASS_HAS_NO_SCORE || c_prev != 0)
	    for(int k = 0; k < num_features; k++) 
	      bout_scores[c_prev] += class_weights[c_prev][k]*tmp_features[k];
	}
        
	// Iterate through all classes c_next, where we are transitioning from a bout of class c_prev to a 
	// bout of class c_next, which is beginning at frame t.  
        for(int c_next = 0; c_next < (t == T ? 1 : num_classes[beh]); c_next++) {
	  if(t < T && (!class_training_count[beh][c_next] ||   // class c_next doesn't appear in the training set
	     (restrict_c_next >= 0 && c_next != restrict_c_next)))  // class c_next disagrees with the partial label
	    continue; 

          // Iterate through all classes c_prev, where we are transitioning from a bout of class c_prev to a
          // bout of class c_next, which is beginning at frame t.  Don't even consider class transition pairs from
          // c_prev to c_next that never occur in the training set
          for(int c_ind = 0; c_ind < (t == T ? num_classes[beh] : class_training_transitions_count[beh][c_next]); c_ind++) { 
            int c_prev = (t == T ? c_ind : class_training_transitions[beh][c_next][c_ind]); 
            if(!class_training_count[beh][c_prev] ||     // class c_prev doesn't appear in the training set
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
            double loss_score = compute_updated_bout_loss(b, y, beh, T, t_p, t, c_prev, fn, gt_bout, dur_gt, loss_fp, loss_fn);
#if USE_NEW_LOSS > 0
	    // compute loss differently
	    //BehaviorBoutSequence *ytemp = (BehaviorBoutSequence*)NewStructuredLabel(b);
            backtrack_optimal_solution(ybar, beh, table, states, duration_weights, t, c_prev, t_p); 
            loss_score = Loss(y, ybar);
#endif
	    double duration_score = get_duration_score(t-t_p, c_prev, beh, duration_weights);
            double score = bout_score + transition_score + unary_score + loss_score + duration_score;       

            // Check if the completed bout has a higher score than all previously examined solutions that 
            // begin a bout of class c at time t
            double f = table[t_p][c_prev] + score;
            if isnan(f)
                int a = 0;
            assert(!isnan(f));
            if(f > table[t][c_next]) {
              table[t][c_next] = f;
              store_solution(states[t][c_next], t_p, t, c_prev, bout_score, transition_score, unary_score, loss_fn, loss_fp, extreme_vals, duration_score);
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
	      for(int i=0; i<num_base_features; i++) {
		// Compute the updated min/max of each bout feature as the min/max over precomputed min/max values
		// for the regions (t_p_p,t_p) and (t_p,t)
		extreme_vals2[0][i] = my_max(states[t_p][c_next].extreme_vals[0][i],b->bout_max_feature_responses[i]);
		extreme_vals2[1][i] = my_min(states[t_p][c_next].extreme_vals[1][i],b->bout_min_feature_responses[i]);
	      }
 
              double bout_score = 0; 
	      if(!NONE_CLASS_HAS_NO_SCORE || c_prev != 0) {
		psi_bout(b, t_p_p, t, beh, -1, tmp_features2, true, !is_first, NULL, extreme_vals2);
		for(int k = 0; k < num_features; k++) 
		  bout_score += class_weights[c_prev][k]*tmp_features2[k];
	      }
	      double loss_fn=0, loss_fp=0;
              double transition_score = t != T ? transition_weights[c_prev][c_next] : 0;   
	      double unary_score = unary_weights[c_prev];
              double loss_score = compute_updated_bout_loss(b, y, beh, T, t_p_p, t, c_prev, fn, gt_bout, dur_gt, loss_fp, loss_fn);
#if USE_NEW_LOSS > 0
	      // compute loss differently
	      //BehaviorBoutSequence *ytemp = (BehaviorBoutSequence*)NewStructuredLabel(b);;
              backtrack_optimal_solution(ybar, beh, table, states, duration_weights, t, c_prev, t_p_p); 
              loss_score = Loss(y, ybar);
#endif
	      double duration_score = get_duration_score(t-t_p_p, c_prev, beh, duration_weights);
              double score = bout_score + transition_score + unary_score + loss_score + duration_score;
	      double f = table[t_p_p][c_prev] + score;
              if(f > table[t][c_next]) {
                table[t][c_next] = f;
                store_solution(states[t][c_next], t_p_p, t, c_prev, bout_score, transition_score, unary_score, loss_fn, loss_fp, extreme_vals2, duration_score);
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
	bool is_first = true;
	int tt = 0, next_duration = 1, t_p = t;
	int last_partial, last_gt;
	while(next_duration >= 0 && t_p >= 0) {
	  // Set durations=NULL to exhaustively consider all durations
	  t_p = get_bout_start_time(beh, NULL, tt, t_p, t, next_duration, last_gt, last_partial, gt_bout, 
				    partial_label_bout, y, y_partial, restrict_c_prev, restrict_c_next, allowable_time_frames);
	
	  // We can quickly discard all solutions where a candidate bout (t_p,t) overlaps a region in the partial
	  // label in which there are multiple bouts of different classes.  Furthermore, we can compute
	  // restrict_c_prev and restrict_c_next, which are restrictions on the possible labels of
	  // c_prev and c_next, respectively
	  if(!check_agreement_with_partial_label(y_partial, beh, t_p, t, partial_label_bout, restrict_c_prev))
	    break;

	  // allowable_time_frames is a computational time saving trick, where we only allow bouts 
	  // to start or end at a subset of allowable time frames
	  if(allowable_time_frames && !allowable_time_frames[t_p])
	    continue; 
       
        
	  // Iterate through all classes c_next, where we are transitioning from a bout of class c_prev to a 
	  // bout of class c_next, which is beginning at frame t.  
	  for(int c_next = 0; c_next < (t == T ? 1 : num_classes[beh]); c_next++) {
	    if(t < T && (!class_training_count[beh][c_next] ||   // class c_next doesn't appear in the training set
			 (restrict_c_next >= 0 && c_next != restrict_c_next)))  // class c_next disagrees with the partial label
	      continue; 
	  
	    // Select c_prev=0, the none behavior
	    for(int c_ind = 0; c_ind < (t == T ? num_classes[beh] : class_training_transitions_count[beh][c_next]); c_ind++) {
	      int c_prev = (t == T ? c_ind : class_training_transitions[beh][c_next][c_ind]); 
	      if(c_prev != 0 || !class_training_count[beh][c_prev] ||     // class c_prev doesn't appear in the training set
		 (restrict_c_prev >= 0 && c_prev != restrict_c_prev))  // class c_prev disagrees with the partial label
		continue; 
	    

	      double loss_fn=0, loss_fp=0;
	      double bout_score = 0;
	      double transition_score = t != T ? transition_weights[c_prev][c_next] : 0;   
	      double unary_score = unary_weights[c_prev];
	      double loss_score = compute_updated_bout_loss(b, y, beh, T, t_p, t, c_prev, fn, gt_bout, dur_gt, loss_fp, loss_fn);
#if USE_NEW_LOSS > 0
	      // compute loss differently
	      //BehaviorBoutSequence *ytemp = (BehaviorBoutSequence*)NewStructuredLabel(b);
	      backtrack_optimal_solution(ybar, beh, table, states, duration_weights, t, c_prev, t_p); 
	      loss_score = Loss(y, ybar);
#endif
	      double duration_score = get_duration_score(t-t_p, c_prev, beh, duration_weights);
	      double score = bout_score + transition_score + unary_score + loss_score + duration_score;       
	    
	      // Check if the completed bout has a higher score than all previously examined solutions that 
	      // begin a bout of class c at time t
	      double f = table[t_p][c_prev] + score;
	      assert(!isnan(f));
	      if(f > table[t][c_next]) {
		table[t][c_next] = f;
		store_solution(states[t][c_next], t_p, t, c_prev, bout_score, transition_score, unary_score, loss_fn, loss_fp, extreme_vals, duration_score);
		//fprintf(stderr, "(t_p=%d, t=%d, c_prev=%d, c_next=%d), f=%f\n", t_p, t, c_prev, c_next, (float)f); 
	      } 
	    }
          } // for(int c_ind = 0; ...),  c_next =...
	}
      }


    } // for(int t = 0; ...)


#if USE_DURATION_COST > 0
    // Backtrack through table and states to extract the optimal solution
    backtrack_optimal_solution(ybar, beh, table, states, duration_weights, T);
    sanity_check_dynamic_programming_solution(beh, b, ybar, y, w, class_weights, transition_weights, unary_weights, duration_weights, table, states, T, y_partial);
#else
    // Backtrack through table and states to extract the optimal solution
    backtrack_optimal_solution(ybar, beh, table, states, T);
    sanity_check_dynamic_programming_solution(beh, b, ybar, y, w, class_weights, transition_weights, unary_weights, table, states, T, y_partial);
#endif

    // Restore modified transition tables, if necessary
    restore_transition_counts(beh, y_partial, old_class_transitions, old_class_transition_counts, old_class_training_counts);
      
    // Cleanup
    free(table);
    free(states);
    free(class_weights);
    free(bout_scores);
    if(fn) free(fn);

  } // for(int beh = 0; ...)

  // Cleanup
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

  for(int beh = 0; beh < behaviors->num; beh++) {
    y_dst->num_bouts[beh] = 0;
    y_dst->bouts[beh] = (BehaviorBout*)realloc(y_dst->bouts[beh], sizeof(BehaviorBout)*(y_src->num_bouts[beh]+2));
    for(int i = 0; i < y_src->num_bouts[beh]; i++) {
      if(y_src->bouts[beh][i].start_frame < t_start || y_src->bouts[beh][i].end_frame > t_end) {
	y_dst->bouts[beh][y_dst->num_bouts[beh]] = y_src->bouts[beh][i];
	if(y_dst->bouts[beh][y_dst->num_bouts[beh]].start_frame < t_start && 
	   y_dst->bouts[beh][y_dst->num_bouts[beh]].end_frame > t_start) {
	  if(y_dst->bouts[beh][y_dst->num_bouts[beh]].end_frame > t_end) {
	    y_dst->bouts[beh][y_dst->num_bouts[beh]++].end_frame = t_start;
	    y_dst->bouts[beh][y_dst->num_bouts[beh]] = y_src->bouts[beh][i];
	    y_dst->bouts[beh][y_dst->num_bouts[beh]].start_frame = t_end;
	  } else
	    y_dst->bouts[beh][y_dst->num_bouts[beh]].end_frame = t_start;
	  changed = true;
	} else if(y_dst->bouts[beh][y_dst->num_bouts[beh]].end_frame > t_end &&
		  y_dst->bouts[beh][y_dst->num_bouts[beh]].start_frame < t_end) {
	  y_dst->bouts[beh][y_dst->num_bouts[beh]].start_frame = t_end;
	  changed = true;
	}
	y_dst->num_bouts[beh]++;
      } else
	changed = true;
    }
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
    while(t_start < T) {
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
	SVM_cached_sample_set_add_sample(set, ybar);
	SVM_cached_sample_set_compute_features(set, trainset->examples[set->i]);
      
	// Use the current sample set to update the weights immediately.  Usually we would never do this, but
	// this allows us to update the model weights more frequently when bout sequences are really
	// long and inference is slow
	Lock();
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
  set->score_gt = w_curr->dot(*set->psi_gt);
  StructuredLabel *ybar = NewStructuredLabel(x);
  retval = Inference(x, ybar, w_curr, NULL, y_gt, 1);
  SVM_cached_sample_set_add_sample(set, ybar);
  SVM_cached_sample_set_compute_features(set, trainset->examples[set->i]);
  if(set->num_samples) {
    SVM_cached_sample s = set->samples[0];
    set->samples[0] = set->samples[set->num_samples-1];
    set->samples[set->num_samples-1] = s;
  }
  delete w_curr;

  return retval;
}



double      SVMBehaviorSequence::Loss(StructuredLabel *y_gt,  StructuredLabel *y_pred) {
	int i;
	double l = 0;
	BehaviorBoutSequence *y = (BehaviorBoutSequence*)y_gt, *y_partial=NULL;
	
	fill_unlabeled_gt_frames(y, y_partial);

	if(behavior >= 0) 
	  return loss2(y, y_pred, behavior, 0);
	else {
		for(i = 0; i < behaviors->num; i++){
			l += loss2(y, y_pred, i, 0);
		}
	}
	
	if(y != y_gt) delete y;
	if(y_partial) delete y_partial;

	return l;
}

double      SVMBehaviorSequence::loss2(StructuredLabel *y_gt,  StructuredLabel *y_pred, int beh, int debug) {
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
	if(y && y->num_bouts && y->num_bouts[0]) 
	  T = my_min(T, y->bouts[0][y->num_bouts[0]-1].end_frame);
	int TT = my_min(T,b->num_frames-1);  
	double cl, l_fn = 0, l_fn2 = 0;

#if HAMMING_LOSS == 1
	while((curr_y < y->num_bouts[beh] || curr_ybar < ybar->num_bouts[beh])) {
		if(curr_y < y->num_bouts[beh] && curr_ybar < ybar->num_bouts[beh]) {
			// Check if the current bout in y and ybar match.  
			dur_ybar = (b->frame_times[my_min(ybar->bouts[beh][curr_ybar].end_frame,TT)] - 
				b->frame_times[ybar->bouts[beh][curr_ybar].start_frame]);
			dur_y = (b->frame_times[my_min(y->bouts[beh][curr_y].end_frame,TT)] -
				b->frame_times[y->bouts[beh][curr_y].start_frame]);
			inter = (b->frame_times[my_min(my_min(y->bouts[beh][curr_y].end_frame, ybar->bouts[beh][curr_ybar].end_frame),TT)] - 
				b->frame_times[my_max(y->bouts[beh][curr_y].start_frame,ybar->bouts[beh][curr_ybar].start_frame)]);
			assert(inter >= 0);
			if(y->bouts[beh][curr_y].behavior != ybar->bouts[beh][curr_ybar].behavior) 
				l_fn += inter;
		} else
			inter = 0;

		if(curr_ybar == ybar->num_bouts[beh] || (curr_y < y->num_bouts[beh] && 
			(y->bouts[beh][curr_y].end_frame <= ybar->bouts[beh][curr_ybar].end_frame))) {
				// Go to the next bout in y adding the appropriate loss based on whether or not
				// the bout was matched to a bout in ybar
				curr_y++;
		} else {
			// Go to the next bout in ybar, adding the appropriate loss based on whether or not
			// the bout was matched to a bout in y
			l += l_fn;
			ybar->bouts[beh][curr_ybar].loss_fn = l_fn/2;
			ybar->bouts[beh][curr_ybar].loss_fp = l_fn/2;
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
	int beh_counts_gt[num_classes[beh]], beh_counts_pred[num_classes[beh]];
	for (int i=0; i<num_classes[beh]; i++) {
		beh_counts_gt[i] = 0;
		beh_counts_pred[i] = 0;
	}
	for (int i=0; i<y->num_bouts[beh]; i++) {
		beh_counts_gt[y->bouts[beh][i].behavior] ++;}
	for (int i=0; i<ybar->num_bouts[beh]; i++) {
		beh_counts_pred[ybar->bouts[beh][i].behavior] ++;}
	for (int i=0; i<num_classes[beh]; i++) {
		beh_counts_gt[i] = my_max(beh_counts_gt[i],1);
		beh_counts_pred[i] = my_max(beh_counts_pred[i],1);
	}
#endif

	while((curr_y < y->num_bouts[beh] || curr_ybar < ybar->num_bouts[beh])) {
		if (curr_ybar == ybar->num_bouts[beh]-2)
			int a = 1;
		if(curr_y < y->num_bouts[beh] && curr_ybar < ybar->num_bouts[beh]) {
			// Check if the current bout in y and ybar match.  
			dur_ybar = (b->frame_times[my_min(ybar->bouts[beh][curr_ybar].end_frame,TT)] - 
				b->frame_times[ybar->bouts[beh][curr_ybar].start_frame]);
			dur_y = (b->frame_times[my_min(y->bouts[beh][curr_y].end_frame,TT)] -
				b->frame_times[y->bouts[beh][curr_y].start_frame]);
			inter = (b->frame_times[my_min(my_min(y->bouts[beh][curr_y].end_frame, ybar->bouts[beh][curr_ybar].end_frame),TT)] - 
				b->frame_times[my_max(y->bouts[beh][curr_y].start_frame,ybar->bouts[beh][curr_ybar].start_frame)]);
			assert(inter >= 0);
			if(y->bouts[beh][curr_y].behavior == ybar->bouts[beh][curr_ybar].behavior) {
				sum_ybar += inter;
				sum_y += inter;
				//if(debug)
#if USE_NEW_LOSS > 0 
				l_fn -= match_false_negative_cost(dur_y, beh, y->bouts[beh][curr_y].behavior)*(dur_y ? (inter/dur_y) : 0) / beh_counts_gt[y->bouts[beh][curr_y].behavior];
#else 
				l_fn -= match_false_negative_cost(dur_y, beh, y->bouts[beh][curr_y].behavior)*(dur_y ? (inter/dur_y) : 0);
# endif
			} else
#if USE_NEW_LOSS > 0
				l_fn2 += match_false_negative_cost(dur_y, beh, y->bouts[beh][curr_y].behavior)*(dur_y ? (inter/dur_y) : 0) / beh_counts_gt[y->bouts[beh][curr_y].behavior];
#else
				l_fn2 += match_false_negative_cost(dur_y, beh, y->bouts[beh][curr_y].behavior)*(dur_y ? (inter/dur_y) : 0);
#endif
		} else
			inter = 0;

		if(curr_ybar == ybar->num_bouts[beh] || (curr_y < y->num_bouts[beh] && 
			(y->bouts[beh][curr_y].end_frame <= ybar->bouts[beh][curr_ybar].end_frame))) {
				// Go to the next bout in y adding the appropriate loss based on whether or not
				// the bout was matched to a bout in ybar
				cl = match_false_negative_cost(dur_y, beh, y->bouts[beh][curr_y].behavior)*(dur_y ? ((dur_y-sum_y)/dur_y) : 1);
#if USE_NEW_LOSS > 0
				cl /= beh_counts_gt[y->bouts[beh][curr_y].behavior];
#endif
				l += cl; 
				if(debug == 2) fprintf(stderr, "y[%d]->%f %f\n", curr_y, cl, l);
				curr_y++;
				sum_y = 0;
		} else {
			// Go to the next bout in ybar, adding the appropriate loss based on whether or not
			// the bout was matched to a bout in y
			cl = match_false_positive_cost(dur_ybar, beh, ybar->bouts[beh][curr_ybar].behavior) *
				(dur_ybar ? ((dur_ybar-sum_ybar)/dur_ybar) : 1); 
#if USE_NEW_LOSS > 0
			cl /= beh_counts_pred[ybar->bouts[beh][curr_ybar].behavior];
#endif
			l += cl;
			if(debug == 2) 
				fprintf(stderr, "ybar[%d]->%f %f\n", curr_ybar, cl, l);
			if(debug) {
				assert(my_abs(cl - ybar->bouts[beh][curr_ybar].loss_fp) < .00001);
				if (my_abs(l_fn2 - ybar->bouts[beh][curr_ybar].loss_fn) > .00001)
					int a = 1;
				assert(my_abs(l_fn2 - ybar->bouts[beh][curr_ybar].loss_fn) < .00001);
			}
			ybar->bouts[beh][curr_ybar].loss_fn = l_fn2;
			ybar->bouts[beh][curr_ybar].loss_fp = cl;
			curr_ybar++;
			sum_ybar = 0;
			l_fn = 0; l_fn2 = 0;
		}
	}
#endif //HAMMING_LOSS == 1
	if(debug)
		assert(my_abs(l-ybar->losses[beh]) < .01);

	return l;
}




Json::Value SVMBehaviorSequence::Save() {
	Json::Value root, num_cl, trans;
	int i, j, beh;
	root["version"] = INST_VERSION;

	for(i = 0; i < behaviors->num; i++)
	  num_cl[i] = this->num_classes[i];
	root["num_classes"] = num_cl;
	root["num_features"] = num_features;
	root["num_base_features"] = num_base_features;
	root["behavior"] = behavior;
	root["feature_diff_frames"] = feature_diff_frames;


	Json::Value tr;
	for(beh = 0; beh < behaviors->num; beh++) {
		Json::Value t;
		for(i = 0; i < num_classes[beh]; i++) {
			Json::Value o, tt;
			o["count"] = class_training_count[beh][i];
			o["transitions_count"] = class_training_transitions_count[beh][i];
			for(j = 0; j < class_training_transitions_count[beh][i]; j++) {
			  tt[j] = class_training_transitions[beh][i][j];
			}
			o["transitions"] = tt;
#if USE_DURATION_COST > 0
			Json::Value limits;
			int idx = 0;
			limits[idx] = (int)min_frame_duration[beh][i];
			idx = 1;
			limits[idx] = (int)max_frame_duration[beh][i];
			o["limits"] = limits;
#endif
			t[i] = o;
		}
		tr[beh] = t;
	}
	root["transitions"] = tr;

	Json::Value params;
	for(i = 0; i < num_base_features; i++) {
		SVMFeatureParams *p = &feature_params[i];
		Json::Value c;
		c["name"] = this->base_feature_names[i];
		c["feature_sample_smoothness_window"] = p->feature_sample_smoothness_window;
		c["num_temporal_levels"] = p->num_temporal_levels;
		c["num_bout_max_thresholds"] = p->num_bout_max_thresholds;
		c["num_bout_min_thresholds"] = p->num_bout_min_thresholds;
		c["num_bout_change_points"] = p->num_bout_change_points;
		c["num_histogram_bins"] = p->num_histogram_bins;
		c["num_histogram_temporal_levels"] = p->num_histogram_temporal_levels;
		c["num_difference_temporal_levels"] = p->num_difference_temporal_levels;
		c["num_harmonic_features"] = p->num_harmonic_features;
		c["use_bout_sum_features"] = p->use_bout_sum_features;
		c["use_bout_ave_features"] = p->use_bout_ave_features;
		c["use_bout_sum_absolute_features"] = p->use_bout_sum_absolute_features;
		c["use_bout_ave_absolute_features"] = p->use_bout_ave_absolute_features;
		c["use_standard_deviation"] = p->use_standard_deviation;
		c["use_sum_variance"] = p->use_sum_variance;
		c["use_bout_max_feature"] = p->use_bout_max_feature;
		c["use_bout_min_feature"] = p->use_bout_min_feature;
		c["use_global_difference_max_ave_features"] = p->use_global_difference_max_ave_features;
		c["use_global_difference_min_ave_features"] = p->use_global_difference_min_ave_features;
		c["use_global_difference_ave_ave_features"] = p->use_global_difference_ave_ave_features;
		c["use_global_difference_max_sum_features"] = p->use_global_difference_max_sum_features;
		c["use_global_difference_min_sum_features"] = p->use_global_difference_min_sum_features;
		c["use_global_difference_ave_sum_features"] = p->use_global_difference_ave_sum_features;
		c["use_bout_change"] = p->use_bout_change;
		c["use_bout_absolute_change"] = p->use_bout_absolute_change;
		c["use_histogram_sum_features"] = p->use_histogram_sum_features;
		c["use_histogram_ave_features"] = p->use_histogram_ave_features;
		c["use_sum_harmonic_features"] = p->use_sum_harmonic_features;
		c["use_ave_harmonic_features"] = p->use_ave_harmonic_features;
		c["use_sum_absolute_harmonic_features"] = p->use_sum_absolute_harmonic_features;
		c["use_ave_absolute_harmonic_features"] = p->use_ave_absolute_harmonic_features;
		c["use_start_sum_absolute_diff_haar_features"] = p->use_start_sum_absolute_diff_haar_features;
		c["use_end_sum_absolute_diff_haar_features"] = p->use_end_sum_absolute_diff_haar_features;
		c["use_start_sum_diff_haar_features"] = p->use_start_sum_diff_haar_features;
		c["use_end_sum_diff_haar_features"] = p->use_end_sum_diff_haar_features;
		c["use_start_ave_absolute_diff_haar_features"] = p->use_start_ave_absolute_diff_haar_features;
		c["use_end_ave_absolute_diff_haar_features"] = p->use_end_ave_absolute_diff_haar_features;
		c["use_start_ave_diff_haar_features"] = p->use_start_ave_diff_haar_features;
		c["use_end_ave_diff_haar_features"] = p->use_end_ave_diff_haar_features;
		params[i] = c;
	}
	root["feature_params"] = params;

	Json::Value mu, gamma, thresh, max_thr, min_thr;

	for(i = 0; i < num_features; i++)
		mu[i] = features_mu[i];
	root["feature_means"] = mu;
	
	for(i = 0; i < num_features; i++)
		gamma[i] = features_gamma[i];
	root["feature_gamma"] = gamma;

	for(i = 0; i < num_base_features; i++) {
		Json::Value a;
		for(j = 0; j < feature_params[i].num_histogram_bins; j++)
			a[j] = histogram_thresholds[i][j];
		thresh[i] = a;
	}
	root["feature_thresholds"] = thresh;

	for(i = 0; i < num_base_features; i++) {
		Json::Value a;
		for(j = 0; j < feature_params[i].num_bout_max_thresholds; j++)
			a[j] = max_thresholds[i][j];
		max_thr[i] = a;
	}
	root["feature_max_thresholds"] = max_thr;

	for(i = 0; i < num_base_features; i++) {
		Json::Value a;
		for(j = 0; j < feature_params[i].num_bout_min_thresholds; j++)
			a[j] = min_thresholds[i][j];
		min_thr[i] = a;
	}
	root["feature_min_thresholds"] = min_thr;

	Json::Value aa(Json::arrayValue), bb(Json::arrayValue);
	for (beh = 0; beh < behaviors->num; beh++) {
	  Json::Value a(Json::arrayValue), b(Json::arrayValue);
	  for(i=0; i<num_classes[beh]; i++){
	    a[i] = false_negative_cost[beh][i];
	    b[i] = false_positive_cost[beh][i];
	  }
	  aa[beh] = a;
	  bb[beh] = a;
	}
	root["false_negative_cost"] = aa;
	root["false_positive_cost"] = bb;
	
	

	return root;
}

bool SVMBehaviorSequence::Load(const Json::Value &root) {
	int i, j, beh;


	assert(!strcmp(root["version"].asString().c_str(), INST_VERSION));

	for(i = 0; i < behaviors->num; i++)
	  this->num_classes[i] = root["num_classes"][i].asInt();
	num_features = root["num_features"].asInt();
	num_base_features = root["num_base_features"].asInt();
	behavior = root["behavior"].asInt();
	feature_diff_frames = root["feature_diff_frames"].asInt();


	class_training_transitions = (int***)malloc(behaviors->num*sizeof(int**));
	class_training_transitions_count = (int**)malloc(behaviors->num*sizeof(int*));
	class_training_count = (int**)malloc(behaviors->num*sizeof(int*));
	for(beh = 0; beh < behaviors->num; beh++) {
		class_training_transitions[beh] = (int**)malloc(num_classes[beh]*sizeof(int*));
		class_training_transitions_count[beh] = (int*)malloc(num_classes[beh]*sizeof(int));
		class_training_count[beh] = (int*)malloc(num_classes[beh]*sizeof(int));
	}

	Json::Value tr = root["transitions"];
	for(beh = 0; beh < behaviors->num; beh++) {
		Json::Value t = tr[beh];
		for(i = 0; i < num_classes[beh]; i++) {
			class_training_transitions[beh][i] = (int*)malloc(num_classes[beh]*sizeof(int));
			Json::Value o = t[i]; 
			Json::Value tt = o["transitions"];
			class_training_transitions_count[beh][i] = o["transitions_count"].asInt();
			class_training_count[beh][i] = o["count"].asInt();
			for(j = 0; j < class_training_transitions_count[beh][i]; j++) {
			  class_training_transitions[beh][i][j] = tt[j].asInt();
			}
#if USE_DURATION_COST > 0
			Json::Value limits = o["limits"];
			int idx = 0;
			min_frame_duration[beh][i] = limits[idx].asInt();
			idx = 1;
			max_frame_duration[beh][i] = limits[idx].asInt();
#endif
		}
	}

	int num_histogram_bins = 0, num_max_thresholds = 0, num_min_thresholds = 0;
	Json::Value params;
	for(i = 0; i < num_base_features; i++) {
		SVMFeatureParams *p = &feature_params[i];
		Json::Value c = root["feature_params"][i];
		strcpy(this->base_feature_names[i], c["name"].asString().c_str());
		p->feature_sample_smoothness_window = c["feature_sample_smoothness_window"].asInt();
		p->num_temporal_levels = c["num_temporal_levels"].asInt();
		p->num_bout_max_thresholds = c["num_bout_max_thresholds"].asInt();
		p->num_bout_min_thresholds = c["num_bout_min_thresholds"].asInt();
		p->num_bout_change_points = c["num_bout_change_points"].asInt();
		p->num_histogram_bins = c["num_histogram_bins"].asInt();
		p->num_histogram_temporal_levels = c["num_histogram_temporal_levels"].asInt();
		p->num_difference_temporal_levels = c["num_difference_temporal_levels"].asInt();
		p->num_harmonic_features = c["num_harmonic_features"].asInt();
		p->use_bout_sum_features = c["use_bout_sum_features"].asBool();
		p->use_bout_ave_features = c["use_bout_ave_features"].asBool();
		p->use_bout_sum_absolute_features = c["use_bout_sum_absolute_features"].asBool();
		p->use_bout_ave_absolute_features = c["use_bout_ave_absolute_features"].asBool();
		p->use_standard_deviation = c["use_standard_deviation"].asBool();
		p->use_sum_variance = c["use_sum_variance"].asBool();
		p->use_bout_max_feature = c["use_bout_max_feature"].asBool();
		p->use_bout_min_feature = c["use_bout_min_feature"].asBool();
		p->use_global_difference_max_ave_features = c["use_global_difference_max_ave_features"].asBool();
		p->use_global_difference_min_ave_features = c["use_global_difference_min_ave_features"].asBool();
		p->use_global_difference_ave_ave_features = c["use_global_difference_ave_ave_features"].asBool();
		p->use_global_difference_max_sum_features = c["use_global_difference_max_sum_features"].asBool();
		p->use_global_difference_min_sum_features = c["use_global_difference_min_sum_features"].asBool();
		p->use_global_difference_ave_sum_features = c["use_global_difference_ave_sum_features"].asBool();
		p->use_bout_change = c["use_bout_change"].asBool();
		p->use_bout_absolute_change = c["use_bout_absolute_change"].asBool();
		p->use_histogram_sum_features = c["use_histogram_sum_features"].asBool();
		p->use_histogram_ave_features = c["use_histogram_ave_features"].asBool();
		p->use_sum_harmonic_features = c["use_sum_harmonic_features"].asBool();
		p->use_ave_harmonic_features = c["use_ave_harmonic_features"].asBool();
		p->use_sum_absolute_harmonic_features = c["use_sum_absolute_harmonic_features"].asBool();
		p->use_ave_absolute_harmonic_features = c["use_ave_absolute_harmonic_features"].asBool();
		p->use_start_sum_absolute_diff_haar_features = c["use_start_sum_absolute_diff_haar_features"].asBool();
		p->use_end_sum_absolute_diff_haar_features = c["use_end_sum_absolute_diff_haar_features"].asBool();
		p->use_start_sum_diff_haar_features = c["use_start_sum_diff_haar_features"].asBool();
		p->use_end_sum_diff_haar_features = c["use_end_sum_diff_haar_features"].asBool();
		p->use_start_ave_absolute_diff_haar_features = c["use_start_ave_absolute_diff_haar_features"].asBool();
		p->use_end_ave_absolute_diff_haar_features = c["use_end_ave_absolute_diff_haar_features"].asBool();
		p->use_start_ave_diff_haar_features = c["use_start_ave_diff_haar_features"].asBool();
		p->use_end_ave_diff_haar_features = c["use_end_ave_diff_haar_features"].asBool();
		num_histogram_bins += p->num_histogram_bins;
		num_max_thresholds += p->num_bout_max_thresholds;
		num_min_thresholds += p->num_bout_min_thresholds;
	}

	Json::Value mu = root["feature_means"], gamma = root["feature_gamma"], thresh = root["feature_thresholds"], 
	  max_thr = root["feature_max_thresholds"], min_thr = root["feature_min_thresholds"];

	features_mu = (double*)malloc(num_features*sizeof(double)*2);
	features_gamma = features_mu + num_features;
	histogram_thresholds = (double**)malloc(3*sizeof(double*)*num_base_features + 
		(num_histogram_bins+num_max_thresholds+num_min_thresholds)*sizeof(double));
	min_thresholds = histogram_thresholds + num_base_features;
	max_thresholds = min_thresholds + num_base_features;
	double *ptr2 = (double*)(max_thresholds + num_base_features);

	for(i = 0; i < num_features; i++)
	  features_mu[i] = mu[i].asDouble();
	
	for(i = 0; i < num_features; i++)
	  features_gamma[i] = gamma[i].asDouble();

	for(i = 0; i < num_base_features; i++) {
		Json::Value a = thresh[i];
		histogram_thresholds[i] = ptr2; ptr2 += feature_params[i].num_histogram_bins;
		for(j = 0; j < feature_params[i].num_histogram_bins; j++)
		  histogram_thresholds[i][j] = a[j].asDouble();
	}

	for(i = 0; i < num_base_features; i++) {
		Json::Value a = max_thr[i];
		max_thresholds[i] = ptr2; ptr2 += feature_params[i].num_histogram_bins;
		for(j = 0; j < feature_params[i].num_bout_max_thresholds; j++)
		  max_thresholds[i][j] = a[j].asDouble();
	}

	for(i = 0; i < num_base_features; i++) {
		Json::Value a = min_thr[i];
		min_thresholds[i] = ptr2; ptr2 += feature_params[i].num_histogram_bins;
		for(j = 0; j < feature_params[i].num_bout_min_thresholds; j++)
		  min_thresholds[i][j] = a[j].asDouble();
	}

	if(root["false_positive_cost"].isArray()) {
	  for (beh = 0; beh < behaviors->num; beh++) {
	    Json::Value a = root["false_negative_cost"][beh], b=root["false_positive_cost"][beh];
	    for(i=0; i<num_classes[beh]; i++){
	      false_negative_cost[beh][i] = a[i].asDouble();
	      false_positive_cost[beh][i] = b[i].asDouble();
	    }
	  }
	}

	if(!feature_names) {
	  feature_names = (char**)realloc(feature_names, sizeof(char*)*num_features);
	  compute_feature_space_size();
	  sizePsi = 0;
	  for(int bh = 0; bh < behaviors->num; bh++) 
	    sizePsi += getPsiSize(this->num_features, this->num_classes[bh]);
	}

	return true;
}



BehaviorBoutFeatures::BehaviorBoutFeatures() {
  partial_label = NULL;
  memory_buffer = NULL;
  features = NULL;
  frame_times = NULL;
  fvec = NULL;
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
	int i;
	int T = num_frames;
	if(full) {
	  long int cache_features_size = 
		num_base_features*sizeof(double*) +        // integral_features
		num_base_features*(T+1)*3*sizeof(double) + // integral_features[i]
		num_base_features*sizeof(double*) +        // integral_sqr_features
		num_base_features*(T+1)*3*sizeof(double) + // integral_sqr_features[i]
		num_base_features*sizeof(double*) +        // smoothed_features
		num_base_features*T*sizeof(double) +       // smoothed_features[i]
		3*num_base_features*sizeof(double) +       // max_feature_responses, min_feature_responses, ave_feature_responses, 
		2*num_base_features*sizeof(double) +       // bout_max_feature_responses, bout_min_feature_responses
		num_base_features*sizeof(int*) +           // histogram_bins
		num_base_features*sizeof(int)*T +          // histogram_bins[i]
		num_base_features*sizeof(double**);        // integral_histogram_features
	  for(i = 0; i < num_base_features; i++) {
		SVMFeatureParams *p = &svm->feature_params[i];
		cache_features_size +=
			sizeof(double*)*p->num_histogram_bins + //integral_histogram_features[i]
			p->num_histogram_bins*(T+1)*3*sizeof(double); //integral_histogram_features[i][j]
			//num_base_features*sizeof(double*)*p->num_histogram_bins + //integral_histogram_features[i]
			//num_base_features*p->num_histogram_bins*(T+1)*3*sizeof(double)+10000; //integral_histogram_features[i][j]
	  }
	  this->memory_buffer = (unsigned char*)malloc(cache_features_size);
	}

	if(!features) {
	  // Extract regular raw features
	  //this->features = (double**)malloc(num_base_features*(sizeof(double*)+(num_base_features+1)*T*sizeof(double)));
	  this->features = (double**)malloc((num_base_features+1)*(sizeof(double*)+(T+1)*sizeof(double)));
	  double *ptr = (double*)(this->features+num_base_features);
	  for(i = 0; i < num_base_features; i++) {
	    this->features[i] = ptr;
	    ptr += T;
	  }
	  this->frame_times = ptr;
	}

	this->partial_label = NULL;

	this->fvec = NULL;
	this->bout_start = this->bout_end = -1;
}


/*
* Precompute certain data structures like integral images, such that bout features can
* be computed more efficiently
*/
void BehaviorBoutFeatures::ComputeCaches(SVMBehaviorSequence *svm) {
	if(!svm->histogram_thresholds) return;
	if(!memory_buffer)
	  AllocateBuffers(svm);

	int i, j, k, T = num_frames;
	double f;
	unsigned char *ptr = memory_buffer;
	int window;

	// integral_features[i] store the running sum value of the i_th feature, such that the sum
	// value of that feature for an arbitrary interval [s,e) can be computed as 
	//   integral_features[e]-integral_features[s]
	// To allow the caller to index into integral_features[i] without worrying about bounds checking,
	// a buffer of size T+1 is added to the beginning and end, which behaves as if the feature in the
	// 0th timestep extends to frame -(T+1) and the (T-1)th pixel extends (T+1) frames into the future
	integral_features = (double**)ptr;
	ptr += num_base_features*sizeof(double*);
	for(i = 0; i < num_base_features; i++) {
		integral_features[i] = ((double*)ptr)+T+1;
		ptr += (T+1)*3*sizeof(double);

		integral_features[i][-T-1] = 0;
		for(j = -T-1; j < 2*(T+1)-1; j++) {
			integral_features[i][j+1] = integral_features[i][j] + 
				features[i][j < 0 ? 0 : (j >= T ? T-1 : j)];
		}
	}

	integral_sqr_features = (double**)ptr;
	ptr += num_base_features*sizeof(double*);
	for(i = 0; i < num_base_features; i++) {
		integral_sqr_features[i] = ((double*)ptr)+T+1;
		ptr += (T+1)*3*sizeof(double);

		integral_sqr_features[i][-T-1] = 0;
		for(j = -T-1; j < 2*(T+1)-1; j++) {
			integral_sqr_features[i][j+1] = integral_sqr_features[i][j] + 
				SQR(features[i][j < 0 ? 0 : (j >= T ? T-1 : j)]);
		}
	}

	// Compute smoothed versions of the raw features
	smoothed_features = (double**)ptr;
	ptr += num_base_features*sizeof(double*);
	for(i = 0; i < num_base_features; i++) {
		smoothed_features[i] = ((double*)ptr);
		ptr += T*sizeof(double);
		window = svm->feature_params[i].feature_sample_smoothness_window;
		for(j = 0; j < T; j++)
			smoothed_features[i][j] = (integral_features[i][j+window+1] - 
			integral_features[i][j-window]) / (2*window+1);
	}

	// ave_feature_responses[i], max_feature_responses[i], and min_feature_responses[i] respectively
	// store the average, maximum, and minimum response of the i_th feature over the entire sequence 1...T
	// The responses at a particular frame are smoothed within a window of size 'window', to add greater
	// robustness to a few outlier frames
	max_feature_responses = (double*)ptr;
	ptr += num_base_features*sizeof(double);
	min_feature_responses = (double*)ptr;
	ptr += num_base_features*sizeof(double);
	ave_feature_responses = (double*)ptr;
	ptr += num_base_features*sizeof(double);
	for(i = 0; i < num_base_features; i++) {
		ave_feature_responses[i] = (integral_features[i][T] -
			integral_features[i][0]) / T;
		max_feature_responses[i] = -INFINITY;
		min_feature_responses[i] = INFINITY;
		for(j = 0; j < T; j++) {
			f = smoothed_features[i][j];
			if(f < min_feature_responses[i])
				min_feature_responses[i] = f;
			if(f > max_feature_responses[i])
				max_feature_responses[i] = f;
		}
	}

	// Compute histogram features.  We discretize the space of values for the i_th feature into num_histogram bins
	// of equal size between min_feature_responses[i] and max_feature_responses[i].  histogram_bins[i][j] stores the
	// index of the assigned bin for the i_th feature at the j_th time step
	histogram_bins = (int**)ptr;
	ptr += num_base_features*sizeof(int*);
	for(i = 0; i < num_base_features; i++, ptr += sizeof(int)*T) {
		histogram_bins[i] = (int*)ptr;
		for(j = 0; j < T; j++) {
			f = features[i][j];
			for(k = 0; k < svm->feature_params[i].num_histogram_bins-1; k++) 
				if(f <= svm->histogram_thresholds[i][k]) 
					break;
			histogram_bins[i][j] = k;
		}
	}

	// Compute integral features on the histogram features (e.g. using the same method as for integral_features),
	// such that the total count of the j_th histogram bin for the i_th feature can be computed using
	//   integral_histogram_features[i][j][e]-integral_histogram_features[i][j][s]
	integral_histogram_features = (double***)ptr;
	ptr += num_base_features*sizeof(double**);
	for(i = 0; i < num_base_features; i++) {
		integral_histogram_features[i] = (double**)ptr;
		ptr += sizeof(double*)*svm->feature_params[i].num_histogram_bins;
		for(j = 0; j < svm->feature_params[i].num_histogram_bins; j++) {
			integral_histogram_features[i][j] = ((double*)ptr)+T+1;
			ptr += (T+1)*3*sizeof(double);
			integral_histogram_features[i][j][-T-1] = 0;
			for(k = -T-1; k < 2*(T+1)-1; k++) {
				integral_histogram_features[i][j][k+1] = integral_histogram_features[i][j][k] + 
					(histogram_bins[i][k < 0 ? 0 : (k >= T ? T-1 : k)] == j ? 1 : 0);
			}
		}
	}

	bout_max_feature_responses = (double*)ptr;
	ptr += num_base_features*sizeof(double);
	bout_min_feature_responses = (double*)ptr;
	ptr += num_base_features*sizeof(double);
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
  this->behaviors = svm->behaviors;
  features = x;
  filled = false;
  num_bouts = NULL;
  bouts = NULL;
  score = loss = slack = 0;
  scores = losses = NULL;
  disable_checks = false;
}

BehaviorBoutSequence::~BehaviorBoutSequence() {
	for(int i = 0; i < behaviors->num; i++)
		if(bouts && bouts[i])
			free(bouts[i]);
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
  if(num_bouts) {
    Json::Value a(Json::arrayValue);
    for(int i = 0; i < behaviors->num; i++)
      a[i] = num_bouts[i];
    r["num_bouts"] = a;
  }

  if(bouts) {
    Json::Value a(Json::arrayValue);
    for(int i = 0; i < behaviors->num; i++) {
      Json::Value b(Json::arrayValue);
      for(int j = 0; j < num_bouts[i]; j++) {
        Json::Value c;
        c["start_frame"] = bouts[i][j].start_frame;
        c["end_frame"] = bouts[i][j].end_frame;
        c["behavior"] = bouts[i][j].behavior;
        b[j] = c;
      }
      a[i] = b;
    }
    r["bouts"] = a;
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
    remove(fname);
    return ret;
  } else {
    ((SVMBehaviorSequence*)s)->init_bout_label(this, NULL);

    if(r["num_bouts"].isArray() && (int)r["bouts"].size() == behaviors->num &&
       (int)r["num_bouts"].size() == behaviors->num) {
      ((SVMBehaviorSequence*)s)->init_bout_label(this, NULL);
      
      for(int i = 0; i < (int)r["num_bouts"].size(); i++) {
	num_bouts[i] = r["num_bouts"][i].asInt();
	this->bouts[i] = (BehaviorBout*)malloc(sizeof(BehaviorBout)*num_bouts[i]);
	for(int j = 0; j < (int)r["bouts"][i].size(); j++) {
	  Json::Value c = r["bouts"][i][j];
	  bouts[i][j].start_frame = c["start_frame"].asInt();
	  bouts[i][j].end_frame = c["end_frame"].asInt();
	  bouts[i][j].behavior = c["behavior"].asInt();
	}
      }
    }
  }
  return true;
}

void SVMBehaviorSequence::DebugFeatures(const char *fname, BehaviorBoutSequence *y, double *ww, int beh) {
   BehaviorBoutFeatures *b = y->features;

  double *tmp_features = (double*)malloc((num_features+3)*sizeof(double));
  double *class_features = (double*)malloc(sizeof(double)*(num_features+behaviors->behaviors[beh].num_values+2));
  double *class_transitions = class_features + num_features;
  double *unary = class_transitions + behaviors->behaviors[beh].num_values;
  double *duration = unary+1;
  int *inds = (int*)malloc(sizeof(int)*(num_features+behaviors->behaviors[beh].num_values+2));
  g_num_features = num_features;
  g_class_features = class_features;
  g_class_transitions = class_transitions;
  double **class_weights, **transition_weights, *unary_weights, *duration_weights;
  init_weight_pointers(beh, ww, class_weights, transition_weights, unary_weights, duration_weights);
  FILE *fout = fopen(fname, "w");

  y->score = 0;
  y->loss = 0;
  for(int i = 0; i < y->num_bouts[beh]; i++) {
    int t_p = y->bouts[beh][i].start_frame, t = y->bouts[beh][i].end_frame, 
      c_prev = y->bouts[beh][i].behavior, 
      c_next = i<y->num_bouts[beh]-1 ? y->bouts[beh][i+1].behavior : -1;
    psi_bout(b, t_p, t, beh, c_prev, tmp_features, true, false);  // compute bout features

    y->bouts[beh][i].bout_score = 0;
    for(int k = 0; k < num_features; k++) {
      class_features[k] = class_weights[c_prev][k]*tmp_features[k];  // update bout score
      y->bouts[beh][i].bout_score += class_features[k];
    }
    for(int j = 0; j < behaviors->behaviors[beh].num_values; j++) 
      class_transitions[j] = 0; 
    class_transitions[c_next] = transition_weights[c_prev][c_next];
		      
    y->bouts[beh][i].transition_score = c_next >=0 ? transition_weights[c_prev][c_next] : 0;   
    *unary = y->bouts[beh][i].unary_score = unary_weights[c_prev];
    *duration = y->bouts[beh][i].duration_score = get_duration_score(t-t_p, c_prev, beh, duration_weights);
    double score = y->bouts[beh][i].bout_score + y->bouts[beh][i].transition_score + y->bouts[beh][i].unary_score + y->bouts[beh][i].duration_score;
    y->score += score;
    y->loss += y->bouts[beh][i].loss_fn + y->bouts[beh][i].loss_fp;
		    
    for(int j = 0; j < num_features+behaviors->behaviors[beh].num_values+2; j++)
      inds[j] = j;
    qsort(inds, num_features+behaviors->behaviors[beh].num_values+2, sizeof(int), cmp_feature_inds);
    fprintf(fout, "<br><br><a name=\"%d\">Behavior=%s score=%f t=(%d-%d)</a>: bout_score=%f transition_score=%f unary_score=%f duration_score=%f loss_fn=%f loss_fp=%f\n", i, behaviors->behaviors[beh].values[c_prev].name, (float)score, t_p, t,
	    (float)y->bouts[beh][i].bout_score, (float)y->bouts[beh][i].transition_score, (float)y->bouts[beh][i].unary_score, (float)y->bouts[beh][i].duration_score, (float)y->bouts[beh][i].loss_fn,  (float)y->bouts[beh][i].loss_fp);
    for(int j = 0; j < num_features+behaviors->behaviors[beh].num_values; j++) {
      if(inds[j] < num_features)
	fprintf(fout, "<br>%lf=%f*%f, %s:%s %s\n", class_features[inds[j]], class_weights[c_prev][inds[j]], tmp_features[inds[j]], behaviors->behaviors[beh].name, behaviors->behaviors[beh].values[c_prev].name, feature_names[inds[j]]);
      else if(inds[j] < num_features+behaviors->behaviors[beh].num_values)
	fprintf(fout, "<br>%lf, %s:%s->%s\n", class_transitions[inds[j]-num_features], behaviors->behaviors[beh].name, behaviors->behaviors[beh].values[c_prev].name, behaviors->behaviors[beh].values[inds[j]-num_features].name);
      else if(inds[j] < num_features+behaviors->behaviors[beh].num_values+1)
	fprintf(fout, "<br>%lf, %s:%s unary\n", *unary, behaviors->behaviors[beh].name, behaviors->behaviors[beh].values[c_prev].name);
      else 
	fprintf(fout, "<br>%lf, %s:%s duration\n", *duration, behaviors->behaviors[beh].name, behaviors->behaviors[beh].values[c_prev].name);
    }
  }
  free(class_weights);
  free(class_features);
  free(tmp_features);
  free(inds);
  fclose(fout);
}

#include "cv.h"
#include "highgui.h"
#define LABEL_BOUTS 0
 void BehaviorBoutSequence::Visualize(BehaviorGroups *groups, int beh, const char *fname, char *html) { 
	BehaviorGroup *group = &groups->behaviors[beh];
	if(!this->num_bouts[beh]) 
		return;

	int h = LABEL_BOUTS ? 50 : 10;
	int T = this->bouts[beh][this->num_bouts[beh]-1].end_frame;
	CvFont font;
	cvInitFont(&font, CV_FONT_VECTOR0, 0.5f, 0.4f, 0, 2);

	IplImage *img = cvCreateImage(cvSize(T, h), IPL_DEPTH_8U, 3);
	cvZero(img);

	// Draw bouts as colored rectangles
	int i;
	for(i = 0; i < this->num_bouts[beh]; i++) {
		int c = this->bouts[beh][i].behavior;
		cvRectangle(img, cvPoint(this->bouts[beh][i].start_frame,0), cvPoint(this->bouts[beh][i].end_frame,h), 
			CV_RGB((group->values[c].color & 0xff0000)>>16,
			(group->values[c].color & 0x00ff00)>>8,
			(group->values[c].color & 0xff)), CV_FILLED);
	}

#if LABEL_BOUTS
	// Draw labels for bouts
	int prev_prev_max_x = -100000, prev_max_x = -100000;
	int last_pos = 2, last_last_pos = 1;
	int y[3], ymin;
	for(i = 0; i < this->num_bouts[beh]; i++) {
		int c = this->bouts[beh][i].behavior;
		int m = (this->bouts[beh][i].start_frame+this->bouts[beh][i].end_frame)/2;

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

	if(fname) {
		char *html_tmp=(char*)malloc(10000000), folder[1000], fname2[1000], fname3[1000];
		if(!strstr(fname, ".png")) sprintf(fname2, "%s.png", fname);
		else strcpy(fname2, fname);
		cvSaveImage(fname2, img);
		cvReleaseImage(&img);

		if(html) {
		  ExtractFolderAndFileName(fname, folder, fname2);    

		  

			sprintf(html_tmp, "<img src=\"%s.png\" usemap=\"#%s\" height=50 />\n<map name=\"%s\">\n", fname2, fname2, fname2);

			char str[10000], alt[1000];
			char *ptr = html_tmp+strlen(html_tmp);
			for(i = 0; i < this->num_bouts[beh]; i++) {
				int c = this->bouts[beh][i].behavior;
				float z = 50.0f/h;
				double score = this->bouts[beh][i].bout_score + this->bouts[beh][i].transition_score + this->bouts[beh][i].unary_score;
#if USE_DURATION_COST > 0
				score += this->bouts[beh][i].duration_score;
#endif
				sprintf(alt, "behavior=%s score=%f bout_score=%f transition_score=%f unary_score=%f duration_score=%f loss_fn=%f loss_fp=%f", group->values[c].name, (float)score,
					(float)this->bouts[beh][i].bout_score, (float)this->bouts[beh][i].transition_score, (float)this->bouts[beh][i].unary_score, (float)this->bouts[beh][i].duration_score, (float)this->bouts[beh][i].loss_fn,  (float)this->bouts[beh][i].loss_fp);
				sprintf(str, "  <area shape=\"rect\" coords=\"%d,%d,%d,%d\" href=\"%s_features.html#%d\" title=\"%s\" />\n", (int)(z*this->bouts[beh][i].start_frame), 0, 
					(int)(z*this->bouts[beh][i].end_frame), (int)(z*h), fname2, i, alt);


				strcpy(ptr, str);
				ptr += strlen(str);
			}
			strcpy(ptr, "</map>");
			strcpy(html, html_tmp);
			free(html_tmp);
		}
	}

}
