#include <stdio.h>
#include <string.h>
#include "svm_behavior_sequence.h"

#define DEBUG 0
//#define ALLOW_SAME_TRANSITIONS 0
#if DEBUG > 0
char *g_currFile; // CSC 20110420: hack to pass current filename for debug purposes
#endif

#define INST_VERSION "0.0.0"


#define compactUnarySize(p_features, p_classes) ((p_features) + (p_classes)    ) * (p_classes)
#define   extraUnarySize(p_features, p_classes) ((p_features) + (p_classes) + 1) * (p_classes)
/*
#ifdef ALLOW_SAME_TRANSITION
#define getPsiSize(p_features, p_classes) compactUnarySize(p_features, p_classes)
#else
#define getPsiSize(p_features, p_classes) extraUnarySize(p_features, p_classes)
#endif
*/
/*
#ifdef ALLOW_SAME_TRANSITIONS  //SJB
#define getPsiSize(p_features, p_classes) extraUnarySize(p_features, p_classes)
#else */
#define getPsiSize(p_features, p_classes) compactUnarySize(p_features, p_classes)
// #endif

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
        eps = .05;
	C = 10.0;
	featureScale = 1;
	method = SPO_DUAL_UPDATE_WITH_CACHE;
	debugLevel = 3;

	strcpy(debugdir, "");
	debug_predictions = debug_weights = debug_features = debug_model = false;

	this->behaviors = behaviors;
	this->behavior = beh;
	this->feature_diff_frames = 10;
#ifdef ALLOW_SAME_TRANSITIONS
	time_approximation = -1;
#else
	time_approximation = 0; // CSC: test EVERY Frame, no log time-saving unless optimality may still be guaranteed!
#endif
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

		// It is assumed here that the label 0 is the unknown label.  The user could define custom
		// class-specific values for these constants .
		// In the future, we could implement this with a per class-pair confusion cost
		false_negative_cost[bh][0] = 1;
		false_positive_cost[bh][0] = 1; // CSC: Move FP/FN costs to Params/BehaviorParams.txt? or to Model/Model.txt?
		for(i = 1; i < num_classes[bh]; i++) {
			false_negative_cost[bh][i] = 1;
			false_positive_cost[bh][i] = 1;
		}

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

void SVMBehaviorSequence::saveBoutFeatures(StructuredDataset *dataset, const char *filename, bool sphered, bool addRandBouts) {
	// EYRUN: Store normalized bout features of training set in a txt file (want to verify that bout features are good enough to separate the classes)
	// Count the total number of bouts in the training set
	int num_bouts = 0;
	int num_rand_factor = 5;
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
		num_bouts_extra += floor(((BehaviorBoutFeatures*)ex[n]->x)->num_frames/30);
	}

	if(addRandBouts)
		num_bouts = num_bouts * (num_rand_factor+1); //num_bouts = num_bouts + num_bouts_extra;

	FILE *featureFile, *featureFileUnsphered;
	featureFile = fopen(filename,"w");
	fprintf(featureFile, "num_features = %d\nnum_bouts = %d", num_features, num_bouts);
	for(int n = 0; n < dataset->num_examples; n++) {
		y = ((BehaviorBoutSequence*)ex[n]->y);
		x = ((BehaviorBoutFeatures*)ex[n]->x);
		for(int beh = 0; beh < behaviors->num; beh++) {
			if(behavior < 0 || beh == behavior) {
				for(int j = 0; j < y->num_bouts[beh]; j++) {
					psi_bout(x, y->bouts[beh][j].start_frame, y->bouts[beh][j].end_frame, beh, y->bouts[beh][j].behavior, tmp_features, sphered, false);
					fprintf(featureFile, "\n%d %d %d %d ", n, y->bouts[beh][j].start_frame, y->bouts[beh][j].end_frame, y->bouts[beh][j].behavior);
					for(int i = 0; i < num_features; i++) {
						fprintf(featureFile, "%f ", tmp_features[i]);
					}
				}
			}
		}
		if(addRandBouts) {
			// Generate random intervals and calculate the bout-features for those
			int num_rand = y->num_bouts[0] * num_rand_factor;
			//num_rand = floor(x->num_frames/30);
			for(int r = 0; r<num_rand; r++) {
				int f_start, f_dur, f_end;
				// generate random start frame between 0 and x->num_frames - 5
				f_start = rand() % (x->num_frames - 5);
				// generate random duration between 4 and 50 (the number of frames we go back)
				f_dur = rand() % 47 + 4;
				f_end = my_min(f_start + f_dur - 1, x->num_frames - 1);

				//f_start = r*30;
				//f_end = (r+1)*30-1;

				// calculate bout score for that interval
				psi_bout(x, f_start, f_end, 0, -1, tmp_features, sphered, false);
				fprintf(featureFile, "\n%d %d %d %d ", n, f_start, f_end, -1);
				for(int i = 0; i < num_features; i++) {
					fprintf(featureFile, "%f ", tmp_features[i]);
				}
			}	
		}	
	}
	fclose(featureFile);
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
	for(n = 0; n < dataset->num_examples; n++)
		for(beh = 0; beh < behaviors->num; beh++)
			if(behavior < 0 || beh == behavior) 
				num_bouts += ((BehaviorBoutSequence*)ex[n]->y)->num_bouts[beh];

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
	histogram_thresholds = (double**)malloc(3*sizeof(double*)*num_base_features + 
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

	feature_names = (char**)malloc(sizeof(char*)*num_features);
	compute_feature_space_size();

	// Now that histogram thresholds are set, we can compute bout-wise features for the training set
	for(n = 0; n < dataset->num_examples; n++) {
		x = (BehaviorBoutFeatures*)ex[n]->x;
		x->ComputeCaches(this);
	}

	// Now that we have bout features, we can compute the mean of each bout-wise feature over the training set
	features_mu = (double*)malloc(num_features*sizeof(double)*2);
	features_gamma = features_mu + num_features;
	for(i = 0; i < num_features; i++) 
		features_mu[i] = features_gamma[i] = 0;
	for(n = 0; n < dataset->num_examples; n++) {
		y = ((BehaviorBoutSequence*)ex[n]->y);
		x = ((BehaviorBoutFeatures*)ex[n]->x);
		for(beh = 0; beh < behaviors->num; beh++) {
			if(behavior < 0 || beh == behavior) {
				for(j = 0; j < y->num_bouts[beh]; j++) {
					psi_bout(x, y->bouts[beh][j].start_frame, y->bouts[beh][j].end_frame, beh, y->bouts[beh][j].behavior, tmp_features, false, false);
					for(i = 0; i < num_features; i++) {
						assert(!isnan(tmp_features[i]));
						features_mu[i] += tmp_features[i];
						assert(!isnan(features_mu[i]));
					}
				}
			}
		}
	}
	for(i = 0; i < num_features; i++) 
		features_mu[i] /= num_bouts;

	// Compute the (inverse of the) standard deviation bout-average feature response
	for(n = 0; n < dataset->num_examples; n++) {
		y = ((BehaviorBoutSequence*)ex[n]->y);
		x = ((BehaviorBoutFeatures*)ex[n]->x);
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
	}
	for(i = 0; i < num_features; i++) {
		if(features_gamma[i] == 0)
			fprintf(stderr, "WARNING: feature %s seems to be the same for every example in the training set\n", feature_names[i]);
		features_gamma[i] = features_gamma[i] ? 1.0/sqrt(features_gamma[i]/num_bouts) : 0;
	}

	saveBoutFeatures(dataset, "train_bout_feat.txt", true, true);
	saveBoutFeatures(dataset, "train_bout_feat_unsphered.txt", false, true);

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
double *SVMBehaviorSequence::psi_bout(BehaviorBoutFeatures *b, int t_start, int t_end, int beh, int c, double *feat, bool normalize, bool fast_update, double extreme_vals[2][NUMFEAT], int start_prev, int end_prev) {
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

	p = &feature_params[0];
	if(p->use_bout_max_feature | p->use_bout_min_feature) { // Eyrun: only update cache if needed (note this assumes that all feature stats are the same as for 0)		
		if(!fast_update)
			b->bout_start = b->bout_end = -1;
		else
			assert(b->bout_end == t_end);

		if(extreme_vals) {
			if (end_prev >= t_end)
				printf("something's wrong...");
			b->UpdateCaches(end_prev, t_end, c);
			for(i=0; i<num_base_features; i++) {
				max_vals[i] = my_max(extreme_vals[0][i],b->bout_max_feature_responses[i]);
				min_vals[i] = my_min(extreme_vals[1][i],b->bout_min_feature_responses[i]);
			}	
		}
		else {
			b->UpdateCaches(t_start, t_end, c);
			for(i=0; i<num_base_features; i++) {
				max_vals[i] = b->bout_max_feature_responses[i];
				min_vals[i] = b->bout_min_feature_responses[i];
			}
		}
	}

	assert(dur > 0);
	double log_dur = log10(dur);
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
			ADD_FEATURE(feat, ind, b->bout_max_feature_responses[i] > max_thresholds[i][l] ? 1 : 0, NUM_BOUT_MAX_THRESHOLDS, i, j);
		for(l = 0; l < p->num_bout_min_thresholds; l++) 
			ADD_FEATURE(feat, ind, b->bout_min_feature_responses[i] < min_thresholds[i][l] ? 1 : 0, NUM_BOUT_MIN_THRESHOLDS, i, j);


		// These are differences in the bout average or sum response for the bout as compared to the global max, min,
		// and average over the entire behavioral sequence
		if(p->use_global_difference_max_sum_features) 
			ADD_FEATURE(feat, ind, sum - b->max_feature_responses[i]*(t_end-t_start), USE_GLOBAL_DIFFERENCE_MAX_SUM_FEATURES, i, 0);
		if(p->use_global_difference_max_ave_features) 
			ADD_FEATURE(feat, ind, ave - b->max_feature_responses[i], USE_GLOBAL_DIFFERENCE_MAX_AVE_FEATURES, i, 0);
		if(p->use_global_difference_min_sum_features) 
			ADD_FEATURE(feat, ind, sum - b->min_feature_responses[i]*(t_end-t_start), USE_GLOBAL_DIFFERENCE_MIN_SUM_FEATURES, i, 0);
		if(p->use_global_difference_min_ave_features) 
			ADD_FEATURE(feat, ind, ave - b->min_feature_responses[i], USE_GLOBAL_DIFFERENCE_MIN_AVE_FEATURES, i, 0);
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
	assert(ind == num_features);

	// Normalize each feature to have mean 0 and standard deviation 1 (as determined by the training set)
	//normalize = false;
	if(normalize) {
		for(i = 0; i < ind; i++)  {
			if(!(abs(feat[i]) <100000000) & isnan(feat[i])){   // EYRUN: inrcreased to 100,000,000 from 10,000,000
				printf("trouble...");
			}
			assert(!isnan(feat[i]));
			assert(abs(feat[i]) < 100000000);
			feat[i] = (feat[i]-features_mu[i])*features_gamma[i];
			//double div = my_max(b->max_feature_responses[i],-b->min_feature_responses[i]);
			//feat[i] = feat[i]/(div < 0.001 ? 1 : div);
			if(!(abs(feat[i]) <100000000) | isnan(feat[i])){
				printf("trouble...");
			}
			assert(!isnan(feat[i]));
			assert(abs(feat[i]) < 100000000);
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
	BehaviorBoutSequence *y = (BehaviorBoutSequence*)yy;
	BehaviorBoutFeatures *b = (BehaviorBoutFeatures*)x;
	int i, j, beh;
	double *tmp_features = (double*)malloc(sizeof(double)*(sizePsi+num_features+NUMFEAT));
	double *all_features = tmp_features+num_features+10;
	double *ptr = all_features, *class_features, *class_transitions;
/*
#ifdef ALLOW_SAME_TRANSITIONS
	double *class_unary;
#endif */

	for(i = 0; i < sizePsi; i++)
		ptr[i] = 0;

	for(beh = 0; beh < behaviors->num; beh++) {
		if(behavior >= 0 && beh != behavior) {
			ptr += getPsiSize(num_features, num_classes[beh]); //num_classes[beh]*(num_classes[beh]+num_features);
			continue;
		}
		class_features = ptr; ptr += num_features*num_classes[beh];        // line 1
		class_transitions = ptr; ptr += num_classes[beh]*num_classes[beh]; // line 2
#ifdef ALLOW_SAME_TRANSITIONS
// 		class_unary = ptr; ptr += num_classes[beh];                        // line 3: these 3 lines have to correspond to getPsiSize
// 
// 		for(i = 0; i < num_classes[beh]; i++)
// 			class_unary[i] = 0;
#endif

		for(i = 0; i < num_classes[beh]*num_features; i++) 
			class_features[i] = 0;

		for(i = 0; i < num_classes[beh]*num_classes[beh]; i++) 
			class_transitions[i] = 0;

		for(i = 0; i < y->num_bouts[beh]; i++) {
/* BEGIN BLOCK */
// 			// Count of the number transitions from class C1 to C2, for all pairs (C1,C2), C!=C2
// 			if(i)
// 				class_transitions[y->bouts[beh][i-1].behavior*num_classes[beh] + y->bouts[beh][i].behavior]++;
// 
// #ifdef ALLOW_SAME_TRANSITIONS
// 			class_unary[y->bouts[beh][i].behavior]++;
// 			/* Note: Same class transition weights were already computed above and are simply NOT OVERWRITTEN by unary costs */
// #ifndef SOLVE_UNARY_COST
// 				; // transition diagonal still counts same class transitions; transition diagonal not used when classifying with this option
// #elif SOLVE_UNARY_COST_CONFLICT == 1
// 				class_transitions[y->bouts[beh][i].behavior*num_classes[beh] + y->bouts[beh][i].behavior] = 0; // transition diagonal is added to unary costs when classifying with this option
// #elif SOLVE_UNARY_COST_CONFLICT == 0
// 				; // transition diagonal still counts same class transitions; unary costs not used when classifying with this option
// 				//class_transitions[y->bouts[beh][i].behavior*num_classes[beh] + y->bouts[beh][i].behavior]++; // transition diagonal is overwritten by unary costs
// #endif
// #else // old compact format using diagonal of transition matrix
// 			// Count of the number of bouts for each class C
// 			class_transitions[y->bouts[beh][i].behavior*num_classes[beh] + y->bouts[beh][i].behavior]++; // transition diagonal is overwritten by unary costs
// #endif

			if (i) {
				// make sure that we don't add transitions for self_transitions
				if (y->bouts[beh][i-1].behavior != y->bouts[beh][i].behavior) {
					class_transitions[y->bouts[beh][i-1].behavior*num_classes[beh] + y->bouts[beh][i].behavior]++; // increase the transition cost
					class_transitions[y->bouts[beh][i].behavior*num_classes[beh] + y->bouts[beh][i].behavior]++;   // increase the unary cost
				}
			}
			else
				class_transitions[y->bouts[beh][i].behavior*num_classes[beh] + y->bouts[beh][i].behavior]++; 

/* END BLOCK */



/* the BLOCK above corresponds to the following lines (copied from above)

#ifdef ALLOW_SAME_TRANSITIONS
		double *unary_costs = (double*)malloc(num_classes[beh]*sizeof(double*));
			for(i = 0; i < num_classes[beh]; i++, ptr++) {
				if(sm->compactUnaryCosts) {
					unary_costs[i] = transition_weights[i][i]; // same transition costs are unary costs
					transition_weights[i][i] = 0; // allowing same class transitions, since no transition cost values were given same transition costs are initialized with 0 (= no cost for same-class transitions besides bout-level-feature costs)
				}
				else // if(sm.extraUnaryCosts)
					unary_costs[i] = *ptr; // unary costs are explicitely given, same transition costs remain
			}
#else
		if(sm->extraUnaryCosts) // we don't allow same transitions, but have read same-transition-weights AND unary weights
			for(i = 0; i < num_classes[beh]; i++, ptr++)
#ifndef SOLVE_UNARY_COST_CONFLICT
				transition_weights[i][i] = *ptr; // default: replace same transition costs by unary costs
#elif SOLVE_UNARY_COST_CONFLICT == 1 // 0... skip unary costs  1... add unary costs to transition costs
				transition_weights[i][i] += *ptr; // optional: add unary costs to same transition costs
#elif SOLVE_UNARY_COST_CONFLICT == 0
				; // keep same transition costs and just skip unary costs
#else
*/


			// The main appearance features used to compute the score of a bout for a given class.
			// Since the total score is the sum of the scores over bouts, we can simply add together
			// the features of the bouts with the same class label
			psi_bout(b, y->bouts[beh][i].start_frame, y->bouts[beh][i].end_frame, beh, y->bouts[beh][i].behavior, tmp_features, true, false);
			for(j = 0; j < num_features; j++) {
				class_features[y->bouts[beh][i].behavior*num_features+j] += tmp_features[j];
			}
		}

	}
	
	SparseVector retval = SparseVector(all_features, sizePsi);
	free(tmp_features);

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

	int num, i, j, beh;       
	char **train_list = load_examples(fname, &num);
	bool computeClassTransitions = class_training_transitions ? false : true;
	StructuredDataset *dataset = new StructuredDataset();

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
					if(i)
						class_training_transitions[beh][y->bouts[beh][i].behavior][y->bouts[beh][i-1].behavior]++;
				}
			}	 
		} 
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

  if(bouts) {
    strcpy(bouts->fname, label_fname);
    bool b = bouts->load(label_fname);
    assert(b);
  }
  strcpy(feature_cache->fname, features_fname);
  bool b = feature_cache->load(features_fname, this, (BehaviorBoutSequence*)ex->y);
  assert(b);
  if(computeFeatures) {
    feature_cache->ComputeCaches(this);
    feature_cache->fvec = label_fname ? Psi(ex->x, ex->y).ptr() : NULL;
  }
  return ex;
}





void SVMBehaviorSequence::OnFinishedIteration(StructuredData *x, StructuredLabel *y, StructuredLabel *ybar) {
		char fname[1000];
		int iter_num = this->t;
		if(strlen(debugdir) && debug_weights) {
			sprintf(fname, "%s/weights_%d.txt", debugdir, iter_num);
			SparseVector *w = GetCurrentWeights();
			double *ww = w->non_sparse<double>(sizePsi);
			print_weights(fname, ww);
			delete [] ww;
			delete w;
		}

		if(strlen(debugdir) && debug_model) {
			sprintf(fname, "%s/model_%d.txt", debugdir, iter_num);
			StructuredSVM::Save(fname, false);
		}

		if(strlen(debugdir) && debug_predictions) {
			sprintf(fname, "%s/index.html", debugdir);
			FILE *fout = fopen(fname, "a");
			assert(fout);
			fprintf(fout, "<br><br><h2>Iteration %d</h2><a href=\"weights_%d.txt\">weights</a>|<a href=\"iter%d.html\">predictions</a>\n", iter_num, iter_num, iter_num);
			fclose(fout);
		}

	// Save a visualization of all predicted bouts
	if(strlen(debugdir) && debug_predictions) {
		char *html=(char*)malloc(10000000), *html_gt=(char*)malloc(10000000), folder[1000], file[1000], fname[1000];
		ExtractFolderAndFileName(((BehaviorBoutFeatures*)x)->GetFileName(), folder, file);
		StripFileExtension(file);
		int beh = this->behavior >= 0 ? this->behavior : 0;
		if(iter_num >= 0) sprintf(fname, "%s/%s_%d", debugdir, file, iter_num);
		else sprintf(fname, "%s/%s", debugdir, file);
		((BehaviorBoutSequence*)ybar)->Visualize(behaviors, beh, fname, html);
		if(y) {
			sprintf(fname, "%s/%s_gt_%d", debugdir, file, iter_num);
			((BehaviorBoutSequence*)y)->Visualize(behaviors, beh, fname, html_gt);
		} 
		if(iter_num >= 0)
			sprintf(fname, "%s/iter%d.html", debugdir, iter_num);
		else
			sprintf(fname, "%s/index.html", debugdir);
		FILE *fout = fopen(fname, "a");
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
}

#if DEBUG > 0
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

#endif



#define getTransitionScore(dest, unary, addTransition, transition) \
	dest = unary; \
	if(addTransition) \
		dest += transition;
    // Generalization of
	// y->bouts[beh][i].transition_score = transition_weights[y->bouts[beh][i].behavior][y->bouts[beh][i].behavior]; // !!! add unary costs here
	// if(i < y->num_bouts[beh]-1) 
    //		y->bouts[beh][i].transition_score += transition_weights[y->bouts[beh][i].behavior][y->bouts[beh][i+1].behavior];


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
				      StructuredLabel *y_partial, StructuredLabel *y_gt) {
	BehaviorBoutSequence *ybar = (BehaviorBoutSequence*)(y_bar);
	BehaviorBoutSequence *y = (BehaviorBoutSequence*)(y_gt);
	BehaviorBoutFeatures *b = (BehaviorBoutFeatures*)x;
	int T = b->num_frames, beh, num_durations = 0, last_gt, last_partial;
	int tt, t, c, t_p, c_p, i, j, cc, restrict_c, restrict_c_p, next_duration;
	int first_frames[num_classes[0]];
	int idx_max=8, idx_min=9, num_stats=18;
	double extreme_vals[num_classes[0]][2][NUMFEAT];
	double *tmp_features = (double*)malloc((num_features+1)*sizeof(double));
	double *tmp_features_ = (double*)malloc((num_features+1)*sizeof(double));
	//double cached_bout_scores[T][num_classes[0]];
	double **class_weights;
	double **transition_weights;
	double dur, inter;
	int *durations = (int*)malloc(sizeof(int)*(T+1)*4);
	int *gt_bout = durations + (T+1)*2;
	int *partial_label_bout = gt_bout + (T+1);
	double *fn;
	double *dur_gt;
	double *ptr, *ptr3;
	double **table;
	BehaviorBout **states;
	BehaviorBout *ptr2;
	double f;
	double bout_score, transition_score, loss_score, loss_fn, loss_fp, max_fn_cost, fp;
	double bout_scores[1000];
	bool bad_label = false;
	int *old_class_transition_counts = NULL, *old_class_training_counts;
	bool is_first, go_back;
	double *ww = w->non_sparse<double>(sizePsi);
	Json::FastWriter writer;

	sprintf(ybar->fname, "%s.pred", b->fname);
/*
	double total_bouts = 0;
	for (int i=0; i<num_classes[0]; i++)
		total_bouts += class_training_count[0][i];
	for (int i=0; i<num_classes[0]; i++) {
		false_negative_cost[0][i] = total_bouts/class_training_count[0][i];
		false_positive_cost[0][i] = total_bouts/class_training_count[0][i];
	}
*/
	// For speed, we may optionally choose to use a geometrically increasing series of bout durations instead
	// of searching through the entire space exhaustively
	//time_approximation = 0;
	if(T < 50) 
		time_approximation = 0;
	if(time_approximation > 0) {
		dur = min_bout_duration;
		while(dur <= T*(1+time_approximation)) {
			if(!num_durations || (int)(dur) != durations[num_durations-1])
				durations[num_durations++] = (int)(dur);
			dur *= (1+time_approximation);
		}
	} else if(time_approximation < 0) {
		for(i = 1; i <= 50; i++) // go only 50 frames back
			durations[num_durations++] = i;
	}else {
		for(i = 1; i <= T; i++)
			durations[num_durations++] = i;
	}


	// Initialize ybar
	ybar->bouts = (BehaviorBout**)realloc(ybar->bouts, (sizeof(BehaviorBout*)+sizeof(int*)+2*sizeof(double))*behaviors->num);
	ybar->num_bouts = (int*)(ybar->bouts+behaviors->num);
	ybar->scores = (double*)(ybar->num_bouts+behaviors->num);
	ybar->losses = (double*)(ybar->scores+behaviors->num);
	ybar->slack = ybar->score = 0;
	ybar->loss = 0;
	ybar->behaviors = behaviors;
	if(y) {
		ybar->features = y->features;
		ybar->behaviors = y->behaviors;
#if DEBUG > 0
		y->score = 0;
#endif
	}

	
	ptr = ww;
	for(beh = 0; beh < behaviors->num; beh++) {
		// Currently optimizes each behavior group independently
		if(behavior >= 0 && beh != behavior) {
			//ptr += num_classes[beh]*(num_classes[beh]+num_features);
			ptr += getPsiSize(this->num_features, this->num_classes[beh]);
			continue;
		}

#if DEBUG > 1
#define MAX_FILENAME 1000
#define DEBUG__FILESEPARATION 3 // 0.. all in one  1.. per destClass  2.. per currClass  2.. per 500 frames
//#define DEBUG__PRINTNUM
#define DEBUG__STDOUT

								char filename[MAX_FILENAME];
//								char g_currFile[MAX_FILENAME]={0};
								char *currFileWithoutPath = getFilenameWithoutPath(g_currFile);
								int *writeHeader;
								FILE **outA;
								int numA = 0;

#if DEBUG__FILESEPARATION == 0
								numA = 1;
#elif DEBUG__FILESEPARATION == 1
								numA = num_classes[beh];
#elif DEBUG__FILESEPARATION == 2
								numA = num_classes[beh];
#elif DEBUG__FILESEPARATION == 3
#define DEBUG__FILESEPARATION_WINDOW 2
								numA = T / DEBUG__FILESEPARATION_WINDOW + 1;
#endif

								int maxStdio = _setmaxstdio(numA+20);
								assert(maxStdio == numA+20);

								outA = (FILE **)malloc(numA*sizeof(FILE*));
								writeHeader = (int*)malloc(numA*sizeof(int));

								for(cc=0; cc < numA; cc++) {
//									sprintf(filename, "%s\\debug_%s_%d_%d__%d.csv", /*debugdir*/"C:\\tmp\\every_line", currFileWithoutPath, DEBUG__FILESEPARATION, beh, cc);
									sprintf(filename, "%s\\debug_%s_%d_%d__%d.csv", debugdir, currFileWithoutPath, DEBUG__FILESEPARATION, beh, cc);
#ifdef DEBUG___STDOUT
									outA[cc] = stdout;
#else
									outA[cc] = fopen(filename, "w"); // comma-separated list
#endif
									assert(outA[cc]);
									writeHeader[cc] = 1;
								}
								unsigned long ul = 0;
								// print header of comma-separated list later, during first iteration
#endif


		class_weights = (double**)malloc(2*num_classes[beh]*sizeof(double*));
		transition_weights = class_weights+num_classes[beh];
		fn = y ? (double*)malloc(sizeof(double)*y->num_bouts[beh]*2) : NULL;
		dur_gt = fn ? fn + y->num_bouts[beh] : NULL;
		table = (double**)malloc((T+1)*(sizeof(double*)+num_classes[beh]*sizeof(double)));
		states = (BehaviorBout**)malloc((T+1)*(sizeof(BehaviorBout*)+num_classes[beh]*sizeof(BehaviorBout)));
		loss_score = loss_fn = loss_fp = 0;


		if(b->partial_label) {
			// When given a user supplied label, it may be the case that some class labels or label 
			// sequence in the partial labelling never appeared in the training set.  We add a couple of
			// checks here to protect against this case (otherwise our algorithm would find no solution)
			old_class_transition_counts = (int*)malloc(num_classes[beh]*sizeof(int)*2);
			old_class_training_counts = old_class_transition_counts + num_classes[beh];
			for(c = 0; c < num_classes[beh]; c++) {
				old_class_transition_counts[c] = class_training_transitions_count[beh][c];
				old_class_training_counts[c] = class_training_count[beh][c];
			}

			for(i = 0; i < b->partial_label->num_bouts[beh]; i++) {
				c = b->partial_label->bouts[beh][i].behavior;
				c_p = i > 0 ? b->partial_label->bouts[beh][i-1].behavior : -1;

				if(!class_training_count[beh][c])
					class_training_count[beh][c]++;
				if(c >= 0 && c_p >= 0) {
					for(j = 0; j < class_training_transitions_count[beh][c]; j++)
						if(class_training_transitions[beh][c][j] == c_p)
							break;
					if(j == class_training_transitions_count[beh][c])
						class_training_transitions[beh][c][class_training_transitions_count[beh][c]++] = c_p;
				}
			}
		}

		// Setup pointers to class_weights and transition weights.  It is assumed that
		// the first class_weightsXnum_features model parameters correspond to the
		// class-feature weights, and the last num_classesXnum_classes weights correspond
		// to the class transition weights
		for(i = 0; i < num_classes[beh]; i++, ptr += num_features) 
			class_weights[i] = ptr;
		for(i = 0; i < num_classes[beh]; i++, ptr += num_classes[beh]) 
			transition_weights[i] = ptr; // treat unary costs differently here?
/* BLOCK BEGIN // EYRUN: comment this out since we're not using predefined unary costs, even with transitions allowed
#ifdef ALLOW_SAME_TRANSITIONS
		double *unary_costs = (double*)malloc(num_classes[beh]*sizeof(double));
			for(i = 0; i < num_classes[beh]; i++, ptr++) {
//				if(sm->compactUnaryCosts) {
//					unary_costs[i] = transition_weights[i][i]; // same transition costs are unary costs
//					transition_weights[i][i] = 0; // allowing same class transitions, since no transition cost values were given same transition costs are initialized with 0 (= no cost for same-class transitions besides bout-level-feature costs)
//				}
//				else // if(sm.extraUnaryCosts)
					unary_costs[i] = *ptr; // unary costs are explicitely given, same transition costs remain
			}
#else
		if(sm->extraUnaryCosts) // we don't allow same transitions, but have read same-transition-weights AND unary weights
			for(i = 0; i < num_classes[beh]; i++, ptr++)
#ifndef SOLVE_UNARY_COST_CONFLICT
				transition_weights[i][i] = *ptr; // default: replace same transition costs by unary costs
#elif SOLVE_UNARY_COST_CONFLICT == 1 // 0... skip unary costs  1... add unary costs to transition costs
				transition_weights[i][i] += *ptr; // optional: add unary costs to same transition costs
#elif SOLVE_UNARY_COST_CONFLICT == 0
				; // keep same transition costs and just skip unary costs
#else
#error invalid values of SOLVE_UNARY_COST_CONFLICT, expected values are 0 or 1.
#endif
		// else: nothing to do (don't allow same transitions, no extra unary costs read)
#endif
BLOCK END */

		// Setup dynamic programming cache tables.  table[t][c] stores the maximal score for any sub-solution
		// to frames 1...t in which a bout of class c begins at time t
		for(i = 0,  ptr3 = (double*)(table+T+1), ptr2 = (BehaviorBout*)(states+T+1); i <= T; 
			i++, ptr3 += num_classes[beh], ptr2 += num_classes[beh]) {
				table[i] = ptr3;
				states[i] = ptr2;
		}

		// If evaluating loss with respect to a ground truth label y, compute the maximal false negative cost,
		// if every bout in y was missed entirely.  This loss will be subtracted off with each iteration of 
		// dynamic programming
		max_fn_cost = 0;
		if(y) {
			for(i = 0; i < y->num_bouts[beh]; i++) {
				dur_gt[i] = (b->frame_times[my_min(y->bouts[beh][i].end_frame,T-1)] -
					b->frame_times[y->bouts[beh][i].start_frame]);
				fn[i] = match_false_negative_cost(dur_gt[i], beh, y->bouts[beh][i].behavior);
				max_fn_cost += fn[i];
			}
		}

		// Base case: initialize scores to max_fn_cost
		for(c = 0; c < num_classes[beh]; c++) {
			table[0][c] = 0;//max_fn_cost;
			states[0][c].start_frame = states[0][c].end_frame = states[0][c].behavior = -1;
			for(i=0; i<NUMFEAT; i++){
				states[0][c].extreme_vals[0][i] = -INFINITY;
				states[0][c].extreme_vals[1][i] = INFINITY;
			}
		}

		gt_bout[0] = 0;
		partial_label_bout[0] = 0;
		for(t = 1; t <= T; t++) { // looping through all possible bout ends // T .. length of annotation window

			//if(t%1000 == 0)
			//	printf("t = %d\n",t);

			//// Initialize cache for all possible bouts ending at time t
			//for (t_p=0; t_p<t; t_p++) { // possible start times
			//	for (c=0; c<num_classes[0]; c++) { // possible bout types
			//		cached_bout_scores[t_p][c] = -INFINITY;
			//	}
			//}

			for(c = 0; c < num_classes[beh]; c++) { // initialize for all classes
				table[t][c] = -INFINITY;
			}

			// When given a groundtruth label y, stores the index of the bout corresponding to this timestep
			// during training invocation, not classifying invocation
			if(y) gt_bout[t] = y->bouts[beh][gt_bout[t-1]].end_frame < t ? gt_bout[t-1]+1 : gt_bout[t-1];

			// When given a manually supplied partial labelling, store the index of the bout corresponding to this timestep
			if(b->partial_label) 
				partial_label_bout[t] = partial_label_bout[t-1] < b->partial_label->num_bouts[beh] && 
				b->partial_label->bouts[beh][partial_label_bout[t-1]].end_frame <= t ? 
				partial_label_bout[t-1]+1 : partial_label_bout[t-1];

			// Suppose a new bout of class c begins at time t.  Compute the optimal score for any labeling
			// through time 0...t that begins a bout of class c in time t.  We can prune the search space, 
			// because given the preceding completed bout (a bout beginning in some time step t_p and of some
			// class c_p) computation of the score function w.r.t. all frames at time t and onwards does not 
			// depend on any frames before t.  We can therefore enumerate all possible completed bouts and 
			// take the one with the highest score
			if(y) last_gt = gt_bout[t];
			if(b->partial_label) last_partial = partial_label_bout[t];

			next_duration = 1;
			tt = 0; t_p = t;
			is_first = true;


			while(next_duration >= 0 && tt < num_durations) { // loop through all possible bout starts (resp. durations)


				// Try all possible durations.  We may choose to add extra durations to make sure we explore the solutions
				// contained in y or b->partial_label
				if(y && (last_gt >= y->num_bouts[beh] || t_p == y->bouts[beh][last_gt].start_frame))
					last_gt--;
				if(b->partial_label && b->partial_label->num_bouts[beh] && 
					(last_partial >= b->partial_label->num_bouts[beh] || t_p == b->partial_label->bouts[beh][last_partial].start_frame))
					last_partial--;

				next_duration = 1;
				t_p = t-durations[tt];
				// BEGIN squeeze in partial duration that is not in traversed durations (some durations are skipped, durations increase geometrically)
// 				if(y && (t_p < 0 ? -1 : gt_bout[t_p]) != last_gt && t_p != y->bouts[beh][last_gt].start_frame) {
// 					t_p = y->bouts[beh][last_gt].start_frame;
// 					next_duration = 0;
// 				}
// 				if(b->partial_label && b->partial_label->num_bouts[beh] && 
// 					(t_p < 0 ? -1 : partial_label_bout[t_p]) != last_partial && t_p != b->partial_label->bouts[beh][last_partial].start_frame) {
// 						t_p = b->partial_label->bouts[beh][last_partial].start_frame;
// 						next_duration = 0;
// 				}
				// END squeeze
				tt += next_duration;
				if(t_p <= 0) {
					t_p = 0;
					next_duration = -1;
				}
				if(t_p >=t | t_p< 0)
					printf("Upcoming Error...\n");
				assert(t_p < t && t_p >= 0);

				restrict_c_p = restrict_c = -1;
				if(b->partial_label) {
					// When a user-supplied partial label is specified, we need to restrict the set of
					// possible classes based on the applicable labelled bouts
					if(partial_label_bout[t] < b->partial_label->num_bouts[beh])
						restrict_c = b->partial_label->bouts[beh][partial_label_bout[t]].behavior;

					bad_label = false;
					for(i = partial_label_bout[t_p]; i < b->partial_label->num_bouts[beh] && 
						b->partial_label->bouts[beh][i].start_frame < t; i++)  {
							j = b->partial_label->bouts[beh][i].behavior;
							if(j && (restrict_c_p < 0 || j != restrict_c_p)) {
								if(restrict_c_p > 0)
									bad_label = true;
								restrict_c_p = j;
							}
					}

					// When given a partial labelling, there is some potential computational savings by cutting 
					// duration search spaces extending past multiple labelled bouts
					if(bad_label) // multiple user-given partial labels withing current start-to-end bout, will conflict for sure for currently selected bout duration
						break;
				}

				// For speed, CURRENTLY ASSUMING ALL CLASSES HAVE THE SAME BASIC FEATURE SPACE.  To get rid of this
				// assumption, move this below the for(c_p...) loop and call psi_bout(b, t_p, t, c_p, tmp_features).  
				// These loops were intentially ordered in a semi-funny way to make sure the computations
				// psi_bout(t_p,t) and w_c*psi() are computed as few times as possible
#ifndef MERGE_SAME_BOUTS
				psi_bout(b, t_p, t, beh, -1, tmp_features, true, !is_first); // compute bout-level features // !!! use startframe of c in case c==c_p
#else // then merge bouts
//				psi_bout(b, last__t_p[c], t, beh, -1, tmp_features, true, !is_first); // tweak t_p resp. t instead
#endif
				
//for(int q = 0; q < num_features; q++) { assert(abs(tmp_features[q]) < 1000000); }

				is_first = false;
				for(c = 0; c < num_classes[beh]; c++) {
					if((t == T && c != 0) || (t < T && (!class_training_count[beh][c])))
						continue;  // class c doesn't appear in the training set					

					go_back = false;
					first_frames[c] = t_p;
#ifdef ALLOW_SAME_TRANSITIONS
					if(states[t_p][c].behavior == c) {
						first_frames[c] = states[t_p][c].start_frame;
						//if(cached_bout_scores[first_frames[c]][c] > -INFINITY) {
						//	bout_scores[c] = cached_bout_scores[first_frames[c]][c];
						//	continue;
						//}
						psi_bout(b, first_frames[c], t, beh, -1, tmp_features_, true, !is_first, states[t_p][c].extreme_vals, states[t_p][c].start_frame, states[t_p][c].end_frame);
						for(i=0; i<num_base_features; i++) {
							extreme_vals[c][0][i] = my_max(states[t_p][c].extreme_vals[0][i],b->bout_max_feature_responses[i]);
							extreme_vals[c][1][i] = my_min(states[t_p][c].extreme_vals[1][i],b->bout_min_feature_responses[i]);
						}
						go_back = true;
					}
					else {
						for(i=0; i<num_base_features; i++) {
							extreme_vals[c][0][i] = b->bout_max_feature_responses[i];
							extreme_vals[c][1][i] = b->bout_min_feature_responses[i];
						}
					}
#endif
					if(class_training_count[beh][c] && (restrict_c_p<=0 || c == restrict_c_p)) {
						if(restrict_behavior_features[beh] && restrict_behavior_features[beh][c]) {
							printf("restrict_behavior_features currently not used\n");
							assert(0);
							// untested
							bout_scores[c] = 0;
							for(int k = 0; k < num_features; k++) 
							  if(!restrict_behavior_features[beh][c][k]) {
							    if(go_back)
							      bout_scores[c] += class_weights[c][k]*tmp_features_[k];
						            else
						              bout_scores[c] += class_weights[c][k]*tmp_features[k];
							  }
						}
						else {
							bout_scores[c] = 0;
							for(int k = 0; k < num_features; k++) {
						          if(go_back)
							    bout_scores[c] += class_weights[c][k]*tmp_features_[k];
						          else
							    bout_scores[c] += class_weights[c][k]*tmp_features[k];
							}
							// function f_b(x,s,e) (8), which is summed to compute f(x,y) (7)
						}
					}
					//cached_bout_scores[first_frames[c]][c] = bout_scores[c];
				}

				for(c = 0; c < num_classes[beh]; c++) { // iterate through all classes (starting at t-next_duration, ending at t)

					if((t == T && c != 0) || (t < T && (!class_training_count[beh][c])))
						continue;  // class c doesn't appear in the training set

					if(restrict_c>0 && restrict_c != c)
						continue;  // class c doesn't agree with a user-supplied partial labeling

					// Consider a bout of class c_p beginning at frame t_p, ending at frame t, and transitioning
					// into a new bout of class c

#ifndef ALLOW_SAME_TRANSITIONS
					for(cc = 0; cc < (t == T ? num_classes[beh] : class_training_transitions_count[beh][c]); cc++) { // iterate through all classes in bouts FOLLOWING the current bout of class c
#else
					for(cc = -1; cc < (t == T ? num_classes[beh] : class_training_transitions_count[beh][c]); cc++) { // iterate through all classes in bouts FOLLOWING the current bout of class c
						if(cc == -1)
							c_p = c;
						else {
#endif
							// Consider only class transition pairs that occur in the training set
							c_p = (t == T ? cc : class_training_transitions[beh][c][cc]); 
						// end if (cc == -1) in #ifdef ALLOW_SAME_TRANSITION
#ifdef ALLOW_SAME_TRANSITIONS
						}
						if(c_p == c && cc != -1)
							continue; // already computed that at iteration -1, don't compute it again
#endif
						if(!class_training_count[beh][c_p])
							continue;
						t_p = first_frames[c_p];
						if(restrict_c_p>0 && restrict_c_p != c_p)
							continue;  // class c_p doesn't agree with a user-supplied partial labeling

						if(t < T || c == 0) {  // t==T,c==0 is a special case to extract the final optimal solution  
							// if t==T, then there is no "following" bout...

							// Compute classification score as a function of the raw input features and the bout's class,
							// start, end end
							bout_score = bout_scores[c_p]; 


// This block is replaced by getTransitionScore BEGIN
							// Compute the score due to transitioning from class c_p to class c
//							transition_score = t != T ? transition_weights[c_p][c] : 0;

							// transition_weights[c_p][c_p] is a special case which is used to stored a unary cost
							// for adding a bout of class c_p.  In general, this should add a penalty that encourages
							// segmentations with fewer total bouts
//#ifdef ALLOW_SAME_TRANSITIONS
//							transition_score += unary_costs[c_p]; // adding "unary" cost (class specific appearance cost)
//#else
//							transition_score += transition_weights[c_p][c_p]; // adding "unary" cost (class specific appearance cost)
//#endif
// This block is replaced by getTransitionScore END

// #ifdef ALLOW_SAME_TRANSITIONS
// 							//getTransitionScore(transition_score, transition_weights[c_p][c_p], unary_costs[c_p], 1, t != T ? transition_weights[c_p][c] : 0)
//                             				getTransitionScore(transition_score, unary_costs[c_p], 1, t != T ? transition_weights[c_p][c] : 0);  // SJB
// #else
							getTransitionScore(transition_score, transition_weights[c_p][c_p], 1, t != T ? transition_weights[c_p][c] : 0);
//#endif


							if(y) { // if invoked in LEARNING mode
								// If a ground truth label y is given, update the componenent of loss(y,ybar) that is 
								// attributable to the completed bout ybar_i=(c_p,t_p,t)

								// computes (5) // !!! use startframe of c in case c==c_p

								dur = b->frame_times[my_min(t,T-1)]-b->frame_times[t_p];
								loss_fp = fp = match_false_positive_cost(dur, beh, c_p); // l^b_fp in (5) // start at maximum FP cost and...
								loss_fn = 0;
								for(j = gt_bout[t_p]; j <= gt_bout[t] && j < y->num_bouts[beh]; j++) {
									// Loop through all ground truth bouts intersecting with ybar_i
									inter = (b->frame_times[my_min(my_min(y->bouts[beh][j].end_frame,t),T-1)] - 
										b->frame_times[my_max(y->bouts[beh][j].start_frame,t_p)]);
									if(inter <= 0) continue;
									if(c_p == y->bouts[beh][j].behavior) {
										//loss_fn -= fn[j]*inter/dur_gt[j]; // loss_fm is a negative number (initialized 0, then subtracting), 
										// the absolute false negative cost is not used but rather the relative FN cost, 
										// where the "common shared" constant is left out, 
										// but since we are only interested in the maximum, that doesn't matter
										loss_fp -= fp*inter/dur; // ... and subtract agreeing frames from that maximum
									}
									else{
										loss_fn += fn[j]*inter/dur_gt[j];
									}
								}
								loss_score = loss_fn + loss_fp;
							} // end if (y)

							// Check if the completed bout has a higher score than all previously examined solutions that 
							// begin a bout of class c at time t
							f = table[t_p][c_p] + bout_score + transition_score + loss_score;
							assert(!isnan(f));

							if(f > table[t][c]) {  
								table[t][c] = f;						
								states[t][c].start_frame =  t_p;
								states[t][c].end_frame =  t;
								states[t][c].behavior = c_p;	
								// states stores start/end/label for every end-frame and label, rest is debug information
								states[t][c].bout_score = bout_score;		
								states[t][c].transition_score = transition_score;
								states[t][c].loss_fn = loss_fn;
								states[t][c].loss_fp = loss_fp;
								for(i=0; i<num_base_features; i++) {
									states[t][c].extreme_vals[0][i] = extreme_vals[c_p][0][i];
									states[t][c].extreme_vals[1][i] = extreme_vals[c_p][1][i];
								}
							} // end if(f > table[t][c]) {  
#if DEBUG > 1
							if(!y) { // do this in classification mode only
								assert(loss_score==0);
								assert(bout_score==bout_scores[c_p]);
								int currI;

#if DEBUG__FILESEPARATION == 0
								currI = 0;
#elif DEBUG__FILESEPARATION == 1
								currI = c_p;
#elif DEBUG__FILESEPARATION == 2
								currI = c;
#elif DEBUG__FILESEPARATION == 3
								currI = t / DEBUG__FILESEPARATION_WINDOW;
#endif

								FILE *out2 = outA[currI];
								if(writeHeader[currI]) {
/*									fprintf(out2, ",,,,,,,,,,,,, ");
									for(i=0; i < num_features; i++)
										fprintf(out2, "%s, %s", base_feature_names[g_feature_map2[i]], base_feature_names[g_feature_map2[i]]);
									fprintf(out2, "\n");*/
#ifdef DEBUG__PRINTNUM
									fprintf(out2, "#, t, t_p, c, c_p, table[t][c], f, table[t_p][c_p], transition_score, bout_score, loss_score, loss_fp, loss_fn, transition_weights[c_p][c], transition_weights[c_p][c_p], ");
#else
									fprintf(out2, "t, t_p, c, c_p, table[t][c], f, table[t_p][c_p], transition_score, bout_score, loss_score, loss_fp, loss_fn, transition_weights[c_p][c], transition_weights[c_p][c_p], ");
#endif
									for(i=0; i < num_features; i++)
										fprintf(out2, "W %s, V %s, ", feature_names[i], feature_names[i]);//bout_feature_names[g_feature_map[i]], bout_feature_names[g_feature_map[i]]);
									fprintf(out2, "backtrack\n");

//									fprintf(out2, ",,,,,,,,,,,,, ");
//									for(i=0; i < num_features; i++)
//										fprintf(out2, "weight %d, value %d", g_feature_map3[i], g_feature_map3[i]);

									writeHeader[currI] = 0;
								}

#ifdef DEBUG__PRINTNUM
								fprintf(out2, "%d, %d, %d, %d, %d, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, ",
										ul++, t, t_p, c, c_p, table[t][c], f, table[t_p][c_p], transition_score, bout_score, loss_score, loss_fp, loss_fn, transition_weights[c_p][c], transition_weights[c_p][c_p]);
#else
								fprintf(out2, "%d, %d, %d, %d, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, ",
										t, t_p, c, c_p, table[t][c], f, table[t_p][c_p], transition_score, bout_score, loss_score, loss_fp, loss_fn, transition_weights[c_p][c], transition_weights[c_p][c_p]);
#endif
								for(i=0; i < num_features; i++)
									fprintf(out2, "%f, %f, ", class_weights[c_p][i], tmp_features[i]);
							    fprintf(out2, "0\n");
								fflush(out2);


#if DEBUG > 3
								sprintf(filename, "%s/debug_f_%d_%d_%d_%d.csv", debugdir, t, t_p, c, c_p);
								FILE *out = fopen(filename, "w");
								assert(out);

								// print comma separated list of summands
								fprintf(out, "'name', 'weight', 'value'\n");
								fprintf(out, "'f [%d][%d]', 1, %f\n", t, c, "equals", f);
								fprintf(out, "'table[%d}[%d]', 1, %f\n", t_p, c_p, table[t_p][c_p]);
								fprintf(out, "'transition_weights[%d][%d]', %f, 1\n", c_p, c, t != T ? transition_weights[c_p][c] : 0);
								fprintf(out, "'transition_weights[%d][%d]', %f, 1\n", c_p, c_p, transition_weights[c_p][c_p]);
								for(i=0; i < num_features; i++)
									fprintf(out, "'%s', %f, %f\n", bout_feature_names[g_feature_map[i]], class_weights[c_p][i], tmp_features[i]);
								fclose(out);
#endif
#if DEBUG > 2
								printf("f = %f { table[%d][%d] {=%f} + transition_score {=%f = %f + %f} + bout_score {=%f =", t, c, f, t_p, c_p, table[t_p][c_p], transition_score, t != T ? transition_weights[c_p][c] : 0, transition_weights[c_p][c_p], bout_scores[c_p]);
								for(i=0; i < num_features; i++)
									printf("%s_score {=%f = weight %f * value %f} + ", bout_feature_names[g_feature_map[i]], class_weights[c_p][i]*tmp_features[i], class_weights[c_p][i], tmp_features[i]);
								printf("))\n");
#endif
							}
#endif							
						} // end if(t < T || c == 0) {
					} // end for(cc = 0; cc < (t == T ? num_classes[beh] : class_training_transitions_count[beh][c]); cc++) {  // iterate through all classes in bouts FOLLOWING the current bout of class c
				} //for(c = 0; c < num_classes[beh]; c++) { // iterate through all classes (starting at t-next_duration, ending at t)
			} // while(next_duration >= 0 && tt < num_durations) { // loop through all possible bout starts (resp. durations)
		} // for(t = 1; t <= T; t++) { // looping through all possible bout ends // T .. length of annotation window

		// First backtrack to count the number of bouts in the optimal solution  
		t = T; c = 0; 
		ybar->num_bouts[beh] = 0;
		while(t >= 0 && states[t][c].start_frame >= 0) { // start at T, determine label and start_frame, check at states[start_frame], ...
			tt = states[t][c].start_frame; 
			c = states[t][c].behavior; 
			t = tt;
			ybar->num_bouts[beh]++; 
		}
		ybar->bouts[beh] = (BehaviorBout*)malloc(sizeof(BehaviorBout)*ybar->num_bouts[beh]);

		// Backtrack one more time to actually store that solution
		t = T; c = 0; 
		ybar->scores[beh] = 0; 
		//ybar->losses[beh] = max_fn_cost;
		ybar->losses[beh] = 0;
		i = ybar->num_bouts[beh]-1;
		while(t >= 0 && states[t][c].start_frame >= 0) { // start at T, determine label and start_frame, check at states[start_frame], ...
			ybar->scores[beh] += states[t][c].bout_score + states[t][c].transition_score;
			ybar->losses[beh] += states[t][c].loss_fn + states[t][c].loss_fp;
			ybar->bouts[beh][i] = states[t][c];
			tt = states[t][c].start_frame; 
			c = states[t][c].behavior; 
			t = tt; 
			i--; 
		}
		ybar->score += ybar->scores[beh];
		assert(my_abs(ybar->scores[beh] + ybar->losses[beh] - table[T][0]) <= .01); // sanity check for dynamic programming

		if(b->partial_label) {
			// Restore modified transition tables, if necessary
			for(c = 0; c < num_classes[beh]; c++) {
				class_training_transitions_count[beh][c] = old_class_transition_counts[c];
				class_training_count[beh][c] = old_class_training_counts[c];
			}
			free(old_class_transition_counts);
		}

		// By subtracting off the score for the groundtruth label, 
		// the score will be positive for a violated constraint
		if(y) { // if invoked in LEARNING mode
			// this block contains either debug output or sanity checks only 

			ybar->slack += ybar->scores[beh] + ybar->losses[beh];
			ybar->loss += ybar->losses[beh];
#if DEBUG > 0
			// Makes sure the computed components of the loss aggregated during dynamic programming are identical
			// to the loss when comparing the sequences y and ybar.  If this fails, something is probably wrong
			// with the computation of the loss during dynamic programming or with the function loss2()
			double l = loss2(y_gt, y_bar, beh, 1);
			assert(my_abs(ybar->losses[beh] - l) < .01);
			SparseVector ybar_psi = Psi(b, ybar);
			double real_score = w->dot(ybar_psi);
			if(my_abs(ybar->score - real_score) >= .1){
				printf("trouble...");
// 				double ybar_bout_score = 0;
// 				double ybar_trans_score = 0;
// 				for(int ibout=0; ibout < ybar->num_bouts[0]; ibout++) {
// 					ybar_bout_score += ybar->bouts[0][ibout].bout_score;
// 					ybar_trans_score += ybar->bouts[0][ibout].transition_score;
// 				}
// 				double ybar_bout_score_real = 0;
// 				double ybar_trans_score_real = 0;
// 				double* w_ns = w->non_sparse<double>(1830);
// 				double* psi_ns = ybar_psi.non_sparse<double>(1830);
// 				for(int i_w=0; i_w<1800; i_w++){
// 					ybar_bout_score_real += w_ns[i_w] * psi_ns[i_w];
// 				}
// 				for(int i_w=1800; i_w<1825; i_w++){
// 					ybar_trans_score_real += w_ns[i_w] * psi_ns[i_w];
// 				}
// 				int bla = 0;
			}
			//assert(my_abs(ybar->score - real_score) < 1);
#endif

			// Compute the scores for the ground truth label y.  The cache tables in dynamic programming should always have
			// at least as high a score as the score yielded by the ground truth labeling.  If this fails, something
			// is probably wrong with the dynamic programming algorithm
			y->scores[beh] = 0;
			j = 1;
			for(i = 0; i < beh; i++) 
				j += num_classes[i]*(num_classes[i]+num_features);
			for(i = 0; i < y->num_bouts[beh]; i++) {
				psi_bout(b, y->bouts[beh][i].start_frame, y->bouts[beh][i].end_frame, beh, -1, tmp_features, true, false);  
				y->bouts[beh][i].bout_score = 0;
				for(int k = 0; k < num_features; k++) 
				  y->bouts[beh][i].bout_score += class_weights[y->bouts[beh][i].behavior][k]*tmp_features[k];
// #ifdef ALLOW_SAME_TRANSITIONS
// 				getTransitionScore(y->bouts[beh][i].transition_score, unary_costs[y->bouts[beh][i].behavior], i < y->num_bouts[beh]-1, transition_weights[y->bouts[beh][i].behavior][y->bouts[beh][i+1].behavior]); // unary costs added, same-transition costs added whenever same transition     is used (y->bouts[beh][i].behavior and y->bouts[beh][i+1].behavior MAY BE the SAME behavior)
// #else
				getTransitionScore(y->bouts[beh][i].transition_score, transition_weights[y->bouts[beh][i].behavior][y->bouts[beh][i].behavior], i < y->num_bouts[beh]-1, transition_weights[y->bouts[beh][i].behavior][y->bouts[beh][i+1].behavior]); // unary costs stored in same transition costs;            same transitions never used (y->bouts[beh][i].behavior and y->bouts[beh][i+1].behavior are ALWAYS DIFFERENT behaviors)
//#endif
//				y->bouts[beh][i].transition_score = transition_weights[y->bouts[beh][i].behavior][y->bouts[beh][i].behavior];
//				if(i < y->num_bouts[beh]-1) 
//					y->bouts[beh][i].transition_score += transition_weights[y->bouts[beh][i].behavior][y->bouts[beh][i+1].behavior];
				y->bouts[beh][i].loss_fn = y->bouts[beh][i].loss_fp = 0;
				y->scores[beh] += y->bouts[beh][i].bout_score + y->bouts[beh][i].transition_score;

#if DEBUG > 0
				// Making sure that ybar_score + ybar_loss >= y_score, so far
				if(y->scores[beh] > (i < y->num_bouts[beh]-1 ? table[y->bouts[beh][i+1].start_frame][y->bouts[beh][i+1].behavior] : table[T][0])) {
					// Something went wrong, it might be informative for debugging (break here using gdb) to test the same dynamic 
					// programming problem but without using loss, then compare y_max to y
					g_table = table; g_states = states; g_y = y;
					assert(0);
				}
#endif
			}
			y->score += y->scores[beh];
#if DEBUG > 0
			if(ybar->score+ybar->loss < y->score) {
				printf("trouble...");
			}
			assert(ybar->score+ybar->loss >= y->score);
#endif
		}

#if DEBUG > 1
								for(cc=0; cc < numA; cc++) {
									fclose(outA[cc]);
								}
								free(outA);
								free(writeHeader);
#endif

		free(class_weights);
/*
#ifdef ALLOW_SAME_TRANSITIONS
		free(unary_costs);
#endif
*/
		free(table);
		free(states);
		if(fn) free(fn);
	} // for(beh = 0; beh < behaviors->num; beh++) {


#if DEBUG > 0
	SparseVector *fvec;
	//LABEL y_max;
	double y_score;
	if(y) {
		// Compare the score computed for y when looping through the bouts in y to the score when evaluating
		// w*psi(x,y).  If this fails, something is probably wrong with the dynamic programming algorithm with
		// respect to the extraction of features or bout scores, or something is wrong with the function psi()
		if(!b->fvec) b->fvec = Psi(x, y).ptr();
		y_score = w->dot(*b->fvec);
		assert(my_abs(y_score - y->score) < 0.1);

		// Further sanity checks that are probably redundant with earlier ones
		ybar->slack -= y->score;
		//assert(ybar->slack >= -0.01);
		assert(ybar->slack >= -1);
		//y_max = classify_struct_example(*x, sm, sparm);
		//assert(((BehaviorBoutSequence*)y_max.data)->score >= y->score);
		fvec = Psi(x, y_bar).ptr();
		double ybar_score = w->dot(*fvec);
		assert(my_abs(ybar->score - ybar_score) < 1);
		delete fvec;

#if DEBUG > 2
		printf("w:"); for(i = 0; i < sizePsi; i++) { printf(" %d:%f", i, (float)ww[i]); } printf("\n");
		char str[100000]; b->fvec->to_string(str); printf("psi: %s\n", str); 
#endif

#if DEBUG > 1
		Json::Value r = y->save(this); printf("y: %s\n", writer.write(r).c_str()); 
		r = ybar->save(this); printf("ybar: %s\n\n", writer.write(r).c_str()); 
#endif
	} 
#endif

	free(tmp_features);
	free(tmp_features_);
	free(durations);
	free(ww);

	return ybar->score + ybar->loss;
}


double      SVMBehaviorSequence::Loss(StructuredLabel *y_gt,  StructuredLabel *y_pred) {
	int i;
	double l = 0;

	if(behavior >= 0) 
	  return loss2(y_gt, y_pred, behavior, 0);
	else {
		for(i = 0; i < behaviors->num; i++)
			l += loss2(y_gt, y_pred, i, 0);
	}
	return l;
}

 double      SVMBehaviorSequence::loss2(StructuredLabel *y_gt,  StructuredLabel *y_pred, int beh, int debug)
{
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
	double cl, l_fn = 0, l_fn2 = 0;

// 	// TEMP EYRUN
// 	l = 0;
// 	for(curr_ybar = 0; curr_ybar < ybar->num_bouts[0]; curr_ybar++)
// 		l += ybar->bouts[0][curr_ybar].loss_fp + ybar->bouts[0][curr_ybar].loss_fn;
// 	return l;

	while(curr_y < y->num_bouts[beh] || curr_ybar < ybar->num_bouts[beh]) {
		if(curr_y < y->num_bouts[beh] && curr_ybar < ybar->num_bouts[beh]) {
			// Check if the current bout in y and ybar match.  
			dur_ybar = (b->frame_times[my_min(ybar->bouts[beh][curr_ybar].end_frame,T-1)] - 
				b->frame_times[ybar->bouts[beh][curr_ybar].start_frame]);
			dur_y = (b->frame_times[my_min(y->bouts[beh][curr_y].end_frame,T-1)] -
				b->frame_times[y->bouts[beh][curr_y].start_frame]);
			inter = (b->frame_times[my_min(my_min(y->bouts[beh][curr_y].end_frame, ybar->bouts[beh][curr_ybar].end_frame),T-1)] - 
				b->frame_times[my_max(y->bouts[beh][curr_y].start_frame,ybar->bouts[beh][curr_ybar].start_frame)]);
			assert(inter >= 0);
			if(y->bouts[beh][curr_y].behavior == ybar->bouts[beh][curr_ybar].behavior) {
				sum_ybar += inter;
				sum_y += inter;
				//if(debug)
				l_fn -= match_false_negative_cost(dur_y, beh, y->bouts[beh][curr_y].behavior)*(dur_y ? (inter/dur_y) : 0);
			} else
				l_fn2 += match_false_negative_cost(dur_y, beh, y->bouts[beh][curr_y].behavior)*(dur_y ? (inter/dur_y) : 0);
		} else
			inter = 0;

		if(curr_ybar == ybar->num_bouts[beh] || (curr_y < y->num_bouts[beh] && 
			(y->bouts[beh][curr_y].end_frame <= ybar->bouts[beh][curr_ybar].end_frame))) {
				// Go to the next bout in y adding the appropriate loss based on whether or not
				// the bout was matched to a bout in ybar
				cl = match_false_negative_cost(dur_y, beh, y->bouts[beh][curr_y].behavior)*(dur_y ? ((dur_y-sum_y)/dur_y) : 1);
				l += cl; 
				if(debug == 2) fprintf(stderr, "y[%d]->%f %f\n", curr_y, cl, l);
				curr_y++;
				sum_y = 0;
		} else {
			// Go to the next bout in ybar, adding the appropriate loss based on whether or not
			// the bout was matched to a bout in y
			cl = match_false_positive_cost(dur_ybar, beh, ybar->bouts[beh][curr_ybar].behavior) *
				(dur_ybar ? ((dur_ybar-sum_ybar)/dur_ybar) : 1);
			l += cl;
			if(debug == 2) 
				fprintf(stderr, "ybar[%d]->%f %f\n", curr_ybar, cl, l);
			if(debug) {
				assert(my_abs(cl - ybar->bouts[beh][curr_ybar].loss_fp) < .00001);
				assert(my_abs(l_fn2 - ybar->bouts[beh][curr_ybar].loss_fn) < .00001);
			}
			ybar->bouts[beh][curr_ybar].loss_fn = l_fn2;
			ybar->bouts[beh][curr_ybar].loss_fp = cl;
			curr_ybar++;
			sum_ybar = 0;
			l_fn = 0; l_fn2 = 0;
		}
	}
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
		}
	}

	int num_histogram_bins = 0, num_max_thresholds = 0, num_min_thresholds = 0;
	Json::Value params;
	for(i = 0; i < num_base_features; i++) {
		SVMFeatureParams *p = &feature_params[i];
		Json::Value c = root["feature_params"][i];
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

	for(beh = 0; beh < behaviors->num; beh++) {
	  compactUnaryCosts = (sizePsi) == compactUnarySize(num_features, num_classes[beh]); /* SJB: removed -1.  CSC: true iff unary costs are stored as same transition costs (and no further weights were given) */
	  extraUnaryCosts = (sizePsi) == extraUnarySize(num_features, num_classes[beh]); /* SJB: removed -1.  CSC: true iff case Unary costs are given in addition to same transition costs; equals !compactUnaryCosts but explicitely given for sanity checks (CHECK THE ASSIGNMENT WHENEVER ADDITIONAL WEIGHTS (other than feature & transition weights) ARE ADDED) */

	  assert(compactUnaryCosts != extraUnaryCosts);
	}

	return true;
}



BehaviorBoutFeatures::BehaviorBoutFeatures() {
  partial_label = NULL;
  memory_buffer = NULL;
  fvec = NULL;
}

BehaviorBoutFeatures::~BehaviorBoutFeatures() {
  if(features) free(features);
  if(fvec) delete fvec;
  if(partial_label) delete partial_label;
}

/*
* Allocate space for feature caches used to compute bout-level features efficiently
*/
void BehaviorBoutFeatures::AllocateBuffers(SVMBehaviorSequence *svm) {
	int i;
	int T = num_frames;
	long int cache_features_size = 
		num_base_features*sizeof(double*) +        // features
		num_base_features*T*sizeof(double) +       // features[i]
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
			//sizeof(double*)*p->num_histogram_bins + //integral_histogram_features[i]
			//p->num_histogram_bins*(T+1)*3*sizeof(double)+10000; //integral_histogram_features[i][j]
			num_base_features*sizeof(double*)*p->num_histogram_bins + //integral_histogram_features[i]
			num_base_features*p->num_histogram_bins*(T+1)*3*sizeof(double)+10000; //integral_histogram_features[i][j]
	}
	this->memory_buffer = (unsigned char*)malloc(cache_features_size);

	// Extract regular raw features
	this->features = (double**)this->memory_buffer;
	this->memory_buffer += num_base_features*sizeof(double*);
	for(i = 0; i < num_base_features; i++) {
		this->features[i] = ((double*)this->memory_buffer);
		this->memory_buffer += T*sizeof(double);
	}
	this->frame_times = (double*)this->memory_buffer;
	this->memory_buffer += T*sizeof(double);

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

	memory_buffer = ptr;
}

bool BehaviorBoutFeatures::load(const Json::Value &r, StructuredSVM *s) {
  strcpy(fname, r.get("fname", "").asString().c_str());
  return load(fname, (SVMBehaviorSequence*)s, NULL);
}

Json::Value BehaviorBoutFeatures::save(StructuredSVM *s) {
  Json::Value r;
  r["fname"] = fname;
  return r;
}


 BehaviorBoutSequence::BehaviorBoutSequence(BehaviorBoutFeatures *x, SVMBehaviorSequence *svm) : StructuredLabel(x) {
  this->behaviors = svm->behaviors;
  features = x;
  num_bouts = NULL;
  bouts = NULL;
  score = loss = slack = 0;
  scores = losses = NULL;
}

BehaviorBoutSequence::~BehaviorBoutSequence() {
	for(int i = 0; i < behaviors->num; i++)
		if(bouts[i])
			free(bouts[i]);
	free(bouts);
}


Json::Value BehaviorBoutSequence::save(StructuredSVM *s) {
  Json::Value r;
  r["fname"] = fname;
  save(fname);
  return r;
}

bool BehaviorBoutSequence::load(const Json::Value &r, StructuredSVM *s) {
  strcpy(fname, r.get("fname", "").asString().c_str());
  return load(fname);
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
		char *html_tmp=(char*)malloc(10000000), folder[1000], fname2[1000];
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
				sprintf(alt, "behavior=%s bout_score=%f transition_score=%f loss_fn=%f loss_fp=%f", group->values[c].name, 
					(float)this->bouts[beh][i].bout_score, (float)this->bouts[beh][i].transition_score, (float)this->bouts[beh][i].loss_fn,  (float)this->bouts[beh][i].loss_fp);
				sprintf(str, "  <area shape=\"rect\" coords=\"%d,%d,%d,%d\" href=\"features_%s.html#%d\" title=\"%s\" />\n", (int)(z*this->bouts[beh][i].start_frame), 0, 
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


