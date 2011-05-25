/***********************************************************************/
/*                                                                     */
/*   svm_struct_api.c                                                  */
/*                                                                     */
/*   Definition of API for attaching implementing SVM learning of      */
/*   structures (e.g. parsing, multi-label classification, HMM)        */ 
/*                                                                     */
/*   Author: Thorsten Joachims                                         */
/*   Date: 03.07.04                                                    */
/*                                                                     */
/*   Copyright (c) 2004  Thorsten Joachims - All rights reserved       */
/*                                                                     */
/*   This software is available for non-commercial use only. It must   */
/*   not be modified and distributed without prior permission of the   */
/*   author. The author is not responsible for implications from the   */
/*   use of this software.                                             */
/*                                                                     */
/***********************************************************************/

#include <stdio.h>
#include <string.h>
#include "svm_struct/svm_struct_common.h"
#include "svm_struct_api_behavior_sequence.h"
#include "svm_struct/svm_struct_learn.h"

#ifdef DEBUG > 0 
char *g_currFile; // CSC 20110420: hack to pass current filename for debug purposes
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


#include "cv.h"
#include "highgui.h"
IplImage *VisualizeBouts(BehaviorBoutSequence *seq, BehaviorGroups *groups, int beh, const char *fname, char *html);

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


SVMBehaviorSequence::SVMBehaviorSequence(int num_feat, struct _BehaviorGroups *behaviors, int beh, SVMFeatureParams *sparams) {
	for(int i = 0; i < behaviors->num; i++) {
		if(beh == -1 || i == beh) {
			restrict_behavior_features[i] = (bool**)my_malloc(behaviors->behaviors[i].num_values*sizeof(bool*));
			memset(restrict_behavior_features[i], 0, sizeof(bool*)*behaviors->behaviors[i].num_values);
			/*
			CSC & SB 20110324: changed beh to i
			*/
		}
	}

	Init(num_feat, behaviors, beh, sparams);
}

SVMBehaviorSequence::SVMBehaviorSequence(struct _BehaviorGroups *behaviors, int beh) {
	Init(0, behaviors, beh);
	for(int i = 0; i < behaviors->num; i++) {
		if(beh == -1 || i == beh) {
			restrict_behavior_features[i] = (bool**)my_malloc(behaviors->behaviors[i].num_values*sizeof(bool*));
			memset(restrict_behavior_features[i], 0, sizeof(bool*)*behaviors->behaviors[i].num_values);
			/*
			CSC & SB 20110324: changed beh to i
			restrict_behavior_features[beh] = (bool**)my_malloc(behaviors->behaviors[beh].num_values*sizeof(bool*));
			memset(restrict_behavior_features[beh], 0, sizeof(bool*)*behaviors->behaviors[beh].num_values);
			*/
		}
	}
}

void SVMBehaviorSequence::Init(int num_feat, struct _BehaviorGroups *behaviors, int beh, SVMFeatureParams *sparams) {
	this->behaviors = behaviors;
	this->behavior = beh;
	this->feature_diff_frames = feature_diff_frames;
#ifdef ALLOW_SAME_TRANSITIONS
	time_approximation = .05;
#else
	time_approximation = 0; // CSC: test EVERY Frame, no log time-saving unless optimality may still be guaranteed!
#endif
	class_training_transitions = NULL;
	class_training_transitions_count = class_training_count = NULL;
	features_mu = features_gamma = NULL;
	histogram_thresholds = min_thresholds = max_thresholds = NULL;
	min_bout_duration = 1;
	feature_names = NULL;
	strcpy(load_const_set_fname, "");
	strcpy(save_const_set_fname, "");
	strcpy(load_constraints_fname, "");
	strcpy(save_constraints_fname, "");

	int i;
	num_base_features = num_feat;
	for(i = 0; i < num_base_features; i++) {
		feature_params[i] = sparams ? sparams[i] : DefaultParams();
		assert(feature_params[i].num_temporal_levels <= MAX_TEMPORAL_LEVELS && 
			feature_params[i].num_histogram_temporal_levels <= MAX_TEMPORAL_LEVELS &&
			feature_params[i].num_difference_temporal_levels <= MAX_TEMPORAL_LEVELS &&
			feature_params[i].num_harmonic_features <= MAX_HARMONIC_LEVELS);
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
* Main training function. Simply calls other routines to read in training parameters, invoke the main
* training algorithm, and save the learned model to disk
*/
int SVMBehaviorSequence::train (int argc, const char* argv[], STRUCT_LEARN_PARM *struct_parm, STRUCTMODEL *structmodel)
{  
	train_main(argc, argv, struct_parm, structmodel, this);

	return 0;
}

int SVMBehaviorSequence::test (int argc, const char* argv[], STRUCT_LEARN_PARM *struct_parm, STRUCTMODEL *structmodel)
{  
	test_main(argc, argv, struct_parm, structmodel, this);

	return 0;
}




#define compactUnarySize(p_features, p_classes) ((p_features) + (p_classes)    ) * (p_classes)
#define   extraUnarySize(p_features, p_classes) ((p_features) + (p_classes) + 1) * (p_classes)
#ifdef ALLOW_SAME_TRANSITION
#define getPsiSize(p_features, p_classes) compactUnarySize(p_features, p_classes)
#else
#define getPsiSize(p_features, p_classes) extraUnarySize(p_features, p_classes)
#endif



/* Called in learning part before anything else is done to allow
*  any initializations that might be necessary. 
* Sets a constant loss term for a bout false positive and false negative
*/
void SVMBehaviorSequence::svm_struct_learn_api_init(int argc, const char* argv[]) {
	int i, beh;

	sizePsi = 0;
	compute_feature_space_size();

	for(beh = 0; beh < behaviors->num; beh++) {
		num_classes[beh] = behaviors->behaviors[beh].num_values;

		false_negative_cost[beh] = (double*)malloc(2*sizeof(double)*num_classes[beh]);
		false_positive_cost[beh] = false_negative_cost[beh]+num_classes[beh];

		// It is assumed here that the label 0 is the unknown label.  The user could define custom
		// class-specific values for these constants .
		// In the future, we could implement this with a per class-pair confusion cost
		false_negative_cost[beh][0] = 0;
		false_positive_cost[beh][0] = 100; // CSC: Move FP/FN costs to Params/BehaviorParams.txt? or to Model/Model.txt?
		for(i = 1; i < num_classes[beh]; i++) {
			false_negative_cost[beh][i] = 100;
			false_positive_cost[beh][i] = 100;
		}

//		sizePsi += (this->num_features+this->num_classes[beh])*this->num_classes[beh];
		sizePsi += getPsiSize(this->num_features, this->num_classes[beh]);
	}
}

/* Called in learning part at the very end to allow any clean-up
*  that might be necessary. 
*/
void        SVMBehaviorSequence::svm_struct_learn_api_exit()
{
}

/* Called in prediction part before anything else is done to allow
* any initializations that might be necessary. 
*/
void        SVMBehaviorSequence::svm_struct_classify_api_init(int argc, const char* argv[])
{
	//compute_feature_space_size();
	svm_struct_learn_api_init(argc, argv);
}

/* Called in prediction part at the very end to allow any clean-up
*  that might be necessary. 
*/
void        SVMBehaviorSequence::svm_struct_classify_api_exit()
{
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

/*
* Compute mean, variance and median statistics on the base features over the training set.
* the mean and variance are used to normalize each feature.  The median statistics are used to
* construct histogram and threshold constants, by placing some fraction of the training set 
* between each threshold
*/
void SVMBehaviorSequence::compute_feature_mean_variance_median_statistics(EXAMPLE *ex, int num_examples) {
	int i, j, k, n, ind, beh, num_bouts = 0, curr = 0;
	double **feat, *ptr, mean, num_histogram_bins = 0, num_max_thresholds = 0, num_min_thresholds = 0, w, target, m;
	BehaviorBoutSequence *y;
	double *tmp_features = (double*)malloc(sizeof(double)*num_features);
	SVMFeatureParams *p;
	BehaviorBoutFeatures *x;


	// Count the total number of bouts in the training set
	for(n = 0; n < num_examples; n++)
		for(beh = 0; beh < behaviors->num; beh++)
			if(behavior < 0 || beh == behavior) 
				num_bouts += ((BehaviorBoutSequence*)ex[n].y.data)->num_bouts[beh];

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
		for(n = 0; n < num_examples; n++) {
			y = ((BehaviorBoutSequence*)ex[n].y.data);
			x = ((BehaviorBoutFeatures*)ex[n].x.data);
			for(beh = 0; beh < behaviors->num; beh++) {
				if(behavior < 0 || beh == behavior) {
					for(j = 0; j < y->num_bouts[beh]; j++) {
						mean = 0;
						for(k = y->bouts[beh][j].start_frame; k < y->bouts[beh][j].end_frame; k++) {
							assert(k >= 0);  // CSC 20110324
							mean += x->features[i][k];
						}
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
		}
	}

	for(i = 0; i < num_base_features; i++) {
		// Compute the bout-level min over all training examples
		curr = 0;
		for(n = 0; n < num_examples; n++) {
			y = ((BehaviorBoutSequence*)ex[n].y.data);
			x = ((BehaviorBoutFeatures*)ex[n].x.data);
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
		for(n = 0; n < num_examples; n++) {
			y = ((BehaviorBoutSequence*)ex[n].y.data);
			x = ((BehaviorBoutFeatures*)ex[n].x.data);
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
	for(n = 0; n < num_examples; n++) {
		x = (BehaviorBoutFeatures*)ex[n].x.data;
		load_behavior_bout_features(x->data, x);
		compute_bout_feature_caches(x);
	}

	// Now that we have bout features, we can compute the mean of each bout-wise feature over the training set
	features_mu = (double*)malloc(num_features*sizeof(double)*2);
	features_gamma = features_mu + num_features;
	for(i = 0; i < num_features; i++) 
		features_mu[i] = features_gamma[i] = 0;
	for(n = 0; n < num_examples; n++) {
		y = ((BehaviorBoutSequence*)ex[n].y.data);
		x = ((BehaviorBoutFeatures*)ex[n].x.data);
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
	for(n = 0; n < num_examples; n++) {
		y = ((BehaviorBoutSequence*)ex[n].y.data);
		x = ((BehaviorBoutFeatures*)ex[n].x.data);
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

	free(feat);
	free(tmp_features);
}

#if DEBUG > 0
#define MAX_FEATURES 1000
int g_feature_map[MAX_FEATURES] = {0}; // if first element is initialized the rest will be zero-filled
int g_feature_map2[MAX_FEATURES] = {0}; // if first element is initialized the rest will be zero-filled
int g_feature_map3[MAX_FEATURES] = {0}; // if first element is initialized the rest will be zero-filled
int g_csc_ind;

#define ADD_FEATURE(feat, ind, val, feature_id, feature2_id, feature3_id) {\
	g_csc_ind = ind; \
	(feat)[(g_csc_ind)] = (val); \
	g_feature_map[(g_csc_ind)] = feature_id; \
	g_feature_map2[(g_csc_ind)] = feature2_id;\
	g_feature_map3[(g_csc_ind)] = feature3_id;\
}
#else
#define ADD_FEATURE(feat, ind, val, feature_id)	(feat)[(ind)] = (val);
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
double *SVMBehaviorSequence::psi_bout(BehaviorBoutFeatures *b, int t_start, int t_end, int beh, int c, double *feat, bool normalize, bool fast_update) {
	if(!feat) feat = (double*)my_malloc(num_features*sizeof(double));

	int i, j, k, l, ind = 0, td, beg;
	double start, temporal_grid_size[MAX_TEMPORAL_LEVELS], harmonic_grid_size[MAX_HARMONIC_LEVELS+1], sum, f, t, ave, sum_sqr;
	SVMFeatureParams *p;
	double dur = t_end-t_start;
	double inv = dur ? 1/dur : 0;

	if(!fast_update)
		b->bout_start = b->bout_end = -1;
	else
		assert(b->bout_end == t_end);

	update_bout_feature_caches(b, t_start, t_end, c);

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
				if(p->use_bout_sum_features) ADD_FEATURE(feat, ind++, f, USE_BOUT_SUM_FEATURES, i, j);
				if(p->use_bout_ave_features) ADD_FEATURE(feat, ind++, f*inv, USE_BOUT_AVE_FEATURES, i, j);
				if(p->use_bout_sum_absolute_features) ADD_FEATURE(feat, ind++, my_abs(f), USE_BOUT_SUM_ABSOLUTE_FEATURES, i, j);
				if(p->use_bout_ave_absolute_features) ADD_FEATURE(feat, ind++, my_abs(f)*inv, USE_BOUT_AVE_ABSOLUTE_FEATURES, i, j);
			}
		}

		// The standard deviation of the feature in the interval (t_start,t_end),  Can also use the total sum variance
		if(p->use_standard_deviation || p->use_sum_variance) {
			sum_sqr = (b->integral_sqr_features[i][t_end]-b->integral_sqr_features[i][t_start]);
			f = sum_sqr - 2*ave*sum + SQR(ave)*dur;
			if(f < 0) f = 0;  // avoid precision-related errors
			if(p->use_sum_variance) ADD_FEATURE(feat, ind++, f, USE_SUM_VARIANCE, i, 0);
			if(p->use_standard_deviation) ADD_FEATURE(feat, ind++, sqrt(f*inv), USE_STANDARD_DEVIATION, i, 0);
		}

		// Add features for the min/max feature response in the current bout, or thresholded versions of the min/max values
		if(p->use_bout_max_feature) 
			ADD_FEATURE(feat, ind++, b->bout_max_feature_responses[i], USE_BOUT_MAX_FEATURE, i, j);
		if(p->use_bout_min_feature) 
			ADD_FEATURE(feat, ind++, b->bout_min_feature_responses[i], USE_BOUT_MIN_FEATURE, i, j);
		for(l = 0; l < p->num_bout_max_thresholds; l++) 
			ADD_FEATURE(feat, ind++, b->bout_max_feature_responses[i] > max_thresholds[i][l] ? 1 : 0, NUM_BOUT_MAX_THRESHOLDS, i, j);
		for(l = 0; l < p->num_bout_min_thresholds; l++) 
			ADD_FEATURE(feat, ind++, b->bout_min_feature_responses[i] < min_thresholds[i][l] ? 1 : 0, NUM_BOUT_MIN_THRESHOLDS, i, j);


		// These are differences in the bout average or sum response for the bout as compared to the global max, min,
		// and average over the entire behavioral sequence
		if(p->use_global_difference_max_sum_features) 
			ADD_FEATURE(feat, ind++, sum - b->max_feature_responses[i]*(t_end-t_start), USE_GLOBAL_DIFFERENCE_MAX_SUM_FEATURES, i, 0);
		if(p->use_global_difference_max_ave_features) 
			ADD_FEATURE(feat, ind++, ave - b->max_feature_responses[i], USE_GLOBAL_DIFFERENCE_MAX_AVE_FEATURES, i, 0);
		if(p->use_global_difference_min_sum_features) 
			ADD_FEATURE(feat, ind++, sum - b->min_feature_responses[i]*(t_end-t_start), USE_GLOBAL_DIFFERENCE_MIN_SUM_FEATURES, i, 0);
		if(p->use_global_difference_min_ave_features) 
			ADD_FEATURE(feat, ind++, ave - b->min_feature_responses[i], USE_GLOBAL_DIFFERENCE_MIN_AVE_FEATURES, i, 0);
		if(p->use_global_difference_ave_sum_features) 
			ADD_FEATURE(feat, ind++, sum - b->ave_feature_responses[i]*(t_end-t_start), USE_GLOBAL_DIFFERENCE_AVE_SUM_FEATURES, i, 0);
		if(p->use_global_difference_ave_ave_features) 
			ADD_FEATURE(feat, ind++, ave - b->ave_feature_responses[i], USE_GLOBAL_DIFFERENCE_AVE_AVE_FEATURES, i, 0);


		// The total change from the beginning of the bout to the end of the bout
		if(p->num_bout_change_points) {
			for(j = 0, t=t_start; j < p->num_bout_change_points; j++, t+=harmonic_grid_size[p->num_bout_change_points]) {
				f = b->features[i][my_max((int)t,t_end-1)] - b->features[i][t_start];
				if(p->use_bout_change) 
					ADD_FEATURE(feat, ind++, f, USE_BOUT_CHANGE, i, j);
				if(p->use_bout_absolute_change) 
					ADD_FEATURE(feat, ind++, my_abs(f), USE_BOUT_ABSOLUTE_CHANGE, i, j);
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
					if(p->use_histogram_sum_features) ADD_FEATURE(feat, ind++, f, USE_HISTOGRAM_SUM_FEATURES, i, j);
					if(p->use_histogram_ave_features) ADD_FEATURE(feat, ind++, f*inv, USE_HISTOGRAM_AVE_FEATURES, i, j);
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
				ADD_FEATURE(feat, ind++, f, USE_START_SUM_DIFF_HAAR_FEATURES, i, 0);
			if(p->use_start_ave_diff_haar_features)
				ADD_FEATURE(feat, ind++, f*inv, USE_START_AVE_DIFF_HAAR_FEATURES, i, 0);
			if(p->use_start_sum_absolute_diff_haar_features) 
				ADD_FEATURE(feat, ind++, my_abs(f), USE_START_SUM_ABSOLUTE_DIFF_HAAR_FEATURES, i, 0);
			if(p->use_start_ave_absolute_diff_haar_features) 
				ADD_FEATURE(feat, ind++, my_abs(f*inv), USE_START_AVE_ABSOLUTE_DIFF_HAAR_FEATURES, i, 0);


			// Difference in average feature response in the region inside the bout (t_start,t_end) and 
			// the region of duration (t_end-t_start)/(2^j) immediately after t_end
			f = -b->integral_features[i][t_end+td] + 2*b->integral_features[i][t_end] - b->integral_features[i][t_end-td];
			if(p->use_end_sum_diff_haar_features)
				ADD_FEATURE(feat, ind++, f, USE_END_SUM_DIFF_HAAR_FEATURES, i, 0);
			if(p->use_end_ave_diff_haar_features)
				ADD_FEATURE(feat, ind++, f*inv, USE_END_AVE_DIFF_HAAR_FEATURES, i, 0);
				//ADD_FEATURE(feat, ind++, f*inv, USE_END_AVE_DIFF_HAAR_FEATURES);
			if(p->use_end_sum_absolute_diff_haar_features) 
				ADD_FEATURE(feat, ind++, my_abs(f), USE_END_SUM_ABSOLUTE_DIFF_HAAR_FEATURES, i, 0);
			if(p->use_end_ave_absolute_diff_haar_features) 
				ADD_FEATURE(feat, ind++, my_abs(f*inv), USE_END_AVE_ABSOLUTE_DIFF_HAAR_FEATURES, i, 0);
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
				ADD_FEATURE(feat, ind++, f, USE_SUM_HARMONIC_FEATURES, i, j);
			if(p->use_ave_harmonic_features)
				ADD_FEATURE(feat, ind++, f*inv, USE_AVE_HARMONIC_FEATURES, i, j);
			if(p->use_sum_absolute_harmonic_features)
				ADD_FEATURE(feat, ind++, my_abs(f), USE_SUM_ABSOLUTE_HARMONIC_FEATURES, i, j);
			if(p->use_ave_absolute_harmonic_features)
				ADD_FEATURE(feat, ind++, my_abs(f*inv), USE_AVE_ABSOLUTE_HARMONIC_FEATURES, i, j);
		}
		assert(ind-beg == p->num_features);
	}
	assert(ind == num_features);

	// Normalize each feature to have mean 0 and standard deviation 1 (as determined by the training set)
	if(normalize) {
		for(i = 0; i < ind; i++)  {
			if(isnan(feat[i])) {
				printf("Upcoming error...\n");
			}
			assert(!isnan(feat[i]));
			feat[i] = (feat[i]-features_mu[i])*features_gamma[i];
			assert(!isnan(feat[i]));
		}
	}

	// For a given behavior class, optionally use only a subset of the available bout features
	if(beh >= 0 && c >= 0 && restrict_behavior_features[beh] && restrict_behavior_features[beh][c]) {
		for(i = 0; i < ind; i++) 
			if(!restrict_behavior_features[beh][c][i])
				feat[i] = 0;
	}

	return feat;
}


/*
* This is called to update the bout feature cache structures (like the bout-wise min and max)
* in a computationally efficient way, e.g. when doing dynamic programming we can update the
* min/max from the bout explored in the previous loop iteration in constant time (as opposed
* to time proportional to the number of frames in the bout)
*/
void SVMBehaviorSequence::update_bout_feature_caches(BehaviorBoutFeatures *b, int t_start, int t_end, int c) {
	/*
	BehaviorBoutFeatures *b: contains b-> start and b->and, which contain

	Two usage scenarios:
	(a) ev
	(b) within 
	*/
	bool reset = false;
	int i, j;
	SVMFeatureParams *p;

	if(t_end == b->bout_end) {
		if(t_start >= b->bout_start) {
			printf("Upcoming Error...\n");
		}
		assert(t_start < b->bout_start);
	} else {
		reset = true;
		b->bout_start = b->bout_end = t_end;
	}
	for(i = 0; i < num_base_features; i++) {
		p = &feature_params[i];
		if(reset) {
			b->bout_max_feature_responses[i] = -INFINITY;
			b->bout_min_feature_responses[i] = INFINITY;
		}    
		for(j = t_start; j < b->bout_start; j++) {
			b->bout_max_feature_responses[i] = my_max(b->smoothed_features[i][j],  b->bout_max_feature_responses[i]);
			b->bout_min_feature_responses[i] = my_min(b->smoothed_features[i][j],  b->bout_min_feature_responses[i]);
		}
	}
	b->bout_start = t_start;
}

SVECTOR     *SVMBehaviorSequence::psi(SPATTERN x, LABEL yy, STRUCTMODEL *sm,
	STRUCT_LEARN_PARM *sparm)
{
	/* Returns a feature vector describing the match between pattern x
	and label y. The feature vector is returned as a list of
	SVECTOR's. Each SVECTOR is in a sparse representation of pairs
	<featurenumber:featurevalue>, where the last pair has
	featurenumber 0 as a terminator. Featurenumbers start with 1 and
	end with sizePsi. Featuresnumbers that are not specified default
	to value 0. As mentioned before, psi() actually returns a list of
	SVECTOR's. Each SVECTOR has a field 'factor' and 'next'. 'next'
	specifies the next element in the list, terminated by a NULL
	pointer. The list can be though of as a linear combination of
	vectors, where each vector is weighted by its 'factor'. This
	linear combination of feature vectors is multiplied with the
	learned (kernelized) weight vector to score label y for pattern
	x. Without kernels, there will be one weight in sm.w for each
	feature. Note that psi has to match
	find_most_violated_constraint_???(x, y, sm) and vice versa. In
	particular, find_most_violated_constraint_???(x, y, sm) finds
	that ybar!=y that maximizes psi(x,ybar,sm)*sm.w (where * is the
	inner vector product) and the appropriate function of the
	loss + margin/slack rescaling method. See that paper for details. */
	BehaviorBoutSequence *y = (BehaviorBoutSequence*)yy.data;
	BehaviorBoutFeatures *b = (BehaviorBoutFeatures*)x.data;
	int i, j, beh;
	double *tmp_features = (double*)malloc(sizeof(double)*(sizePsi+num_features+20));
	SVECTOR *retval;
	double *all_features = tmp_features+num_features+10;
	double *ptr = all_features+1, *class_features, *class_transitions;
#ifdef ALLOW_SAME_TRANSITIONS
	double *class_unary;
#endif

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
		class_unary = ptr; ptr += num_classes[beh];                        // line 3: these 3 lines have to correspond to getPsiSize

		for(i = 0; i < num_classes[beh]; i++)
			class_unary[i] = 0;
#endif

		for(i = 0; i < num_classes[beh]*num_features; i++) 
			class_features[i] = 0;

		for(i = 0; i < num_classes[beh]*num_classes[beh]; i++) 
			class_transitions[i] = 0;

		for(i = 0; i < y->num_bouts[beh]; i++) {
/* BEGIN BLOCK */
			// Count of the number transitions from class C1 to C2, for all pairs (C1,C2), C!=C2
			if(i)
				class_transitions[y->bouts[beh][i-1].behavior*num_classes[beh] + y->bouts[beh][i].behavior]++;

#ifdef ALLOW_SAME_TRANSITIONS
			class_unary[y->bouts[beh][i].behavior]++;
			/* Note: Same class transition weights were already computed above and are simply NOT OVERWRITTEN by unary costs */
#ifndef SOLVE_UNARY_COST
				; // transition diagonal still counts same class transitions; transition diagonal not used when classifying with this option
#elif SOLVE_UNARY_COST_CONFLICT == 1
				class_transitions[y->bouts[beh][i].behavior*num_classes[beh] + y->bouts[beh][i].behavior] = 0; // transition diagonal is added to unary costs when classifying with this option
#elif SOLVE_UNARY_COST_CONFLICT == 0
				; // transition diagonal still counts same class transitions; unary costs not used when classifying with this option
#endif
#else // old compact format using diagonal of transition matrix
			// Count of the number of bouts for each class C
			class_transitions[y->bouts[beh][i].behavior*num_classes[beh] + y->bouts[beh][i].behavior]++; // transition diagonal is overwritten by unary costs
#endif
/* END BLOCK */



/* the BLOCK above corresponds to the following lines (copied from above)

#ifdef ALLOW_SAME_TRANSITIONS
		double *unary_costs = (double*)my_malloc(num_classes[beh]*sizeof(double*));
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

	retval = create_svector_n_r(all_features, sizePsi, NULL, 1, 0);
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

void SVMBehaviorSequence::print_features(const char *fname, EXAMPLE *ex, int n, bool normalized) {
	double *tmp_features = (double*)my_malloc(sizeof(double)*num_features);
	int i, j, beh;
	BehaviorBoutSequence *y;
	BehaviorBoutFeatures *x;
	FILE *fout = fopen(fname, "w");
	assert(fout);

	for(i = 0; i < n; i++) {
		y = (BehaviorBoutSequence*)ex[i].y.data;
		x = (BehaviorBoutFeatures*)ex[i].x.data;
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
	for(i = 0; i < num_base_features; i++) {
		p = &feature_params[i];
		curr = num_features;

		for(j = 0; j < p->num_temporal_levels; j++) {
			for(k = 0; k < (1<<j); k++) {
				if(p->use_bout_sum_features) { sprintf(str, "bout_sum_%.3f<=t<%.3f", k/(float)(1<<j), (k+1)/(float)(1<<j)); set_feature_name(num_features++, i, str); }
				if(p->use_bout_ave_features) { sprintf(str, "bout_average_%.3f<=t<%.3f", k/(float)(1<<j), (k+1)/(float)(1<<j)); set_feature_name(num_features++, i, str); }
				if(p->use_bout_sum_absolute_features) { sprintf(str, "bout_sum_abs_%.3f<=t<%.3f", k/(float)(1<<j), (k+1)/(float)(1<<j)); set_feature_name(num_features++, i, str); }
				if(p->use_bout_ave_absolute_features) { sprintf(str, "bout_average_abs_%.3f<=t<%.3f", k/(float)(1<<j), (k+1)/(float)(1<<j)); set_feature_name(num_features++, i, str); }
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

/*
* Precompute certain data structures like integral images, such that bout features can
* be computed more efficiently
*/
void SVMBehaviorSequence::compute_bout_feature_caches(BehaviorBoutFeatures *feature_cache) {
	int i, j, k, T = feature_cache->num_frames;
	double f;
	unsigned char *ptr = feature_cache->memory_buffer;
	int window;

	// integral_features[i] store the running sum value of the i_th feature, such that the sum
	// value of that feature for an arbitrary interval [s,e) can be computed as 
	//   integral_features[e]-integral_features[s]
	// To allow the caller to index into integral_features[i] without worrying about bounds checking,
	// a buffer of size T+1 is added to the beginning and end, which behaves as if the feature in the
	// 0th timestep extends to frame -(T+1) and the (T-1)th pixel extends (T+1) frames into the future
	feature_cache->integral_features = (double**)ptr;
	ptr += num_base_features*sizeof(double*);
	for(i = 0; i < num_base_features; i++) {
		feature_cache->integral_features[i] = ((double*)ptr)+T+1;
		ptr += (T+1)*3*sizeof(double);

		feature_cache->integral_features[i][-T-1] = 0;
		for(j = -T-1; j < 2*(T+1)-1; j++) {
			feature_cache->integral_features[i][j+1] = feature_cache->integral_features[i][j] + 
				feature_cache->features[i][j < 0 ? 0 : (j >= T ? T-1 : j)];
		}
	}

	feature_cache->integral_sqr_features = (double**)ptr;
	ptr += num_base_features*sizeof(double*);
	for(i = 0; i < num_base_features; i++) {
		feature_cache->integral_sqr_features[i] = ((double*)ptr)+T+1;
		ptr += (T+1)*3*sizeof(double);

		feature_cache->integral_sqr_features[i][-T-1] = 0;
		for(j = -T-1; j < 2*(T+1)-1; j++) {
			feature_cache->integral_sqr_features[i][j+1] = feature_cache->integral_sqr_features[i][j] + 
				SQR(feature_cache->features[i][j < 0 ? 0 : (j >= T ? T-1 : j)]);
		}
	}

	// Compute smoothed versions of the raw features
	feature_cache->smoothed_features = (double**)ptr;
	ptr += num_base_features*sizeof(double*);
	for(i = 0; i < num_base_features; i++) {
		feature_cache->smoothed_features[i] = ((double*)ptr);
		ptr += T*sizeof(double);
		window = feature_params[i].feature_sample_smoothness_window;
		for(j = 0; j < T; j++)
			feature_cache->smoothed_features[i][j] = (feature_cache->integral_features[i][j+window+1] - 
			feature_cache->integral_features[i][j-window]) / (2*window+1);
	}

	// ave_feature_responses[i], max_feature_responses[i], and min_feature_responses[i] respectively
	// store the average, maximum, and minimum response of the i_th feature over the entire sequence 1...T
	// The responses at a particular frame are smoothed within a window of size 'window', to add greater
	// robustness to a few outlier frames
	feature_cache->max_feature_responses = (double*)ptr;
	ptr += num_base_features*sizeof(double);
	feature_cache->min_feature_responses = (double*)ptr;
	ptr += num_base_features*sizeof(double);
	feature_cache->ave_feature_responses = (double*)ptr;
	ptr += num_base_features*sizeof(double);
	for(i = 0; i < num_base_features; i++) {
		feature_cache->ave_feature_responses[i] = (feature_cache->integral_features[i][T] -
			feature_cache->integral_features[i][0]) / T;
		feature_cache->max_feature_responses[i] = -INFINITY;
		feature_cache->min_feature_responses[i] = INFINITY;
		for(j = 0; j < T; j++) {
			f = feature_cache->smoothed_features[i][j];
			if(f < feature_cache->min_feature_responses[i])
				feature_cache->min_feature_responses[i] = f;
			if(f > feature_cache->max_feature_responses[i])
				feature_cache->max_feature_responses[i] = f;
		}
	}

	// Compute histogram features.  We discretize the space of values for the i_th feature into num_histogram bins
	// of equal size between min_feature_responses[i] and max_feature_responses[i].  histogram_bins[i][j] stores the
	// index of the assigned bin for the i_th feature at the j_th time step
	feature_cache->histogram_bins = (int**)ptr;
	ptr += num_base_features*sizeof(int*);
	for(i = 0; i < num_base_features; i++, ptr += sizeof(int)*T) {
		feature_cache->histogram_bins[i] = (int*)ptr;
		for(j = 0; j < T; j++) {
			f = feature_cache->features[i][j];
			for(k = 0; k < feature_params[i].num_histogram_bins-1; k++) 
				if(f <= histogram_thresholds[i][k]) 
					break;
			feature_cache->histogram_bins[i][j] = k;
		}
	}

	// Compute integral features on the histogram features (e.g. using the same method as for integral_features),
	// such that the total count of the j_th histogram bin for the i_th feature can be computed using
	//   integral_histogram_features[i][j][e]-integral_histogram_features[i][j][s]
	feature_cache->integral_histogram_features = (double***)ptr;
	ptr += num_base_features*sizeof(double**);
	for(i = 0; i < num_base_features; i++) {
		feature_cache->integral_histogram_features[i] = (double**)ptr;
		ptr += sizeof(double*)*feature_params[i].num_histogram_bins;
		for(j = 0; j < feature_params[i].num_histogram_bins; j++) {
			feature_cache->integral_histogram_features[i][j] = ((double*)ptr)+T+1;
			ptr += (T+1)*3*sizeof(double);
			feature_cache->integral_histogram_features[i][j][-T-1] = 0;
			for(k = -T-1; k < 2*(T+1)-1; k++) {
				feature_cache->integral_histogram_features[i][j][k+1] = feature_cache->integral_histogram_features[i][j][k] + 
					(feature_cache->histogram_bins[i][k < 0 ? 0 : (k >= T ? T-1 : k)] == j ? 1 : 0);
			}
		}
	}

	feature_cache->bout_max_feature_responses = (double*)ptr;
	ptr += num_base_features*sizeof(double);
	feature_cache->bout_min_feature_responses = (double*)ptr;
	ptr += num_base_features*sizeof(double);

	feature_cache->memory_buffer = ptr;
}


/*
* Allocate space for feature caches used to compute bout-level features efficiently
*/
BehaviorBoutFeatures *SVMBehaviorSequence::create_behavior_bout_feature_cache(void *d, bool compute_cache) {
	int i;
	int T = num_frames(d);
	SVMFeatureParams *p;
	unsigned int cache_features_size = 
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
		p = &feature_params[i];
		cache_features_size +=
			num_base_features*sizeof(double*)*p->num_histogram_bins + //integral_histogram_features[i]
			num_base_features*sizeof(double)*p->num_histogram_bins*(T+1)*3*sizeof(double)+10000; //integral_histogram_features[i][j]
	}
	BehaviorBoutFeatures *feature_cache = (BehaviorBoutFeatures*)my_malloc(sizeof(BehaviorBoutFeatures) + cache_features_size);
	feature_cache->memory_buffer = (unsigned char*)(feature_cache+1);

	// Extract regular raw features
	feature_cache->features = (double**)feature_cache->memory_buffer;
	feature_cache->memory_buffer += num_base_features*sizeof(double*);
	for(i = 0; i < num_base_features; i++) {
		feature_cache->features[i] = ((double*)feature_cache->memory_buffer);
		feature_cache->memory_buffer += T*sizeof(double);
	}
	feature_cache->frame_times = (double*)feature_cache->memory_buffer;
	feature_cache->memory_buffer += T*sizeof(double);

	feature_cache->num_frames = T;
	feature_cache->partial_label = NULL;
	feature_cache->data = d;

	load_behavior_bout_features(d, feature_cache);
	feature_cache->fvec = NULL;
	feature_cache->bout_start = feature_cache->bout_end = -1;

	if(compute_cache)
		compute_bout_feature_caches(feature_cache);

	return feature_cache;
}



SAMPLE      SVMBehaviorSequence::read_struct_examples(char *file, STRUCT_LEARN_PARM *sparm)
{
	/* Reads struct examples and returns them in sample. The number of
	examples must be written into sample.n */
	SAMPLE   sample;  /* sample */
	EXAMPLE  *examples;
	long     n; /* number of examples */
	int num, i, j, beh;       
	char **train_list = load_examples(file, &num);
	void *d;
	BehaviorBoutSequence *bouts;
	bool computeClassTransitions = class_training_transitions ? false : true;
	BehaviorBoutFeatures *feature_cache;

	examples=(EXAMPLE *)my_malloc(sizeof(EXAMPLE)*num);

	// Keep track of the number of transitions between each pair of classes
	if(computeClassTransitions) {
		class_training_transitions = (int***)my_malloc(behaviors->num*sizeof(int**));
		class_training_transitions_count = (int**)my_malloc(behaviors->num*sizeof(int*));
		class_training_count = (int**)my_malloc(behaviors->num*sizeof(int*));
		for(beh = 0; beh < behaviors->num; beh++) {
			class_training_transitions[beh] = (int**)my_malloc(num_classes[beh]*sizeof(int*));
			class_training_transitions_count[beh] = (int*)my_malloc(num_classes[beh]*sizeof(int));
			class_training_count[beh] = (int*)my_malloc(num_classes[beh]*sizeof(int));
			for(i = 0; i < num_classes[beh]; i++) {
				class_training_transitions[beh][i] = (int*)my_malloc(num_classes[beh]*sizeof(int));
				class_training_count[beh][i] = class_training_transitions_count[beh][i] = 0;
				for(j = 0; j < num_classes[beh]; j++) {
					class_training_transitions[beh][i][j] = 0;
				}
			}
		}
	}

	n = 0;
	for(j = 0; j < num; j++) {
		d = load_training_example(train_list[j]);
		if(!d) continue;

		bouts = create_behavior_bout_sequence(d, behaviors, false);
		bouts->features = feature_cache = create_behavior_bout_feature_cache(d, false);
		strcpy(bouts->features->fname, train_list[j]);
		examples[n].x.data = feature_cache;
		examples[n].y.data = bouts;
		strcpy(examples[n].labelname, getLabelName(d));

		if(computeClassTransitions) {
			for(beh = 0; beh < behaviors->num; beh++) {
				for(i = 0; i < bouts->num_bouts[beh]; i++) {
					class_training_count[beh][bouts->bouts[beh][i].behavior]++;
					if(i)
						class_training_transitions[beh][bouts->bouts[beh][i].behavior][bouts->bouts[beh][i-1].behavior]++;
				}
			}	 
		} 
		n++;
	}
	if(computeClassTransitions) {
		if(struct_verbosity>=0)
			printf("Training examples %d behavior sequences\n", (int)n);


		// Compute all feature cache data structures for all training examples.  This requires first computing some
		// statistics of the training set features for normalization purposes
		compute_feature_mean_variance_median_statistics(examples, n);

		// Compute features for each training bout, normalized to have (0,1) mean and standard deviation
		for(j = 0; j < n; j++) {
			BehaviorBoutFeatures *behavior_bout = ((BehaviorBoutFeatures*)examples[j].x.data);
			behavior_bout->fvec = psi(examples[j].x, examples[j].y, NULL, sparm);
		}

		if(strlen(sparm->debugdir) && sparm->debug_features) {
			char fname[1000];
			CreateDirectoryIfNecessary(sparm->debugdir);
			sprintf(fname, "%s/features_unnormalized.txt", sparm->debugdir);
			print_features(fname, examples, n, false);
			sprintf(fname, "%s/features_normalized.txt", sparm->debugdir);
			print_features(fname, examples, n, true);
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
	}else {
		for(j = 0; j < n; j++) {
			BehaviorBoutFeatures *x = (BehaviorBoutFeatures*)examples[j].x.data;
			load_behavior_bout_features(x->data, x);
			compute_bout_feature_caches(x);
		}
	}


	sample.n=n;
	sample.examples=examples;
	return(sample);
}


void        SVMBehaviorSequence::init_struct_model(SAMPLE sample, STRUCTMODEL *sm, 
	STRUCT_LEARN_PARM *sparm, LEARN_PARM *lparm, 
	KERNEL_PARM *kparm)
{
	/* Initialize structmodel sm. The weight vector w does not need to be
	initialized, but you need to provide the maximum size of the
	feature space in sizePsi. This is the maximum number of different
	weights that can be learned. Later, the weight vector w will
	contain the learned weights for the model. */


	sm->sizePsi=sizePsi;
	if(struct_verbosity>=0)
		printf("Training set properties: %d features, %d behavior groups\n",
		num_features,behaviors->num);
	if(struct_verbosity>=0)
		printf("Size of Psi: %ld\n",sm->sizePsi);


	if(strlen(sparm->debugdir)) {
		char fname[1000];
		CreateDirectoryIfNecessary(sparm->debugdir);
		if(sparm->debug_predictions) {
			sprintf(fname, "%s/index.html", sparm->debugdir);
			FILE *fout = fopen(fname, "w");
			assert(fout);
			fprintf(fout, "<html></body><h2>Iteration %d</h2>\n", 0);
			fclose(fout);
			sprintf(fname, "%s/iter%d.html", sparm->debugdir, 1);
			fout = fopen(fname, "w");
			if(fout) fclose(fout);
		}
	}
}



CONSTSET    SVMBehaviorSequence::init_struct_constraints(SAMPLE sample, STRUCTMODEL *sm, 
	STRUCT_LEARN_PARM *sparm)
{
	/* Initializes the optimization problem. Typically, you do not need
	to change this function, since you want to start with an empty
	set of constraints. However, if for example you have constraints
	that certain weights need to be positive, you might put that in
	here. The constraints are represented as lhs[i]*w >= rhs[i]. lhs
	is an array of feature vectors, rhs is an array of doubles. m is
	the number of constraints. The function returns the initial
	set of constraints. */
	CONSTSET c;
	BehaviorBoutSequence *ybar;
	char fname[400];
	EXAMPLE *ex;
	SVECTOR *fy, *fybar;
	LABEL yybar;
	double lossval, factor;

	/* normal case: start with empty set of constraints */
	c.lhs=NULL;
	c.rhs=NULL;
	c.m=0;

	if(strlen(load_const_set_fname)) {
		// Load a set of constraints from disk
		c = load_const_set(save_const_set_fname);
	} else if(load_constraints_fname) {
		// We can use the labels ybar, returned from find_most_violated_constraint from some previous run 
		// of the algorithm.  This lets us quickly synthesize a bunch of constraints
		// that are likely very informative without doing any significant computation.
		// This is helpful for situations where we would like to try learning with different feature sets
		// or different parameters.  The behavior segmentations searched through when training for a different
		// set of parameters are likely still informative.
		FILE *fin = fopen(load_constraints_fname, "r");
		if(fin) {
			int last_iter = -1, iter, num_ex, last_num_ex=-1;
			double *lhs_n = create_nvector(sm->sizePsi);
			double rhs = 0;
			SVECTOR *lhs = NULL;
			int num = 0;
			clear_nvector(lhs_n,sm->sizePsi);

			while((ybar=read_bout_sequence(fin, fname, &iter, &num_ex)) != NULL) {
				if((ex=find_example(sample, fname)) != NULL) {
					if(iter != last_iter && last_iter != -1) {
						if(num == /*sample.n*/ num_ex) {
							// Add a constraint, which is the sum over each training example
							c.lhs = (DOC**)realloc(c.lhs, (c.m+1)*sizeof(DOC*));
							lhs=create_svector_n_r(lhs_n,sm->sizePsi,NULL,1.0,COMPACT_ROUNDING_THRESH);
							c.lhs[c.m] = create_example(c.m,0,1,1,lhs);
							c.rhs = (double*)realloc(c.rhs, (c.m+1)*sizeof(double));
							c.rhs[c.m++] = rhs/num;
							fprintf(stderr, "Adding cached constraint %d\n", last_iter);
						} else
							fprintf(stderr, "Skipping constraint in iteration %d, expected %d training examples found %d\n", last_iter, sample.n, num);
						clear_nvector(lhs_n,sm->sizePsi);
						rhs = 0;
						num = 0;
					} else if(last_iter != -1)
						assert(num_ex == last_num_ex);
					last_num_ex = num_ex;
					last_iter = iter;

					// Compute psi(x,y)-psi(x,ybar) and loss(y,ybar), and use it to update the constraint for this iteration
					yybar.data = ybar;
					fybar=psi(ex->x,yybar,sm,sparm);
					lossval=loss(ex->y,yybar,sparm);
					if(sparm->loss_type == SLACK_RESCALING)
						factor=lossval/sample.n;
					else                 /* do not rescale vector for */
						factor=1.0/sample.n;      /* margin rescaling loss type */
					fy=psi(ex->x,ex->y,sm,sparm);
					mult_svector_list(fy,factor);
					mult_svector_list(fybar,-factor);
					append_svector_list(fybar,fy);   /* compute fy-fybar */

					add_list_n_ns(lhs_n,fybar,1.0);
					rhs += lossval; 
					free_svector(fy);
					num++;
				} else
					fprintf(stderr, "Example %s not found\n", fname);
				free_behavior_bout_sequence(ybar, behaviors->num);
			}
			fclose(fin);
		}
	}

	return(c);
}

EXAMPLE *SVMBehaviorSequence::find_example(SAMPLE s, const char *fname) {
	int i, j;
	const char *f;
	bool found;

	for(i = 0; i < s.n; i++) {
		// String comparison, ignoring spaces
		f = ((BehaviorBoutFeatures*)s.examples[i].x.data)->fname;
		if(strlen(f) != strlen(fname))
			continue;
		found = true;
		for(j = 0; j < (int)strlen(fname); j++) {
			if(fname[j] != f[j] && (!isspace(fname[j]) && !isspace(f[j]))) {
				found = false;
				break;
			}
		}
		if(found)
			return &s.examples[i];
	}
	return NULL;
}


void SVMBehaviorSequence::on_finished_iteration(CONSTSET c, STRUCTMODEL *sm, 
	STRUCT_LEARN_PARM *sparm, int iter_num) {
		if(strlen(save_const_set_fname))
			save_const_set(c, save_const_set_fname);

		sparm->iter = iter_num;

		char fname[1000];
		if(strlen(sparm->debugdir) && sparm->debug_weights) {
			sprintf(fname, "%s/weights_%d.txt", sparm->debugdir, iter_num);
			print_weights(fname, sm->svm_model->lin_weights+1);
		}

		if(strlen(sparm->debugdir) && sparm->debug_model) {
			sprintf(fname, "%s/model_%d.txt", sparm->debugdir, iter_num);
			write_struct_model(fname, sm, sparm);
		}

		if(strlen(sparm->debugdir) && sparm->debug_predictions) {
			sprintf(fname, "%s/index.html", sparm->debugdir);
			FILE *fout = fopen(fname, "a");
			assert(fout);
			fprintf(fout, "<br><br><h2>Iteration %d</h2><a href=\"weights_%d.txt\">weights</a>|<a href=\"iter%d.html\">predictions</a>\n", iter_num, iter_num, iter_num);
			fclose(fout);
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
LABEL       SVMBehaviorSequence::inference_via_dynamic_programming(SPATTERN *x, STRUCTMODEL *sm, 
	STRUCT_LEARN_PARM *sparm, LABEL *yy)
{
	BehaviorBoutSequence *y = (BehaviorBoutSequence*)(yy ? yy->data : NULL);
	LABEL yybar;
	unsigned int sz = sizeof(BehaviorBoutSequence)+behaviors->num*(sizeof(int)+sizeof(BehaviorBout*)+2*sizeof(double));
	BehaviorBoutSequence *ybar = (BehaviorBoutSequence*)my_malloc(sz);
	BehaviorBoutFeatures *b = (BehaviorBoutFeatures*)x->data;
	int T = b->num_frames, beh, num_durations = 0, last_gt, last_partial;
	int tt, t, c, t_p, c_p, i, j, cc, restrict_c, restrict_c_p, next_duration;
	double *tmp_features = (double*)my_malloc((num_features+1)*sizeof(double));
	double **class_weights;
	double **transition_weights;
	double dur, inter;
	int *durations = (int*)my_malloc(sizeof(int)*(T+1)*4);
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
	bool is_first;


	// For speed, we may optionally choose to use a geometrically increasing series of bout durations instead
	// of searching through the entire space exhaustively
	if(time_approximation > 0) {
		dur = min_bout_duration;
		while(dur <= T*(1+time_approximation)) {
			if(!num_durations || (int)(dur) != durations[num_durations-1])
				durations[num_durations++] = (int)(dur);
			dur *= (1+time_approximation);
		}
	} else {
		for(i = 1; i <= T; i++)
			durations[num_durations++] = i;
	}


	// Initialize ybar
	memset(ybar, 0, sz);
	ybar->bouts = (BehaviorBout**)(ybar+1);
	ybar->num_bouts = (int*)(ybar->bouts+behaviors->num);
	ybar->scores = (double*)(ybar->num_bouts+behaviors->num);
	ybar->losses = (double*)(ybar->scores+behaviors->num);
	ybar->slack = ybar->score = -sm->svm_model->b;
	ybar->loss = 0;
	if(y) {
		ybar->features = y->features;
		ybar->behaviors = y->behaviors;
#if DEBUG > 0
		y->score = -sm->svm_model->b;
#endif
	}


	ptr = sm->svm_model->lin_weights+1;
	for(beh = 0; beh < behaviors->num; beh++) {
		// Currently optimizes each behavior group independently
		if(behavior >= 0 && beh != behavior) {
			ptr += num_classes[beh]*(num_classes[beh]+num_features);
			continue;
		}

#if DEBUG > 0
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

								outA = (FILE **)my_malloc(numA*sizeof(FILE*));
								writeHeader = (int*)my_malloc(numA*sizeof(int));

								for(cc=0; cc < numA; cc++) {
//									sprintf(filename, "%s\\debug_%s_%d_%d__%d.csv", /*sparm->debugdir*/"C:\\tmp\\every_line", currFileWithoutPath, DEBUG__FILESEPARATION, beh, cc);
									sprintf(filename, "%s\\debug_%s_%d_%d__%d.csv", sparm->debugdir, currFileWithoutPath, DEBUG__FILESEPARATION, beh, cc);
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


		class_weights = (double**)my_malloc(2*num_classes[beh]*sizeof(double*));
		transition_weights = class_weights+num_classes[beh];
		fn = y ? (double*)my_malloc(sizeof(double)*y->num_bouts[beh]*2) : NULL;
		dur_gt = fn ? fn + y->num_bouts[beh] : NULL;
		table = (double**)my_malloc((T+1)*(sizeof(double*)+num_classes[beh]*sizeof(double)));
		states = (BehaviorBout**)my_malloc((T+1)*(sizeof(BehaviorBout*)+num_classes[beh]*sizeof(BehaviorBout)));
		loss_score = loss_fn = loss_fp = 0;


		if(b->partial_label) {
			// When given a user supplied label, it may be the case that some class labels or label 
			// sequence in the partial labelling never appeared in the training set.  We add a couple of
			// checks here to protect against this case (otherwise our algorithm would find no solution)
			old_class_transition_counts = (int*)my_malloc(num_classes[beh]*sizeof(int)*2);
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
#ifdef ALLOW_SAME_TRANSITIONS
		double *unary_costs = (double*)my_malloc(num_classes[beh]*sizeof(double*));
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
#error invalid values of SOLVE_UNARY_COST_CONFLICT, expected values are 0 or 1.
#endif
		// else: nothing to do (don't allow same transitions, no extra unary costs read)
#endif

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
			table[0][c] = max_fn_cost;
			states[0][c].start_frame = states[0][c].end_frame = states[0][c].behavior = -1;
		}

		gt_bout[0] = 0;
		partial_label_bout[0] = 0;
		for(t = 1; t <= T; t++) { // looping through all possible bout ends // T .. length of annotation window


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
				if(y && (t_p < 0 ? -1 : gt_bout[t_p]) != last_gt && t_p != y->bouts[beh][last_gt].start_frame) {
					t_p = y->bouts[beh][last_gt].start_frame;
					next_duration = 0;
				}
				if(b->partial_label && b->partial_label->num_bouts[beh] && 
					(t_p < 0 ? -1 : partial_label_bout[t_p]) != last_partial && t_p != b->partial_label->bouts[beh][last_partial].start_frame) {
						t_p = b->partial_label->bouts[beh][last_partial].start_frame;
						next_duration = 0;
				}
				// END squeeze
				tt += next_duration;
				if(t_p <= 0) {
					t_p = 0;
					next_duration = -1;
				}
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
				psi_bout(b, t_p, t, beh, -1, tmp_features, true, !is_first); // compute bout-level features // !!! use startframe of c in case c==c_p
				is_first = false;
				for(c = 0; c < num_classes[beh]; c++) {
					if(class_training_count[beh][c] && (restrict_c_p<=0 || c == restrict_c_p)) {
						if(restrict_behavior_features[beh] && restrict_behavior_features[beh][c]) {
							printf("restrict_behavior_features currently not used\n");
							assert(0);
							bout_scores[c] = sprod_nn_map(class_weights[c], tmp_features, num_features, restrict_behavior_features[beh][c]); // untested
						}
						else
							bout_scores[c] = sprod_nn(class_weights[c], tmp_features, num_features); // function f_b(x,s,e) (8), which is summed to compute f(x,y) (7)
					}
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
#ifdef ALLOW_SAME_TRANSITIONS
							getTransitionScore(transition_score, transition_weights[c_p][c_p], unary_costs[c_p], 1, t != T ? transition_weights[c_p][c] : 0)
#else
							getTransitionScore(transition_score, transition_weights[c_p][c_p], 1, t != T ? transition_weights[c_p][c] : 0)
#endif


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
									if(c_p == y->bouts[beh][j].behavior && inter > 0) {
										loss_fn -= fn[j]*inter/dur_gt[j]; // loss_fm is a negative number (initialized 0, then subtracting), 
										// the absolute false negative cost is not used but rather the relative FN cost, 
										// where the "common shared" constant is left out, 
										// but since we are only interested in the maximum, that doesn't matter
										loss_fp -= fp*inter/dur; // ... and subtract agreeing frames from that maximum
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
							} // end if(f > table[t][c]) {  
#if DEBUG > 0
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
								sprintf(filename, "%s/debug_f_%d_%d_%d_%d.csv", sparm->debugdir, t, t_p, c, c_p);
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
		ybar->bouts[beh] = (BehaviorBout*)my_malloc(sizeof(BehaviorBout)*ybar->num_bouts[beh]);

		// Backtrack one more time to actually store that solution
		t = T; c = 0; 
		ybar->scores[beh] = 0; 
		ybar->losses[beh] = max_fn_cost;
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
		yybar.data = ybar;
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
			assert(my_abs(ybar->losses[beh] - loss2(*yy, yybar, sparm, beh, 1)) < .01);
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
				y->bouts[beh][i].bout_score = sprod_nn(class_weights[y->bouts[beh][i].behavior], tmp_features, num_features);
#ifdef ALLOW_SAME_TRANSITIONS
				getTransitionScore(y->bouts[beh][i].transition_score,                                   unary_costs[y->bouts[beh][i].behavior], i < y->num_bouts[beh]-1, transition_weights[y->bouts[beh][i].behavior][y->bouts[beh][i+1].behavior]); // unary costs added, same-transition costs added whenever same transition     is used (y->bouts[beh][i].behavior and y->bouts[beh][i+1].behavior MAY BE the SAME behavior)
#else
				getTransitionScore(y->bouts[beh][i].transition_score, transition_weights[y->bouts[beh][i].behavior][y->bouts[beh][i].behavior], i < y->num_bouts[beh]-1, transition_weights[y->bouts[beh][i].behavior][y->bouts[beh][i+1].behavior]); // unary costs stored in same transition costs;            same transitions never used (y->bouts[beh][i].behavior and y->bouts[beh][i+1].behavior are ALWAYS DIFFERENT behaviors)
#endif
//				y->bouts[beh][i].transition_score = transition_weights[y->bouts[beh][i].behavior][y->bouts[beh][i].behavior];
//				if(i < y->num_bouts[beh]-1) 
//					y->bouts[beh][i].transition_score += transition_weights[y->bouts[beh][i].behavior][y->bouts[beh][i+1].behavior];
				y->bouts[beh][i].loss_fn = y->bouts[beh][i].loss_fp = 0;
				y->scores[beh] += y->bouts[beh][i].bout_score + y->bouts[beh][i].transition_score;

#if DEBUG > 0
				if(y->scores[beh] > (i < y->num_bouts[beh]-1 ? table[y->bouts[beh][i+1].start_frame][y->bouts[beh][i+1].behavior] : table[T][0])) {
					// Something went wrong, it might be informative for debugging (break here using gdb) to test the same dynamic 
					// programming problem but without using loss, then compare y_max to y
					g_table = table; g_states = states; g_y = y;
					LABEL yy_max = classify_struct_example(*x, sm, sparm);
					BehaviorBoutSequence *y_max = (BehaviorBoutSequence*)(yy_max.data);
					assert(0);
				}
#endif
			}
			y->score += y->scores[beh];
#if DEBUG > 0
			assert(ybar->score+ybar->loss >= y->score);
#endif
		}

#if DEBUG > 0
								for(cc=0; cc < numA; cc++) {
									fclose(outA[cc]);
								}
								free(outA);
								free(writeHeader);
#endif

		free(class_weights);
		free(table);
		free(states);
		if(fn) free(fn);
	} // for(beh = 0; beh < behaviors->num; beh++) {


#if DEBUG > 0
	DOC *doc;
	SVECTOR *fvec;
	//LABEL y_max;
	double y_score;
	if(y) {
		// Compare the score computed for y when looping through the bouts in y to the score when evaluating
		// w*psi(x,y).  If this fails, something is probably wrong with the dynamic programming algorithm with
		// respect to the extraction of features or bout scores, or something is wrong with the function psi()
		y_score = (sprod_ns(sm->svm_model->lin_weights, b->fvec)*b->fvec->factor - sm->svm_model->b);
		assert(my_abs(y_score - y->score) < .01);

		// Further sanity checks that are probably redundant with earlier ones
		ybar->slack -= y->score;
		assert(ybar->slack >= -0.01);
		//y_max = classify_struct_example(*x, sm, sparm);
		//assert(((BehaviorBoutSequence*)y_max.data)->score >= y->score);
		fvec = psi(*x, yybar, sm, sparm);
		doc=create_example(1,0,1,1,fvec);
		assert(my_abs(ybar->score - (sprod_ns(sm->svm_model->lin_weights, fvec)*fvec->factor - sm->svm_model->b)) < .1);
		assert(my_abs(ybar->score-classify_example(sm->svm_model,doc)) < .1);
		free_example(doc,0);
		free_svector(fvec);

#if DEBUG > 2
		printf("w:"); for(i = 1; i <= sizePsi; i++) { printf(" %d:%f", i, (float)sm->svm_model->lin_weights[i]); } printf("\n");
		printf("psi:"); SWORD *ptr=b->fvec->words; while(ptr->wnum) { printf(" %d:%f", ptr->wnum, (float)ptr->weight); ptr++; } printf("\n");
#endif

#if DEBUG > 1
		printf("y: "); write_label(stdout, *yy); printf("\n");
		printf("ybar: "); write_label(stdout, yybar); printf("\n");
		printf("\n\n");
#endif
	} 
#endif

	free(tmp_features);
	free(durations);

	return(yybar);
}

void SVMBehaviorSequence::write_label(LABEL y, FILE *fout, int iter, int num_examples) {
	write_bout_sequence((BehaviorBoutSequence*)y.data, fout, iter, num_examples);
}

void SVMBehaviorSequence::on_finished_find_most_violated_constraint(LABEL *ybar, LABEL *y, int iter, STRUCT_LEARN_PARM *sparm, const char *ename) {
	// Save a visualization of all predicted bouts
	if(strlen(sparm->debugdir) && sparm->debug_predictions) {
		char *html=(char*)malloc(10000000), *html_gt=(char*)malloc(10000000), folder[1000], file[1000], fname[1000];
		ExtractFolderAndFileName(ename, folder, file);
		StripFileExtension(file);
		int beh = this->behavior >= 0 ? this->behavior : 0;
		if(iter >= 0) sprintf(fname, "%s/%s_%d", sparm->debugdir, file, iter);
		else sprintf(fname, "%s/%s", sparm->debugdir, file);
		VisualizeBouts((BehaviorBoutSequence*)ybar->data, behaviors, beh, fname, html);
		if(y) {
			sprintf(fname, "%s/%s_gt_%d", sparm->debugdir, file, iter);
			VisualizeBouts((BehaviorBoutSequence*)y->data, behaviors, beh, fname, html_gt);
		} 
		if(iter >= 0)
			sprintf(fname, "%s/iter%d.html", sparm->debugdir, iter);
		else
			sprintf(fname, "%s/index.html", sparm->debugdir);
		FILE *fout = fopen(fname, "a");
		assert(fout);
		BehaviorBoutSequence *yybar = (BehaviorBoutSequence*)ybar->data;
		if(y) {
			BehaviorBoutSequence *yy = (BehaviorBoutSequence*)y->data;
			fprintf(fout, "<br><br>%s: %s, score=%f, loss=%f<br>%s<br>%s: ground truth, score=%f, loss=%f<br>%s\n", file, iter >= 0 ? "most violated" : "best_score", (float)yybar->score, (float)yybar->loss, html, file, (float)yy->score, (float)yy->loss, html_gt);
		} else 
			fprintf(fout, "<br><br>%s: best score, score=%f, loss=%f<br>%s\n", file, (float)yybar->score, (float)yybar->loss, html, file);
		fclose(fout);
		free(html);
		free(html_gt);
	}
}

LABEL SVMBehaviorSequence::read_label(FILE *fin, char *fname, int *iter, int *num_examples) {
	LABEL y;
	y.data = read_bout_sequence(fin, fname, iter, num_examples);
	return y;
}

void SVMBehaviorSequence::write_bout_sequence(BehaviorBoutSequence *b, FILE *fout, int iter, int num_examples) {
	int beh, i;
	char fname[10000];

	strcpy(fname, b->features->fname);
	for(i = 0; i < (int)strlen(fname); i++)
		if(isspace(fname[i]))
			fname[i] = '_';
	fprintf(fout, "%s %d %d ", fname, iter, num_examples);
	for(beh = 0; beh < behaviors->num; beh++) {
		fprintf(fout, "%d ", b->num_bouts[beh]);
		for(i = 0; i < b->num_bouts[beh]; i++) 
			fprintf(fout, ", %d %d %d", b->bouts[beh][i].behavior, b->bouts[beh][i].start_frame, b->bouts[beh][i].end_frame);
		fprintf(fout, ";");
	}
	fprintf(fout, "\n");
} 

BehaviorBoutSequence *SVMBehaviorSequence::read_bout_sequence(FILE *fin, char *fname, int *iter, int *num_examples) {
	int i, beh;
	unsigned int sz = sizeof(BehaviorBoutSequence)+behaviors->num*(sizeof(int)+sizeof(BehaviorBout*)+2*sizeof(double));
	BehaviorBoutSequence *y = (BehaviorBoutSequence*)my_malloc(sz);
	char fn[400];
	char *f = fname ? fname : fn;

	memset(y, 0, sz);
	y->bouts = (BehaviorBout**)(y+1);
	y->num_bouts = (int*)(y->bouts+behaviors->num);

	if(fscanf(fin, "%s %d %d ", f, iter, num_examples) < 3) { free(y); return NULL; }
	for(beh = 0; beh < behaviors->num; beh++) {
		if(fscanf(fin, "%d ", &y->num_bouts[beh]) <= 0) { free(y); return NULL; }
		y->bouts[beh] = (BehaviorBout*)my_malloc(sizeof(BehaviorBout)*y->num_bouts[beh]);
		for(i = 0; i < y->num_bouts[beh]; i++) 
			assert(fscanf(fin, ", %d %d %d", &y->bouts[beh][i].behavior, &y->bouts[beh][i].start_frame, &y->bouts[beh][i].end_frame) == 3);
		fscanf(fin, ";");
	}
	fscanf(fin,"%*[^\n]\n");

	return y;
}



LABEL       SVMBehaviorSequence::classify_struct_example(SPATTERN x, STRUCTMODEL *sm, 
	STRUCT_LEARN_PARM *sparm)
{
	/* Finds the label yhat for pattern x that scores the highest
	according to the linear evaluation function in sm, especially the
	weights sm.w. The returned label is taken as the prediction of sm
	for the pattern x. The weights correspond to the features defined
	by psi() and range from index 1 to index sm->sizePsi. If the
	function cannot find a label, it shall return an empty label as
	recognized by the function empty_label(y). */
	return inference_via_dynamic_programming(&x, sm, sparm, NULL);
}


LABEL       SVMBehaviorSequence::find_most_violated_constraint_slackrescaling(SPATTERN x, LABEL y, 
	STRUCTMODEL *sm, 
	STRUCT_LEARN_PARM *sparm) {
		/* Finds the label ybar for pattern x that that is responsible for
		the most violated constraint for the slack rescaling
		formulation. For linear slack variables, this is that label ybar
		that maximizes

		argmax_{ybar} loss(y,ybar)*(1-w*psi(x,y)+w*psi(x,ybar)) 

		Note that ybar may be equal to y (i.e. the max is 0), which is
		different from the algorithms described in
		[Tschantaridis/05]. Note that this argmax has to take into
		account the scoring function in sm, especially the weights sm.w,
		as well as the loss function, and whether linear or quadratic
		slacks are used. The weights in sm.w correspond to the features
		defined by psi() and range from index 1 to index
		sm->sizePsi. Most simple is the case of the zero/one loss
		function. For the zero/one loss, this function should return the
		highest scoring label ybar (which may be equal to the correct
		label y), or the second highest scoring label ybar, if
		Psi(x,ybar)>Psi(x,y)-1. If the function cannot find a label, it
		shall return an empty label as recognized by the function
		empty_label(y). */

		return inference_via_dynamic_programming(&x, sm, sparm, &y);
}


LABEL       SVMBehaviorSequence::find_most_violated_constraint_marginrescaling(SPATTERN x, LABEL y, 
	STRUCTMODEL *sm, 
	STRUCT_LEARN_PARM *sparm)
{
	/* Finds the label ybar for pattern x that that is responsible for
	the most violated constraint for the margin rescaling
	formulation. For linear slack variables, this is that label ybar
	that maximizes

	argmax_{ybar} loss(y,ybar)+w*psi(x,ybar)-w*psi(x,y)

	Note that ybar may be equal to y (i.e. the max is 0), which is
	different from the algorithms described in
	[Tschantaridis/05]. Note that this argmax has to take into
	account the scoring function in sm, especially the weights sm.w,
	as well as the loss function, and whether linear or quadratic
	slacks are used. The weights in sm.w correspond to the features
	defined by psi() and range from index 1 to index
	sm->sizePsi. Most simple is the case of the zero/one loss
	function. For the zero/one loss, this function should return the
	highest scoring label ybar (which may be equal to the correct
	label y), or the second highest scoring label ybar, if
	Psi(x,ybar)>Psi(x,y)-1. If the function cannot find a label, it
	shall return an empty label as recognized by the function
	empty_label(y). */
	return inference_via_dynamic_programming(&x, sm, sparm, &y);
}

int         SVMBehaviorSequence::empty_label(LABEL y)
{
	/* Returns true, if y is an empty label. An empty label might be
	returned by find_most_violated_constraint_???(x, y, sm) if there
	is no incorrect label that can be found for x, or if it is unable
	to label x at all */
	return ((BehaviorBoutSequence*)y.data)->slack <= 0;
}





double      SVMBehaviorSequence::loss(LABEL y, LABEL ybar, STRUCT_LEARN_PARM *sparm) {
	int i;
	double l = 0;

	if(behavior >= 0) 
		return loss2(y, ybar, sparm, behavior, 0);
	else {
		for(i = 0; i < behaviors->num; i++)
			l += loss2(y, ybar, sparm, i, 0);
	}
	return l;
}

double      SVMBehaviorSequence::loss2(LABEL yy, LABEL yybar, STRUCT_LEARN_PARM *sparm, int beh, int debug)
{
	/* loss for correct label y and predicted label ybar. The loss for
	y==ybar has to be zero. sparm->loss_function is set with the -l option. */

	BehaviorBoutSequence *y = (BehaviorBoutSequence*)yy.data;
	BehaviorBoutSequence *ybar = (BehaviorBoutSequence*)yybar.data;
	int curr_y=0, curr_ybar=0;  // current bout in y and ybar
	double l = 0;               // total loss so far
	double dur_y=0, dur_ybar=0; // duration in frames of the current bout
	double inter;               // intersection time between bouts
	double sum_y = 0, sum_ybar = 0;
	BehaviorBoutFeatures *b = y->features;
	int T = b->num_frames;
	double cl, l_fn = 0, l_fn2 = 0;

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
				assert(my_abs(l_fn - ybar->bouts[beh][curr_ybar].loss_fn) < .00001);
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

int         SVMBehaviorSequence::finalize_iteration(double ceps, int cached_constraint,
	SAMPLE sample, STRUCTMODEL *sm,
	CONSTSET cset, double *alpha, 
	STRUCT_LEARN_PARM *sparm)
{
	/* This function is called just before the end of each cutting plane iteration. ceps is the amount by which the most violated constraint found in the current iteration was violated. cached_constraint is true if the added constraint was constructed from the cache. If the return value is FALSE, then the algorithm is allowed to terminate. If it is TRUE, the algorithm will keep iterating even if the desired precision sparm->epsilon is already reached. */
	return(0);
}

void        SVMBehaviorSequence::print_struct_learning_stats(SAMPLE sample, STRUCTMODEL *sm,
	CONSTSET cset, double *alpha, 
	STRUCT_LEARN_PARM *sparm)
{
	/* This function is called after training and allows final touches to
	the model sm. But primarly it allows computing and printing any
	kind of statistic (e.g. training error) you might want. */

	/* Replace SV with single weight vector */
	MODEL *model=sm->svm_model;
	if(model->kernel_parm.kernel_type == K_LINEAR) {
		if(struct_verbosity>=1) {
			printf("Compacting linear model..."); fflush(stdout);
		}
		sm->svm_model=compact_linear_model(model);
		sm->w=sm->svm_model->lin_weights; /* short cut to weight vector */
		free_model(model,1);
		if(struct_verbosity>=1) {
			printf("done\n"); fflush(stdout);
		}
	} 
}

void        SVMBehaviorSequence::print_struct_testing_stats(SAMPLE sample, STRUCTMODEL *sm,
	STRUCT_LEARN_PARM *sparm, 
	STRUCT_TEST_STATS *teststats)
{
	/* This function is called after making all test predictions in
	svm_struct_classify and allows computing and printing any kind of
	evaluation (e.g. precision/recall) you might want. You can use
	the function eval_prediction to accumulate the necessary
	statistics for each prediction. */
}

void        SVMBehaviorSequence::eval_prediction(long exnum, EXAMPLE ex, LABEL ypred, 
	STRUCTMODEL *sm, STRUCT_LEARN_PARM *sparm, 
	STRUCT_TEST_STATS *teststats)
{
	/* This function allows you to accumlate statistic for how well the
	predicition matches the labeled example. It is called from
	svm_struct_classify. See also the function
	print_struct_testing_stats. */
	if(exnum == 0) { /* this is the first time the function is
					 called. So initialize the teststats */
	}
}

void        SVMBehaviorSequence::write_struct_model(const char *file, STRUCTMODEL *sm, 
	STRUCT_LEARN_PARM *sparm)
{
	/* Writes structural model sm to file file. */
	/* Writes structural model sm to file file. */
	FILE *modelfl;
	long j,i,sv_num;
	MODEL *model=sm->svm_model;
	SVECTOR *v;
	int beh;
	SVMFeatureParams *p;

	if ((modelfl = fopen (file, "w")) == NULL)
	{ perror (file); exit (1); }
	fprintf(modelfl,"SVM-Behavior Version %s\n",INST_VERSION);
	for(i = 0; i < behaviors->num; i++)
		fprintf(modelfl,"%d ", this->num_classes[i]);
	fprintf(modelfl,"# number of classes\n");
	fprintf(modelfl,"%d # number of features\n",
		this->num_features);
	fprintf(modelfl,"%d # number of base features\n",
		this->num_base_features);
	fprintf(modelfl,"%d # behavior group\n", behavior);
	fprintf(modelfl,"%d # feature diff frames\n", feature_diff_frames);
	fprintf(modelfl,"%d # loss function\n",
		sparm->loss_function);
	fprintf(modelfl,"%ld # kernel type\n",
		model->kernel_parm.kernel_type);
	fprintf(modelfl,"%ld # kernel parameter -d \n",
		model->kernel_parm.poly_degree);
	fprintf(modelfl,"%.8g # kernel parameter -g \n",
		model->kernel_parm.rbf_gamma);
	fprintf(modelfl,"%.8g # kernel parameter -s \n",
		model->kernel_parm.coef_lin);
	fprintf(modelfl,"%.8g # kernel parameter -r \n",
		model->kernel_parm.coef_const);
	fprintf(modelfl,"%s# kernel parameter -u \n",model->kernel_parm.custom);
	fprintf(modelfl,"%.8g # Regularization C -r \n", sparm->C);
	fprintf(modelfl,"%ld # highest feature index \n",model->totwords);
	fprintf(modelfl,"%ld # number of training documents \n",model->totdoc);


	for(beh = 0; beh < behaviors->num; beh++) {
		fprintf(modelfl, beh ? "; behavior_group %d: " : "behavior_group %d: ", beh);
		for(i = 0; i < num_classes[beh]; i++) {
			fprintf(modelfl, i ? ", value %d: raw_count=%d transition_count=%d" : "value %d: raw_count=%d transition_count=%d", 
				(int)i, class_training_count[beh][i], class_training_transitions_count[beh][i]);
			for(j = 0; j < class_training_transitions_count[beh][i]; j++)
				fprintf(modelfl, " %d", class_training_transitions[beh][i][j]);
		}
	}
	fprintf(modelfl, " # allowable class transitions\n");
	for(i = 0; i < num_base_features; i++) {
		p = &feature_params[i];
		fprintf(modelfl, FORMAT__BOUT_FEATURE_PARAMS, 
			p->feature_sample_smoothness_window, p->num_temporal_levels, p->num_bout_max_thresholds, 
			p->num_bout_min_thresholds, p->num_bout_change_points, p->num_histogram_bins, p->num_histogram_temporal_levels, 
			p->num_difference_temporal_levels, p->num_harmonic_features, p->use_bout_sum_features?1:0, p->use_bout_ave_features?1:0,
			p->use_bout_sum_absolute_features?1:0, p->use_bout_ave_absolute_features?1:0, 
			p->use_standard_deviation?1:0, p->use_sum_variance?1:0, p->use_bout_max_feature?1:0, p->use_bout_min_feature?1:0,
			p->use_global_difference_max_ave_features?1:0, p->use_global_difference_min_ave_features?1:0, p->use_global_difference_ave_ave_features?1:0,
			p->use_global_difference_max_sum_features?1:0, p->use_global_difference_min_sum_features?1:0, p->use_global_difference_ave_sum_features?1:0,
			p->use_bout_change?1:0, p->use_bout_absolute_change?1:0, p->use_histogram_sum_features?1:0,
			p->use_histogram_ave_features?1:0, p->use_sum_harmonic_features?1:0, p->use_ave_harmonic_features?1:0, p->use_sum_absolute_harmonic_features?1:0,
			p->use_ave_absolute_harmonic_features?1:0, p->use_start_sum_absolute_diff_haar_features?1:0,
			p->use_end_sum_absolute_diff_haar_features?1:0, p->use_start_sum_diff_haar_features?1:0, p->use_end_sum_diff_haar_features?1:0, 
			p->use_start_ave_absolute_diff_haar_features?1:0, p->use_end_ave_absolute_diff_haar_features?1:0, p->use_start_ave_diff_haar_features?1:0, 
			p->use_end_ave_diff_haar_features?1:0, (int)i);
	}
	for(i = 0; i < num_features; i++)
		fprintf(modelfl, "%lf ", features_mu[i]);
	fprintf(modelfl, " # mean feature responses\n");
	for(i = 0; i < num_features; i++)
		fprintf(modelfl, "%lf ", features_gamma[i]);
	fprintf(modelfl, " # feature normalization factors\n");
	for(i = 0; i < num_base_features; i++)
		for(j = 0; j < feature_params[i].num_histogram_bins; j++)
			fprintf(modelfl, "%lf ", histogram_thresholds[i][j]);
	fprintf(modelfl, " # histogram thresholds\n");
	for(i = 0; i < num_base_features; i++)
		for(j = 0; j < feature_params[i].num_bout_max_thresholds; j++)
			fprintf(modelfl, "%lf ", max_thresholds[i][j]);
	fprintf(modelfl, " # bout max thresholds\n");
	for(i = 0; i < num_base_features; i++)
		for(j = 0; j < feature_params[i].num_bout_min_thresholds; j++)
			fprintf(modelfl, "%lf ", min_thresholds[i][j]);
	fprintf(modelfl, " # bout min thresholds\n");



	sv_num=1;
	for(i=1;i<model->sv_num;i++) {
		for(v=model->supvec[i]->fvec;v;v=v->next) 
			sv_num++;
	}
	fprintf(modelfl,"%ld # number of support vectors plus 1 \n",sv_num);
	fprintf(modelfl,"%.8g # threshold b, each following line is a SV (starting with alpha*y)\n",model->b);

	for(i=1;i<model->sv_num;i++) {
		for(v=model->supvec[i]->fvec;v;v=v->next) {
			fprintf(modelfl,"%.32g ",model->alpha[i]*v->factor);
			fprintf(modelfl,"qid:%ld ",v->kernel_id);
			for (j=0; (v->words[j]).wnum; j++) {
				fprintf(modelfl,"%ld:%.8g ",
					(long)(v->words[j]).wnum,
					(double)(v->words[j]).weight);
			}
			if(v->userdefined)
				fprintf(modelfl,"#%s\n",v->userdefined);
			else
				fprintf(modelfl,"#\n");
			/* NOTE: this could be made more efficient by summing the
			alpha's of identical vectors before writing them to the
			file. */
		}
	}
	fclose(modelfl);
}

bool SVMBehaviorSequence::ReadFeatureParam(FILE *modelfl, SVMFeatureParams *p) {
	int num;
	int b[30];
	if((num=fscanf(modelfl, FORMAT__BOUT_FEATURE_PARAMS, &p->feature_sample_smoothness_window, &p->num_temporal_levels, &p->num_bout_max_thresholds, 
		&p->num_bout_min_thresholds, &p->num_bout_change_points, &p->num_histogram_bins, &p->num_histogram_temporal_levels, 
		&p->num_difference_temporal_levels, &p->num_harmonic_features, &b[0], &b[1], &b[28], &b[29], &b[26], &b[27], &b[2], &b[3], &b[4], &b[5], &b[6], &b[7], &b[8], &b[9], 
		&b[10], &b[11], &b[12], &b[13], &b[14], &b[15], &b[16], &b[17], &b[18], &b[19], &b[20], &b[21], &b[22], &b[23], &b[24], &b[25]))!=39) {
			fprintf(stderr, "ERROR parsing feature params: only read first %d numbers\n", num);
			return false;
	}
	p->use_bout_sum_features=b[0]; p->use_bout_ave_features=b[1]; p->use_bout_max_feature=b[2]; p->use_bout_min_feature=b[3];
	p->use_global_difference_max_ave_features=b[4]; p->use_global_difference_min_ave_features=b[5]; p->use_global_difference_ave_ave_features=b[6];
	p->use_global_difference_max_sum_features=b[7]; p->use_global_difference_min_sum_features=b[8]; p->use_global_difference_ave_sum_features=b[9];
	p->use_bout_change=b[10]; p->use_bout_absolute_change=b[11]; p->use_histogram_sum_features=b[12];
	p->use_histogram_ave_features=b[13]; p->use_sum_harmonic_features=b[14]; p->use_ave_harmonic_features=b[15]; p->use_sum_absolute_harmonic_features=b[16];
	p->use_ave_absolute_harmonic_features=b[17]; p->use_start_sum_absolute_diff_haar_features=b[18];
	p->use_end_sum_absolute_diff_haar_features=b[19]; p->use_start_sum_diff_haar_features=b[20]; p->use_end_sum_diff_haar_features=b[21];
	p->use_start_ave_absolute_diff_haar_features=b[22]; p->use_end_ave_absolute_diff_haar_features=b[23]; p->use_start_ave_diff_haar_features=b[24]; 
	p->use_end_ave_diff_haar_features=b[25];  p->use_standard_deviation=b[26]; p->use_sum_variance=b[27];

	p->use_bout_sum_absolute_features=b[28]; p->use_bout_ave_absolute_features=b[29]; 

	return true;
}


STRUCTMODEL SVMBehaviorSequence::read_struct_model(const char *file, STRUCT_LEARN_PARM *sparm) {
	/* Reads structural model sm from file file. This function is used
	only in the prediction module, not in the learning module. */
	FILE *modelfl;
	STRUCTMODEL sm;
	int beh_tmp, i_tmp, beh, j;
	long i=0,queryid,slackid;
	double costfactor;
	long max_sv,max_words,ll,wpos;
	char *line,*comment;
	SWORD *words;
	char version_buffer[100], *ptr, line2[10000];
	MODEL *model;
	SVMFeatureParams *p;
	int num_histogram_bins = 0, num_max_thresholds = 0, num_min_thresholds = 0;
	double tmp;

	nol_ll(file,&max_sv,&max_words,&ll); /* scan size of model file */
	max_words+=2;
	ll+=2;

	words = (SWORD *)my_malloc(sizeof(SWORD)*(max_words+10));
	line = (char *)my_malloc(sizeof(char)*ll);
	model = (MODEL *)my_malloc(sizeof(MODEL));

	if ((modelfl = fopen (file, "r")) == NULL)
	{ perror (file); exit (1); }

	fscanf(modelfl,"SVM-Behavior Version %s\n",version_buffer);
	if(strcmp(version_buffer, INST_VERSION)) {
		perror ("Version of model-file does not match version of svm_struct_classify!"); 
		exit (1); 
	}
	fgets(line2, 9999, modelfl);
	while((ptr=strtok(i ? NULL : line2, " ")) != NULL) 
		sscanf(ptr, "%d", &this->num_classes[i++]); 
	fscanf(modelfl,"%d%*[^\n]\n", &this->num_features);  
	fscanf(modelfl,"%d%*[^\n]\n", &this->num_base_features);  
	fscanf(modelfl,"%d%*[^\n]\n", &this->behavior);
	fscanf(modelfl,"%d%*[^\n]\n", &this->feature_diff_frames);
	fscanf(modelfl,"%d%*[^\n]\n", &sparm->loss_function);  
	fscanf(modelfl,"%ld%*[^\n]\n", &model->kernel_parm.kernel_type);  
	fscanf(modelfl,"%ld%*[^\n]\n", &model->kernel_parm.poly_degree);
	fscanf(modelfl,"%lf%*[^\n]\n", &model->kernel_parm.rbf_gamma);
	fscanf(modelfl,"%lf%*[^\n]\n", &model->kernel_parm.coef_lin);
	fscanf(modelfl,"%lf%*[^\n]\n", &model->kernel_parm.coef_const);
	fscanf(modelfl,"%[^#]%*[^\n]\n", model->kernel_parm.custom);
	fscanf(modelfl,"%lf%*[^\n]\n", &tmp);

	fscanf(modelfl,"%ld%*[^\n]\n", &model->totwords);
	fscanf(modelfl,"%ld%*[^\n]\n", &model->totdoc);


	class_training_transitions = (int***)my_malloc(behaviors->num*sizeof(int**));
	class_training_transitions_count = (int**)my_malloc(behaviors->num*sizeof(int*));
	class_training_count = (int**)my_malloc(behaviors->num*sizeof(int*));
	for(beh = 0; beh < behaviors->num; beh++) {
		assert(fscanf(modelfl, beh ? "; behavior_group %d: " : "behavior_group %d: ", &beh_tmp) == 1); 
		assert(beh_tmp == beh);

		class_training_transitions[beh] = (int**)my_malloc(num_classes[beh]*sizeof(int*));
		class_training_transitions_count[beh] = (int*)my_malloc(num_classes[beh]*sizeof(int));
		class_training_count[beh] = (int*)my_malloc(num_classes[beh]*sizeof(int));
		for(i = 0; i < num_classes[beh]; i++) {      
			assert(fscanf(modelfl, i ? ", value %d: raw_count=%d transition_count=%d" : "value %d: raw_count=%d transition_count=%d", 
				&i_tmp, &class_training_count[beh][i], &class_training_transitions_count[beh][i]) == 3); 
			assert(i_tmp == i);

			class_training_transitions[beh][i] = (int*)my_malloc(num_classes[beh]*sizeof(int));
			for(j = 0; j < class_training_transitions_count[beh][i]; j++)
				assert(fscanf(modelfl, " %d", &class_training_transitions[beh][i][j]) == 1);
		}
	}
	fscanf(modelfl,"%*[^\n]\n");


	for(i = 0; i < num_base_features; i++) {
		p = &feature_params[i];
		assert(ReadFeatureParam(modelfl, p));
		num_histogram_bins += p->num_histogram_bins;
		num_max_thresholds += p->num_bout_max_thresholds;
		num_min_thresholds += p->num_bout_min_thresholds;
	}


	features_mu = (double*)malloc(num_features*sizeof(double)*2);
	features_gamma = features_mu + num_features;
	histogram_thresholds = (double**)malloc(3*sizeof(double*)*num_base_features + 
		(num_histogram_bins+num_max_thresholds+num_min_thresholds)*sizeof(double));
	min_thresholds = histogram_thresholds + num_base_features;
	max_thresholds = min_thresholds + num_base_features;
	double *ptr2 = (double*)(max_thresholds + num_base_features);
	for(i = 0; i < num_features; i++)
		assert(fscanf(modelfl, "%lf", &features_mu[i]) == 1);
	fscanf(modelfl,"%*[^\n]\n");
	for(i = 0; i < num_features; i++)
		assert(fscanf(modelfl, "%lf", &features_gamma[i]));
	fscanf(modelfl,"%*[^\n]\n");
	for(i = 0; i < num_base_features; i++) {
		histogram_thresholds[i] = ptr2;  ptr2 += feature_params[i].num_histogram_bins;
		for(j = 0; j < feature_params[i].num_histogram_bins; j++)
			fscanf(modelfl, "%lf ", &histogram_thresholds[i][j]);
	}
	fscanf(modelfl,"%*[^\n]\n");
	for(i = 0; i < num_base_features; i++) {
		max_thresholds[i] = ptr2;  ptr2 += feature_params[i].num_histogram_bins;
		for(j = 0; j < feature_params[i].num_bout_max_thresholds; j++)
			fscanf(modelfl, "%lf ", &max_thresholds[i][j]);
	}
	fscanf(modelfl,"%*[^\n]\n");
	for(i = 0; i < num_base_features; i++) {
		min_thresholds[i] = ptr2;  ptr2 += feature_params[i].num_histogram_bins;
		for(j = 0; j < feature_params[i].num_bout_min_thresholds; j++)
			fscanf(modelfl, "%lf ", &min_thresholds[i][j]);
	}
	fscanf(modelfl,"%*[^\n]\n");

	fscanf(modelfl,"%ld%*[^\n]\n", &model->sv_num);
	fscanf(modelfl,"%lf%*[^\n]\n", &model->b);


	model->supvec = (DOC **)my_malloc(sizeof(DOC *)*model->sv_num);
	model->alpha = (double *)my_malloc(sizeof(double)*model->sv_num);
	model->index=NULL;
	model->lin_weights=NULL;

	for(i=1;i<model->sv_num;i++) {
		fgets(line,(int)ll,modelfl);
		if(!parse_document(line,words,&(model->alpha[i]),&queryid,&slackid,
			&costfactor,&wpos,max_words,&comment)) {
				printf("\nParsing error while reading model file in SV %ld!\n%s",
					i,line);
				exit(1);
		}
		model->supvec[i] = create_example(-1,0,0,0.0,
			create_svector(words,comment,1.0));
		model->supvec[i]->fvec->kernel_id=queryid;
	}
	fclose(modelfl);
	free(line);
	free(words);

	feature_names = (char**)malloc(sizeof(char*)*num_features);
	memset(feature_names, 0, sizeof(char*)*num_features);
	compute_feature_space_size();

	if(verbosity>=1) {
		fprintf(stdout, " (%d support vectors read) ",(int)(model->sv_num-1));
	}
	sm.svm_model=compact_linear_model(model);
	sm.w=sm.svm_model->lin_weights; /* short cut to weight vector */
	free_model(model,1);
	//sm.svm_model=model;

	sizePsi = sm.sizePsi= sm.svm_model->totwords;
	sm.w=NULL;

	beh = 0; // HACK!! Should iterate through all behaviors here
//	sm.compactUnaryCosts = (sizePsi-1) == (num_features + num_classes[beh]) * num_classes[beh]; /* CSC: true iff unary costs are stored as same transition costs (and no further weights were given) */
//	sm.extraUnaryCosts = (sizePsi-1) == (num_features + num_classes[beh] + 1) * num_classes[beh]; /* CSC: true iff case Unary costs are given in addition to same transition costs; equals !compactUnaryCosts but explicitely given for sanity checks (CHECK THE ASSIGNMENT WHENEVER ADDITIONAL WEIGHTS (other than feature & transition weights) ARE ADDED) */
	sm.compactUnaryCosts = (sizePsi-1) == compactUnarySize(num_features, num_classes[beh]); /* CSC: true iff unary costs are stored as same transition costs (and no further weights were given) */
	sm.extraUnaryCosts = (sizePsi-1) == extraUnarySize(num_features, num_classes[beh]); /* CSC: true iff case Unary costs are given in addition to same transition costs; equals !compactUnaryCosts but explicitely given for sanity checks (CHECK THE ASSIGNMENT WHENEVER ADDITIONAL WEIGHTS (other than feature & transition weights) ARE ADDED) */
	
	assert(sm.compactUnaryCosts != sm.extraUnaryCosts);

	return(sm);
}

void        SVMBehaviorSequence::write_label(FILE *fp, LABEL yy)
{
	BehaviorBoutSequence *bouts = (BehaviorBoutSequence*)yy.data;
	int i, beh;

	fprintf(fp,"score=%lf loss=%lf", bouts->score, bouts->loss);
	for(beh = 0; beh < behaviors->num; beh++) {
		fprintf(fp,"num_bouts[%d]=%d score[%d]=%lf loss[%d]=%lf", beh, bouts->num_bouts[beh], beh, bouts->scores[beh], beh, bouts->losses[beh]);
		for(i = 0; i < bouts->num_bouts[beh]; i++) 
			fprintf(fp," (%d:%d,%d,%d)", beh, bouts->bouts[beh][i].behavior, bouts->bouts[beh][i].start_frame, bouts->bouts[beh][i].end_frame);
	}
} 

void free_behavior_bout_sequence(BehaviorBoutSequence *b, int num) {
	int i;
	for(i = 0; i < num; i++)
		if(b->bouts[i])
			free(b->bouts[i]);
	free(b);
}

void        SVMBehaviorSequence::free_pattern(SPATTERN x) {
	/* Frees the memory of x. */
	BehaviorBoutFeatures *b = (BehaviorBoutFeatures*)x.data;
	if(b->data) free_data(b->data);
	if(b->partial_label) free_behavior_bout_sequence(b->partial_label, behaviors->num);
	free_svector(b->fvec);
	free(b);
}



void        SVMBehaviorSequence::free_label(LABEL y) {
	/* Frees the memory of y. */
	free_behavior_bout_sequence((BehaviorBoutSequence*)y.data, behaviors->num);
}

void        SVMBehaviorSequence::free_struct_model(STRUCTMODEL sm) 
{
	/* Frees the memory of model. */
	/* if(sm.w) free(sm.w); */ /* this is free'd in free_model */
	if(sm.svm_model) free_model(sm.svm_model,1);
	/* add free calls for user defined data here */
}

void        SVMBehaviorSequence::free_struct_sample(SAMPLE s)
{
	/* Frees the memory of sample s. */
	int i;
	for(i=0;i<s.n;i++) { 
		free_pattern(s.examples[i].x);
		free_label(s.examples[i].y);
	}
	free(s.examples);
}

void        SVMBehaviorSequence::print_struct_help()
{
	/* Prints a help text that is appended to the common help text of
	svm_struct_learn. */
	printf("         --* string  -> custom parameters that can be adapted for struct\n");
	printf("                        learning. The * can be replaced by any character\n");
	printf("                        and there can be multiple options starting with --.\n");
}

void         SVMBehaviorSequence::parse_struct_parameters(STRUCT_LEARN_PARM *sparm)
{
	/* Parses the command line parameters that start with -- */
	int i;

	for(i=0;(i<sparm->custom_argc) && ((sparm->custom_argv[i])[0] == '-');i++) {
		switch ((sparm->custom_argv[i])[2]) 
		{ 
		case 'a': i++; /* strcpy(learn_parm->alphafile,argv[i]); */ break;
		case 'e': i++; /* sparm->epsilon=atof(sparm->custom_argv[i]); */ break;
		case 'k': i++; /* sparm->newconstretrain=atol(sparm->custom_argv[i]); */ break;
		default: printf("\nUnrecognized option %s!\n\n",sparm->custom_argv[i]);
			exit(0);
		}
	}
}

void        SVMBehaviorSequence::print_struct_help_classify()
{
	/* Prints a help text that is appended to the common help text of
	svm_struct_classify. */
	printf("         --* string -> custom parameters that can be adapted for struct\n");
	printf("                       learning. The * can be replaced by any character\n");
	printf("                       and there can be multiple options starting with --.\n");
}

void         SVMBehaviorSequence::parse_struct_parameters_classify(STRUCT_LEARN_PARM *sparm)
{
	/* Parses the command line parameters that start with -- for the
	classification module */
	int i;

	for(i=0;(i<sparm->custom_argc) && ((sparm->custom_argv[i])[0] == '-');i++) {
		switch ((sparm->custom_argv[i])[2]) 
		{ 
			/* case 'x': i++; strcpy(xvalue,sparm->custom_argv[i]); break; */
		default: printf("\nUnrecognized option %s!\n\n",sparm->custom_argv[i]);
			exit(0);
		}
	}
}

#define LABEL_BOUTS 0
IplImage *VisualizeBouts(BehaviorBoutSequence *seq, BehaviorGroups *groups, int beh, const char *fname, char *html) { 
	BehaviorGroup *group = &groups->behaviors[beh];
	if(!seq->num_bouts[beh]) 
		return NULL;

	int h = LABEL_BOUTS ? 50 : 10;
	int T = seq->bouts[beh][seq->num_bouts[beh]-1].end_frame;
	CvFont font;
	CvSize sz;
	cvInitFont(&font, CV_FONT_VECTOR0, 0.5f, 0.4f, 0, 2);

	IplImage *img = cvCreateImage(cvSize(T, h), IPL_DEPTH_8U, 3);
	cvZero(img);

	// Draw bouts as colored rectangles
	int i;
	for(i = 0; i < seq->num_bouts[beh]; i++) {
		int c = seq->bouts[beh][i].behavior;
		cvRectangle(img, cvPoint(seq->bouts[beh][i].start_frame,0), cvPoint(seq->bouts[beh][i].end_frame,h), 
			CV_RGB((group->values[c].color & 0xff0000)>>16,
			(group->values[c].color & 0x00ff00)>>8,
			(group->values[c].color & 0xff)), CV_FILLED);
	}

#if LABEL_BOUTS
	// Draw labels for bouts
	int prev_prev_max_x = -100000, prev_max_x = -100000;
	int last_pos = 2, last_last_pos = 1;
	int y[3], ymin;
	for(i = 0; i < seq->num_bouts[beh]; i++) {
		int c = seq->bouts[beh][i].behavior;
		int m = (seq->bouts[beh][i].start_frame+seq->bouts[beh][i].end_frame)/2;

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
		sprintf(fname2, "%s.png", fname);
		cvSaveImage(fname2, img);
		cvReleaseImage(&img);

		if(html) {
			ExtractFolderAndFileName(fname, folder, fname2);    
			sprintf(html_tmp, "<img src=\"%s.png\" usemap=\"#%s\" height=50 />\n<map name=\"%s\">\n", fname2, fname2, fname2);

			char str[10000], alt[1000];
			char *ptr = html_tmp+strlen(html_tmp);
			for(i = 0; i < seq->num_bouts[beh]; i++) {
				int c = seq->bouts[beh][i].behavior;
				float z = 50.0f/h;
				sprintf(alt, "behavior=%s bout_score=%f transition_score=%f loss_fn=%f loss_fp=%f", group->values[c].name, 
					(float)seq->bouts[beh][i].bout_score, (float)seq->bouts[beh][i].transition_score, (float)seq->bouts[beh][i].loss_fn,  (float)seq->bouts[beh][i].loss_fp);
				sprintf(str, "  <area shape=\"rect\" coords=\"%d,%d,%d,%d\" href=\"features_%s.html#%d\" title=\"%s\" />\n", (int)(z*seq->bouts[beh][i].start_frame), 0, 
					(int)(z*seq->bouts[beh][i].end_frame), (int)(z*h), fname2, i, alt);


				strcpy(ptr, str);
				ptr += strlen(str);
			}
			strcpy(ptr, "</map>");
			strcpy(html, html_tmp);
			free(html_tmp);
		}
	}

	return img;
}


