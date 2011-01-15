#include "svm_struct_api_blob_behavior_sequence.h"

#define BEH(f) ((!build_partial_label || (blobs->frames[f].is_good && (blobs->frames[f].is_manual&2)==2)) ? blobs->frames[f].behaviors[beh] : 0)


SVMFeatureParams g_blob_heavior_params[F_NUM_GOOD_GLOBAL_FEATURES] = { 
  SVMBlobBehaviorSequence::AspectRatioParams(),
  SVMBlobBehaviorSequence::VelocityParams(),
  SVMBlobBehaviorSequence::VelocityParams(),
  SVMBlobBehaviorSequence::MaxCurvatureParams(),
  SVMBlobBehaviorSequence::MaxCurvatureParams(),
  SVMBlobBehaviorSequence::CurveParams(),
  SVMBlobBehaviorSequence::CurveParams(),
  SVMBlobBehaviorSequence::CurveParams(),
  SVMBlobBehaviorSequence::CurveParams(),
  SVMBlobBehaviorSequence::CurveParams(),
  SVMBlobBehaviorSequence::CurveParams()
};

SVMFeatureParams SVMBlobBehaviorSequence::AspectRatioParams() {
  SVMFeatureParams p;
  memset(&p, 0, sizeof(SVMFeatureParams));
  p.feature_sample_smoothness_window = 5;
  
  // Compare the aspect ratio at each frame in the bout to the global min/max/average
  p.use_global_difference_max_sum_features = p.use_global_difference_min_sum_features = 
    p.use_global_difference_ave_sum_features = true;

  // Standard deviation and sum variance of the aspect ratio over the entire bout
  p.use_standard_deviation = p.use_sum_variance = true;

  // Use a 6D-histogram in feature space over the entire bout as a feature
  p.num_histogram_bins = 6;  
  p.num_histogram_temporal_levels = 1;
  p.use_histogram_sum_features = true;

  // Compare the aspect ratio at each frame in the bout to the aspect ratio in the frames immediately proceeding or following this bout
  p.num_difference_temporal_levels = 2;
  p.use_start_sum_absolute_diff_haar_features = p.use_end_sum_absolute_diff_haar_features = 
    p.use_start_sum_diff_haar_features = p.use_end_sum_diff_haar_features = true;

  return p;
}

SVMFeatureParams SVMBlobBehaviorSequence::VelocityParams() {
  SVMFeatureParams p;
  memset(&p, 0, sizeof(SVMFeatureParams));
  p.feature_sample_smoothness_window = 0;  // the velocity is already smoothed

  // Average velocity in the bout
  p.num_temporal_levels = 1;   
  p.use_bout_ave_features = true;

  // Total forward distance travelled over the duration of the bout
  p.use_bout_sum_features = true;

  // Standard deviation and sum variance of the velocity over the total bout
  p.use_standard_deviation = p.use_sum_variance = true;

  // Use a 6D-histogram in feature space over the entire bout as a feature
  p.num_histogram_bins = 6;  
  p.num_histogram_temporal_levels = 1;
  p.use_histogram_sum_features = true;

  // Compare the aspect ratio at each frame in the bout to the aspect ratio in the frames immediately proceeding or following this bout
  p.num_difference_temporal_levels = 2;
  p.use_start_sum_absolute_diff_haar_features = p.use_end_sum_absolute_diff_haar_features = 
    p.use_start_sum_diff_haar_features = p.use_end_sum_diff_haar_features = true;

  return p;
}

SVMFeatureParams SVMBlobBehaviorSequence::MaxCurvatureParams() {
  SVMFeatureParams p;
  memset(&p, 0, sizeof(SVMFeatureParams));
  p.feature_sample_smoothness_window = 5;

  // Bout-average curvature
  p.num_temporal_levels = 1;
  p.use_bout_ave_features = true;
  
  // Standard deviation and sum variance of the curvature over the total bout
  p.use_standard_deviation = p.use_sum_variance = true;

  // Bout maximum and minimum curvature
  p.use_bout_max_feature = p.use_bout_min_feature = true;

  // Use a 16D-histogram in feature space over the entire bout as a feature
  p.num_histogram_bins = 16;  
  p.num_histogram_temporal_levels = 1;
  p.use_histogram_sum_features = true;

  // Compare the curvature at each frame in the bout to the average curvature in the frames immediately proceeding or following this bout
  p.num_difference_temporal_levels = 2;
  p.use_start_sum_absolute_diff_haar_features = p.use_end_sum_absolute_diff_haar_features = 
    p.use_start_sum_diff_haar_features = p.use_end_sum_diff_haar_features = true;

  return p;
}

SVMFeatureParams SVMBlobBehaviorSequence::CurveParams() {
  return MaxCurvatureParams();
}


SVMBlobBehaviorSequence::SVMBlobBehaviorSequence(FitParams *p, struct _BehaviorGroups *behaviors, int beh, SVMFeatureParams *sparams) :
  SVMBehaviorSequence(F_NUM_GOOD_GLOBAL_FEATURES, behaviors, beh, sparams ? sparams : g_blob_heavior_params)
{
  params = *p;
}

int SVMBlobBehaviorSequence::num_frames(void *b) {
  return ((BlobSequence*)b)->num_frames;
}

char **SVMBlobBehaviorSequence::load_examples(const char *fname, int *num) {
  return load_train_list(fname, num);
}

/*
 * Use a sequence of behavior bouts to set the per-frame behavior labels in b
 */
void SVMBlobBehaviorSequence::load_from_bout_sequence(BehaviorBoutSequence *y, void *b) {
  int beh, i, j;
  BlobSequence *blobs = (BlobSequence*)b;

  for(beh = 0; beh < behaviors->num; beh++) {
    for(i = 0; i < y->num_bouts[beh]; i++) {
      for(j = y->bouts[beh][i].start_frame; j < y->bouts[beh][i].end_frame; j++) {
	if(!(blobs->frames[j].is_manual & 1))
	  blobs->frames[j].behaviors[beh] = y->bouts[beh][i].behavior;
      }
    }
  }
}

/*
 * Convert per frame behavior labels into a sequence of behavior bouts
 */
BehaviorBoutSequence *SVMBlobBehaviorSequence::create_behavior_bout_sequence(void *b, BehaviorGroups *behaviors, bool build_partial_label) {
  int i, j, beh, num = behaviors->num;
  unsigned int sz = sizeof(BehaviorBoutSequence)+(sizeof(int)+sizeof(BehaviorBout*)+2*sizeof(double))*num;
  BehaviorBoutSequence *s = (BehaviorBoutSequence*)my_malloc(sz);
  BlobSequence *blobs = (BlobSequence*)b;

  memset(s, 0, sz);
  s->bouts = (BehaviorBout**)(s+1);
  s->num_bouts = (int*)(s->bouts+num);
  s->scores = (double*)(s->num_bouts+num);
  s->losses = (double*)(s->scores+num);
  s->behaviors = behaviors;
  s->score = s->loss = s->slack = 0;

  for(beh = 0; beh < num; beh++) {
    i = 0;
    while(i < blobs->num_frames) {
      j = i;
      while(j < blobs->num_frames && BEH(j) == BEH(i))
	j++;
      if(BEH(i) >= 0) {
	s->bouts[beh] = (BehaviorBout*)realloc(s->bouts[beh], sizeof(BehaviorBout)*(s->num_bouts[beh]+1));
      
	s->bouts[beh][s->num_bouts[beh]].start_frame = i;
	s->bouts[beh][s->num_bouts[beh]].end_frame = j;
	s->bouts[beh][s->num_bouts[beh]].behavior = BEH(i);
	s->bouts[beh][s->num_bouts[beh]].bout_score = s->bouts[beh][s->num_bouts[beh]].transition_score = 0;
	s->bouts[beh][s->num_bouts[beh]].loss_fn = s->bouts[beh][s->num_bouts[beh]].loss_fp = 0;
	s->num_bouts[beh]++;
      }
      i = j;
    }
  }

  return s;
}

extern const char *g_global_feature_names[F_NUM_GOOD_GLOBAL_FEATURES];
const char *SVMBlobBehaviorSequence::get_base_feature_name(int ind) {
  return g_global_feature_names[ind];
}
 
/*
 * Store the blob-level (spine tracker output) features in b into the feature_cache->features
 */
void SVMBlobBehaviorSequence::load_behavior_bout_features(void *b, BehaviorBoutFeatures *feature_cache) {
  BlobSequence *blobs = (BlobSequence*)b;
  int i, j, T = blobs->num_frames;

  for(i = 0; i < T; i++) 
    if(!blobs->frames[i].features) 
      blobs->frames[i].features = (double*)my_malloc(2*sizeof(double)*((2*F_NUM_FEATURES)*(params.num_spine_lines+1)+F_NUM_GLOBAL_FEATURES));
  
  compute_deterministic_spine_attributes_blob_sequence(blobs, &params, 0);
  compute_all_global_features(blobs, &params, 0);

  // Extract regular raw features
  for(i = 0; i < num_base_features; i++) 
    for(j = 0; j < T; j++) 
      feature_cache->features[i][j] = GLOBAL_FEATURES(&blobs->frames[j])[i];
    
  for(j = 0; j < T; j++)
    feature_cache->frame_times[j] = blobs->frames[j].frame_time;
}

/*
 * Read a blob annotation file from disk.  Returns NULL if the specified file has no annotated frames
 */
void *SVMBlobBehaviorSequence::load_training_example(const char *fname, BehaviorGroups *behaviors) {
  Blob *b;
  BlobSequence *blobs = load_blob_sequence(fname, behaviors);  
  int has_frames = 0, i;

  for(i = 0; i < blobs->num_frames; i++) {
    b = &blobs->frames[i];
    if(b->is_good && (b->is_manual & 3) == 3 && b->num_model_pts == params.num_spine_lines+1) {
      has_frames = 1;
      break;
    }
  }
  if(!has_frames) {
    free_blob_sequence(blobs);
    return NULL;
  }

  return blobs;
}

/*
 * Read a blob annotation file from disk.  Returns NULL if the specified file has no annotated frames
 */
void SVMBlobBehaviorSequence::save_example(void *b, void *d, const char *fname) {
  load_from_bout_sequence((BehaviorBoutSequence*)d, b);
  save_blob_sequence(fname, (BlobSequence*)b, behaviors, params.num_spine_lines+1, params.num_orientations) ;
}

/*
 * free a blob
 */
void SVMBlobBehaviorSequence::free_data(void *b) {
  free_blob_sequence((BlobSequence*)b);
}
