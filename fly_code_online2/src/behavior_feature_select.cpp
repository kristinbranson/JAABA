#include "svm_behavior_sequence.h"

#define MAX_MEMORY_BYTES 2e9


typedef struct _DecisionStump {
  int feat;
  double thresh;

  double score;
  double *weights;
  double *scores;
  double *sample_features;
  int num_samples;
} DecisionStump;

typedef struct _ExampleWeights {
  double weight;
  double alpha;
  double f;
  int ex_ind;
  int samp_ind;
  int bout_ind;
  bool is_gt;
  BehaviorBout bout;
} ExampleWeights;

typedef enum {
  WS_MINIMIZE_DUAL,
  WS_FUNCTIONAL_GRADIENT_DESCENT
} WeightingScheme;


int ExampleWeights_cmp_f(const void *a, const void *b) {
  double d = ((ExampleWeights*)a)->f - ((ExampleWeights*)b)->f;
  //if(((ExampleWeights*)a)->alpha > 0 && ((ExampleWeights*)b)->alpha == 0) return -1;
  //else if((ExampleWeights*)b)->alpha > 0 && (ExampleWeights*)a)->alpha == 0) return 1;
  return d < 0 ? -1 : (d > 0 ? 1 : 0);
}

// Create a sample for each predicted bout and ground truth bout, and weight each sample
ExampleWeights *SVMBehaviorSequence::ComputeExampleWeights(StructuredDataset *trainset, int *numEx) {
  int i, j, k, numTotal = 0, num = 0;
  for(i = 0; i < trainset->num_examples; i++) {
    if(!trainset->examples[i]->set) continue;
    UncondenseSamples(trainset->examples[i]->set);
    BehaviorBoutSequence *y_gt = (BehaviorBoutSequence*)trainset->examples[i]->y;
    numTotal += y_gt->num_bouts;
    for(j = 0; j < trainset->examples[i]->set->num_samples; j++) {
      BehaviorBoutSequence *y = (BehaviorBoutSequence*)trainset->examples[i]->set->samples[j].ybar;
      numTotal += y->num_bouts;
    }
  }
  ExampleWeights *r = new ExampleWeights[numTotal];

  WeightingScheme weighting_scheme = WS_MINIMIZE_DUAL;  // TODO: make this a parameter
  if(weighting_scheme == WS_MINIMIZE_DUAL) {
    // Weighting scheme will result in choosing the next feature that has maximum potential to 
    // reduce the dual objective
    for(i = 0, num = 0; i < trainset->num_examples; i++) {
      int num_old = num;
      if(!trainset->examples[i]->set) continue;
      BehaviorBoutSequence *y_gt = (BehaviorBoutSequence*)trainset->examples[i]->y, *y_partial = NULL;
      fill_unlabeled_gt_frames(y_gt, y_partial);
      for(k = 0; k < y_gt->num_bouts; k++) {
	if((!NONE_CLASS_HAS_NO_SCORE || y_gt->bouts[k].behavior != 0)) {
	  r[num].bout = y_gt->bouts[k];
	  r[num].ex_ind = i;
	  r[num].samp_ind = -1;
	  r[num].bout_ind = k;
	  r[num].is_gt = true;
	  r[num].alpha = -trainset->examples[i]->set->alpha;
	  r[num++].weight = -trainset->examples[i]->set->alpha;
	}
      }
      for(j = 0; j < trainset->examples[i]->set->num_samples; j++) {
	BehaviorBoutSequence *y = (BehaviorBoutSequence*)trainset->examples[i]->set->samples[j].ybar;
	for(k = 0; k < y->num_bouts; k++) {
	  if((!NONE_CLASS_HAS_NO_SCORE || y->bouts[k].behavior != 0)) {
	    r[num].bout = y->bouts[k];
	    r[num].ex_ind = i;
	    r[num].samp_ind = j;
	    r[num].bout_ind = k;
	    r[num].is_gt = false;
	    r[num].alpha = trainset->examples[i]->set->samples[j].alpha;
	    r[num++].weight = trainset->examples[i]->set->samples[j].alpha;
	  }
	}
      }
      if(y_gt != trainset->examples[i]->y)
	delete y_gt;
      if(y_partial)
	delete y_partial;
      //fprintf(stderr, " %d:%d:%d", i, num-num_old, trainset->examples[i]->set->num_samples);
    }
  } else if(weighting_scheme == WS_FUNCTIONAL_GRADIENT_DESCENT) {
    // Weighting scheme will result in choosing the next feature that results in the steepest descent
    // of the primal objective
    for(i = 0, num = 0; i < trainset->num_examples; i++) {
      if(!trainset->examples[i]->set) continue;
      SVM_cached_sample_set_compute_features(trainset->examples[i]->set, trainset->examples[i]);
      double best = -1;
      BehaviorBoutSequence *y_gt = (BehaviorBoutSequence*)trainset->examples[i]->y, *y_partial = NULL;
      fill_unlabeled_gt_frames(y_gt, y_partial);
      int num2 = num;
      for(k = 0; k < y_gt->num_bouts; k++) {
	if((!NONE_CLASS_HAS_NO_SCORE || y_gt->bouts[k].behavior != 0)) {
	  r[num2].bout = y_gt->bouts[k];
	  r[num2].ex_ind = i;
	  r[num2].samp_ind = j;
	  r[num2].bout_ind = k;
	  r[num2].is_gt = true;
	  r[num2].alpha = -trainset->examples[i]->set->alpha;
	  r[num2++].weight = -1;
	}
      }
      for(j = 0; j < trainset->examples[i]->set->num_samples; j++) {
	double slack = sum_w->dot(*trainset->examples[i]->set->samples[j].psi) - 
	  trainset->examples[i]->set->score_gt*sum_w_scale;
	if(slack > best && slack > 0) {
	  BehaviorBoutSequence *y = (BehaviorBoutSequence*)trainset->examples[i]->set->samples[j].ybar;
	  for(k = 0; k < y->num_bouts; k++) {
	    if(!NONE_CLASS_HAS_NO_SCORE || y->bouts[k].behavior != 0) {
	      r[num2].bout = y->bouts[k];
	      r[num2].ex_ind = i;
	      r[num2].samp_ind = j;
	      r[num2].bout_ind = k;
	      r[num2].is_gt = false;
	      r[num2].alpha = trainset->examples[i]->set->samples[j].alpha;
	      r[num2++].weight = 1;
	    }
	  }
	  best = slack;
	}
      }
      if(best > 0) num=num2;
      if(y_gt != trainset->examples[i]->y)
	delete y_gt;
      if(y_partial)
	delete y_partial;
      OnFinishedIteration(trainset->examples[i]->x, trainset->examples[i]->y);
    }
  }

  *numEx = num;

  return r;
}




// Consider adding a new bout feature, which will add behaviors->num_values new features (one for each class).  The new bout feature 
// will be a decision stump applied to a bout feature drawn from the array bout_expansion_features.  The score of a candidate bout
// feature and threshold will be the sum of squares score over the induced features, and we pick the best one
DecisionStump *SVMBehaviorSequence::SelectBestFeature(ExampleWeights *samples, int num_samples, double **f) {
  int last_example = -1;
  ExampleWeights *v = new ExampleWeights[num_samples];
  double *scores = new double[behaviors->num_values*2];
  double *weights = scores + behaviors->num_values;
  int chunk_size = (int)(MAX_MEMORY_BYTES/sizeof(double)/num_samples);

  DecisionStump *retval = (DecisionStump*)malloc(sizeof(DecisionStump)+(behaviors->num_values*2+num_samples)*sizeof(double));
  retval->scores = (double*)(retval+1);
  retval->weights = retval->scores + behaviors->num_values;
  retval->sample_features = retval->weights + behaviors->num_values;
  retval->thresh = -1;
  retval->feat = -1;
  retval->score = -1;
  retval->num_samples = num_samples;
    
  // Loop through each feature, and compute the highest scoring decision stump for all possible thresholds
  for(int i = 0; i < num_bout_expansion_features; i+=chunk_size) {
    int n = my_min(chunk_size, num_bout_expansion_features-i);

    // extract features
    for(int j = 0; j < num_samples; j++) {
      if(samples[j].ex_ind != last_example && last_example != -1)
	OnFinishedIteration(trainset->examples[last_example]->x, trainset->examples[last_example]->y);
      last_example = samples[j].ex_ind;
      bool noFeat = ((!samples[j].is_gt && !trainset->examples[samples[j].ex_ind]->set->samples[samples[j].samp_ind].psi) || 
		      (samples[j].is_gt && !trainset->examples[samples[j].ex_ind]->set->psi_gt));

      if(!f[j] || chunk_size < num_bout_expansion_features || noFeat) {
	BehaviorBoutFeatures *x = (BehaviorBoutFeatures*)trainset->examples[samples[j].ex_ind]->x;
	if(!x->memory_buffer)
	  x->ComputeCaches(this);
	if(noFeat)
	  SVM_cached_sample_set_compute_features(trainset->examples[samples[j].ex_ind]->set, trainset->examples[samples[j].ex_ind]);
	if(f[j]) delete [] f[j];
	f[j] = psi_bout(x, samples[j].bout.start_frame, samples[j].bout.end_frame, samples[j].bout.behavior, NULL, false, false, 
			bout_expansion_features+i, n);
      }
    }
      
    for(int ii=0; ii < n; ii++) {
      for(int j = 0; j < num_samples; j++) {
	v[j] = samples[j];
	v[j].f = f[j][ii];
      }

      // Sort samples by feature value, such that if we consider a threshold t, for 
      // some value of i, \forall_{i<=j} v[i].f<t and \forall_{i>j} v[i].f>=t
      qsort(v, num_samples, sizeof(ExampleWeights), ExampleWeights_cmp_f);
      double score = 0;
      for(int k=0; k < behaviors->num_values; k++)
	scores[k] = weights[k] = 0;

      // Choose the decision stump threshold on feature i+ii that gives the highest score
      int updatedBest = -1;
      for(int j = 0; j < num_samples-1; j++) {
	score -= SQR(scores[v[j].bout.behavior]);
	scores[v[j].bout.behavior] += v[j].weight;   // scores[b] = sum_{i=1}^n ((f_psi<t)-(f_psi_gt<t))*samples[i].weight
	score += SQR(scores[v[j].bout.behavior]);
	weights[v[j].bout.behavior] -= v[j].alpha;   // weights[b] = sum_{i=1}^n ((f_psi<t)-(f_psi_gt<t))*samples[i].alpha
	if(score > retval->score && v[j].f != v[j+1].f) {  
	  updatedBest = j;
	  retval->score = score;
	  retval->feat = i+ii;

	  // choose a threshold t such that \forall_{i<=j} v[i].f<t and \forall_{i>j} v[i].f>=t
	  retval->thresh = (v[j].f+v[j+1].f)/2;  
	}
      }
      if(updatedBest >= 0) {
	for(int j = 0; j < behaviors->num_values; j++) 
	  retval->weights[j] = retval->scores[j] = 0;
	for(int j = 0; j <= updatedBest; j++) {
	  retval->scores[v[j].bout.behavior] += v[j].weight;
	  retval->weights[v[j].bout.behavior] -= v[j].alpha;
	}
	for(int j = 0; j < num_samples; j++) 
	  retval->sample_features[j] = f[j][ii];
      }
    }
  }

  delete [] v;
  delete [] scores;

  return retval;
}


void SVMBehaviorSequence::AugmentFeatureSpace() {
  PauseWorkerThreads(true, true, false);
  int num_samples=0, num_samples_old = -1;
  double **f = NULL;
  do {
    fprintf(stderr, "Augment feature space...");
    ExampleWeights *samples = ComputeExampleWeights(trainset, &num_samples);
    if(num_samples_old >= 0) assert(num_samples == num_samples_old);
    num_samples_old = num_samples;
    if(!f) {
      f = new double*[num_samples];
      memset(f, 0, sizeof(double*)*num_samples);
    }
    DecisionStump *s = SelectBestFeature(samples, num_samples, f);
    fprintf(stderr, "Adding feature %s<%f with weight %f\n", bout_expansion_features[s->feat].name, (float)s->thresh, s->weights[1]);
    AppendNewFeature(s, samples, num_samples);
    free(s);
    OptimizeAllConstraints(40);
    delete [] samples;
  } while(sum_dual/n > maxDual);
  for(int j = 0; j < num_samples; j++) 
    delete [] f[j];
  delete [] f;
  PauseWorkerThreads(false, false, false);
}

void SVMBehaviorSequence::AppendNewFeature(DecisionStump *s, ExampleWeights *samples, int num_samples) {
  bout_features = (BoutFeature*)realloc(bout_features, (num_bout_features+1)*sizeof(BoutFeature));
  char name[1000];
  sprintf(name, "%s<%.3f", bout_expansion_features[s->feat].name, bout_expansion_features[s->feat].thresh);
  bout_features[num_bout_features] = bout_expansion_features[s->feat];
  bout_features[num_bout_features].name = StringCopy(name);
  bout_features[num_bout_features].thresh = s->thresh; 
  bout_features[num_bout_features].num_thresholds = 1;
  bout_features[num_bout_features].mu = 0;
  bout_features[num_bout_features].gamma = 1;

  for(int ii = 0; ii < behaviors->num_values; ii++) {
    double w = s->weights[ii];
    int ind = (num_bout_features+1)*ii + num_bout_features;
    sum_w->InsertAt(ind, w);
    sum_w_sqr += SQR(w);
    regularization_error = sum_w_sqr/SQR(sum_w_scale)*lambda/2;
    sum_dual = -sum_w_sqr/(2*sum_w_scale) + sum_alpha_loss;
  
    // Augment existing cached feature spaces 
    int ss = 0;
    for(int i = 0; i < trainset->num_examples; i++) {
      SVM_cached_sample_set *set = trainset->examples[i]->set;
      if(!set) continue;
      if(set->psi_gt) {
	double f = 0;
	BehaviorBoutSequence *y_gt = (BehaviorBoutSequence*)trainset->examples[i]->y;
	for(int k = 0; k < y_gt->num_bouts; k++) 
	  if(!NONE_CLASS_HAS_NO_SCORE || y_gt->bouts[k].behavior != 0) {
	    if(y_gt->bouts[k].behavior == ii) 
	      f += s->sample_features[ss] < s->thresh;
	    ss++;
	  }
	set->psi_gt->InsertAt(ind, f);
      } else {
	if(set->psi_gt) delete set->psi_gt;
	set->psi_gt = NULL;
      }
      
      for(int j = 0; j < set->num_samples; j++) {
	if(set->samples[j].psi) {
	  double f = 0;
	  BehaviorBoutSequence *y = (BehaviorBoutSequence*)set->samples[j].ybar;
	  for(int k = 0; k < y->num_bouts; k++) 
	    if(!NONE_CLASS_HAS_NO_SCORE || y->bouts[k].behavior != 0) {
	      if(y->bouts[k].behavior == ii) 
		f += s->sample_features[ss] < s->thresh;
	      ss++;
	    }
	  set->samples[j].psi->InsertAt(ind, f);
	} else {
	  if(set->samples[j].psi) delete set->samples[j].psi;
	  set->samples[j].psi = NULL;
	  set->samples[j].sqr = 0;
	}
      }
      SVM_cached_sample_set_recompute_caches(set);
   
    }
    assert(ss == s->num_samples);
  }

  num_bout_features++;
  sizePsi += behaviors->num_values;
  u_i_buff = (double*)realloc(u_i_buff, sizeof(double)*sizePsi);
  for(int i = 0; i < sizePsi; i++) u_i_buff[i] = 0;
}





