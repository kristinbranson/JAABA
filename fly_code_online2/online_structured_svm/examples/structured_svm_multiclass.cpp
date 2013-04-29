#include "structured_svm_multiclass.h"

/**
 * @file structured_svm_multiclass.cpp
 * @brief Simple example of how to use this structured SVM package: implements routines for loss-sensitive multiclass SVM training and classification
 */



/**
 * @example structured_svm_multiclass.cpp
 *
 * This is an example of how to use the structured learning API to implement a custom structured learner.  This
 * example implements a multiclass SVM learner and classification with custom loss function
 *
 * Example usage:
 * - Train using a fixed dataset without running in server mode, outputting the learned model to learned_model.txt.  classes.txt
 *    defines the confusion cost between classes and is in the format of MulticlassStructuredSVM::Load, and train.txt is in the
 *    format of MulticlassStructuredSVM::LoadDataset
\htmlonly <div style="padding: 0.5em 1em; border-top: 1px solid #ddd; border-bottom: 1px solid #ddd; background-color: #eaeafa;">
$ examples/bin/release_static/structured_svm_multiclass.out -p classes.txt -d train.txt -o learned_model.txt
</div> \endhtmlonly
 * - Evaluate performance on a testset:
\htmlonly <div style="padding: 0.5em 1em; border-top: 1px solid #ddd; border-bottom: 1px solid #ddd; background-color: #eaeafa;">
$ examples/bin/release_static/structured_svm_multiclass.out -p learned_model.txt -t test.txt test.txt.predictions
</div> \endhtmlonly
 *
 */

#define DONT_STORE_SAMPLES


StructuredLabel *MulticlassStructuredSVM::NewStructuredLabel(StructuredData *x) { return new MulticlassStructuredLabel(x); }

StructuredData *MulticlassStructuredSVM::NewStructuredData() { return new MulticlassStructuredData; }


MulticlassStructuredSVM::MulticlassStructuredSVM() {
  num_classes = 0;
  num_features = 0;
  classConfusionCosts = NULL;
  eps = .01;
  C = 5000.0;
  window = 1000;
  canScaleW = true;
  mergeSamples = false;
  alphas = NULL;
  alphas_alloc_size = 0;
  numMultiSampleIterations = 10;
}

MulticlassStructuredSVM::~MulticlassStructuredSVM() {
  if(classConfusionCosts) 
    free(classConfusionCosts);
  if(alphas) {
    for(int i = 0; i < alphas_alloc_size; i++)
      if(alphas[i]) 
	delete [] alphas[i];
    free(alphas);
  }
}

double MulticlassStructuredSVM::Inference(StructuredData *x, StructuredLabel *ybar, SparseVector *w,
					  StructuredLabel *y_partial, StructuredLabel *y_gt, double w_scale) {
  int bestclass=-1, first=1;
  double score,bestscore=-1;

  MulticlassStructuredLabel *m_ybar = (MulticlassStructuredLabel*)ybar;
  MulticlassStructuredLabel *m_y_partial = y_partial ? (MulticlassStructuredLabel*)y_partial : NULL;

  // Loop through every possible class y and compute its score <w,Psi(x,y)>
  for(int class_id = 1; class_id <= num_classes; class_id++) {
    // By default, compute ybar = max_y <w,Psi(x,y)>, but it y_partial is non-null,
    // only consider labels that agree with y_partial 
    m_ybar->class_id=class_id;
    if(y_partial && m_y_partial->class_id != class_id)
      score = -INFINITY;
    else 
      score = w->dot(Psi(x, ybar))*w_scale;
   
    // If y_gt is non-null, compute ybar = max_y <w,Psi(x,y)>+Loss(y_gt,y) 
    if(y_gt)  
      score += Loss(y_gt, ybar);

    if(score > bestscore || first) {
      bestscore=score;
      bestclass=class_id;
      first=0;
    }
  }
  m_ybar->class_id = bestclass;

  return bestscore;
}

void MulticlassStructuredSVM::CreateSamples(struct _SVM_cached_sample_set *set, StructuredData *x, StructuredLabel *y_gt) {
  if(!set->num_samples) {
    // Add a sample set that includes all classes 
    for(int i = 0; i < num_classes; i++) {
      if(i+1 == ((MulticlassStructuredLabel*)y_gt)->class_id) continue;
      MulticlassStructuredLabel *ybar = (MulticlassStructuredLabel*)NewStructuredLabel(x);
      ybar->class_id = i+1;
      SVM_cached_sample *sample = SVM_cached_sample_set_add_sample(set, ybar);

      // Optionally set these things, so they don't have to be computed later
      sample->psi = Psi(x, ybar).ptr();
      sample->loss = Loss(y_gt, ybar);
      sample->sqr = 2*set->psi_gt_sqr;   // <set->psi_gt-sample->psi,set->psi_gt-sample->psi>
      sample->dot_psi_gt_psi = 0;  // <set->psi_gt,sample->psi>
    }
  }
}

double MulticlassStructuredSVM::ImportanceSample(StructuredData *x, SparseVector *w, StructuredLabel *y_gt, struct _SVM_cached_sample_set *set, double w_scale) {
  int first=1;
  double score,bestscore=0;
  int cl = ((MulticlassStructuredLabel*)y_gt)->class_id;

  CreateSamples(set, x, y_gt);

  for(int i = 0, j = 0; i < num_classes; i++) {
    if(i+1 == ((MulticlassStructuredLabel*)y_gt)->class_id) continue;
    SVM_cached_sample *sample = &set->samples[j++];
    score = w->dot(*sample->psi)*w_scale + sample->loss;
    sample->slack = score-set->score_gt;
    sample->dot_w = (sample->slack - sample->loss)*sum_w_scale;
    if(score > bestscore || first) {
      bestscore=score;
      first=0;
      cl = i+1;
    }
  }
  set->dot_sum_w_psi_gt = set->score_gt*sum_w_scale;
  qsort(set->samples, set->num_samples, sizeof(SVM_cached_sample), SVM_cached_sample_cmp);

#ifdef DONT_STORE_SAMPLES
  if(set->i+1 > alphas_alloc_size) {
    alphas = (double**)realloc(alphas, sizeof(double*)*(set->i+1));
    for(int i = alphas_alloc_size; i < set->i+1; i++)
      alphas[i] = NULL;
    alphas_alloc_size = set->i+1;
  }
#endif
  return bestscore;
}


// Iteratively optimize alpha parameters for a set of samples 's'.  Implements 'Multi-Sample Dual Update Step' of writeup
void MulticlassStructuredSVM::MultiSampleUpdate(SVM_cached_sample_set *set, StructuredExample *ex, int R) {
  VFLOAT dot_w_u = 0;   // dot_w_u = <w,u>
  VFLOAT L_i = 0;       // L_i = <w,-u_i/(lambda*t)> + \sum_ybar alpha_{i,ybar} loss(y_i,ybar)
  int i, r;
  double *alpha = new double[set->num_samples]; 

#ifdef DONT_STORE_SAMPLES
  CreateSamples(set, ex->x, ex->y);
  if(alphas[set->i])
    for(i = 0; i < set->num_samples; i++) 
      set->samples[i].alpha = alphas[set->i][((MulticlassStructuredLabel*)set->samples[i].ybar)->class_id-1];
#endif

  if(!set->num_samples && !set->alpha) return;

  for(i = 0; i < set->num_samples; i++) {
    alpha[i] = set->samples[i].alpha;
    dot_w_u += alpha[i]*set->samples[i].dot_w/sum_w_scale;
  }
  L_i = set->D_i+dot_w_u;

  VFLOAT D_i = set->D_i;
  VFLOAT sum_alpha = set->alpha;

  for(r = 0; r < R; r++) {
    double s_u = 1;
    for(i = 0; i < set->num_samples; i++) {
      SVM_cached_sample *s = &set->samples[i];

      VFLOAT dot_u_v = set->psi_gt_sqr*(alpha[i]*s_u+sum_alpha);  // <u_i,v>
      VFLOAT dot = s->dot_w - set->psi_gt_sqr*((alpha[i]*s_u-s->alpha) + sum_alpha-set->alpha);

      VFLOAT scale=1, new_sum_alpha;
      VFLOAT dalpha = (dot + s->loss*(sum_w_scale)) / my_max(s->sqr,.0000000001);
      if(sum_alpha < 1 || dalpha < 0) {
	// alpha expand: solve for the optimal amount 'dalpha' to increase s->alpha
	// and then scale all set->samples[:].alpha (the value of 'dalpha' and 'scale' 
	// that maximizes the increase in the dual objective). 
	dalpha = my_min(1-alpha[i]*s_u, my_max(-alpha[i]*s_u,dalpha));
	if(set->u_i_sqr && sum_alpha) {
	  scale = 1 + (L_i*sum_w_scale - dalpha*dot_u_v) / set->u_i_sqr;
	  scale = my_min((1-dalpha)/sum_alpha, my_max(0, scale));
	  if(alpha[i]*s_u*scale+dalpha < 0) dalpha = -alpha[i]*s_u*scale;
	}
	new_sum_alpha = sum_alpha*scale + dalpha;
      } else {
	// alpha swap: solve for the optimal amount 'dalpha' to increase s->alpha
	// while scaling down all set->samples[:].alpha, such that we preserve sum_k{s->samples[:].alpha}=1
	// (choose the value of 'dalpha' that maximizes the increase in the dual objective)
	VFLOAT e = dot/(sum_w_scale) + s->loss;
	VFLOAT sqr = set->u_i_sqr + 2*dot_u_v + s->sqr;
	dalpha = (e-L_i)*(sum_w_scale) / my_max(sqr,.00000000001);
	dalpha = my_min(1-alpha[i]*s_u, my_max(-alpha[i]*s_u,dalpha));
	scale = 1-dalpha;
	//assert(scale > 0 && scale <= 1);
	new_sum_alpha = 1;
      }
      assert(scale >= 0 && new_sum_alpha >= -0.000000001 && new_sum_alpha <= 1.000000001);
      new_sum_alpha = my_min(1,my_max(new_sum_alpha,0));

      if(dalpha != 0 || scale != 1) {
	s_u *= scale;
	if(s_u < .01 || s_u > 10) {
	  for(int j = 0; j < set->num_samples; j++) alpha[j] *= s_u;
	  s_u = 1;
	}
	alpha[i] += dalpha/s_u;

	// Keep track of L_i, D_i, u_i_sqr, dot_w_u, dot_u_psi_gt using inexpensive online updates
	sum_alpha = new_sum_alpha;
        dot_w_u = scale*dot_w_u + (dalpha*dot - scale*(scale-1)*set->u_i_sqr - (2*scale*dalpha-dalpha)*dot_u_v - s->sqr*SQR(dalpha)) / (sum_w_scale);
	set->u_i_sqr = SQR(scale)*set->u_i_sqr + 2*scale*dalpha*dot_u_v + s->sqr*SQR(dalpha);
	set->dot_u_psi_gt = scale*set->dot_u_psi_gt + dalpha*(s->dot_psi_gt_psi - set->psi_gt_sqr);
	D_i = scale*D_i + dalpha*s->loss;
	L_i = dot_w_u + D_i;
	assert(!isnan(L_i));
      }
    }
    if(s_u != 1) 
      for(i = 0; i < set->num_samples; i++)
	alpha[i] *= s_u;
  }

  set->u_i_sqr = 0;
  double d_sum_w_sqr = 0;
  D_i = 0;
  for(i = 0; i < set->num_samples; i++) {
    double dalpha = alpha[i] - set->samples[i].alpha;
    if(dalpha) {
      d_sum_w_sqr += -2*dalpha*(set->samples[i].dot_w+set->dot_sum_w_psi_gt) + set->psi_gt_sqr*SQR(dalpha);
      *sum_w -= (*set->samples[i].psi*dalpha);
    }
    set->samples[i].alpha = alpha[i];
    set->u_i_sqr += set->psi_gt_sqr*(SQR(alpha[i]) + SQR(sum_alpha));
    D_i += alpha[i]*set->samples[i].loss;
  }
  set->dot_u_psi_gt = -sum_alpha*set->psi_gt_sqr;
  double dalpha = sum_alpha-set->alpha;
  if(dalpha) {
    d_sum_w_sqr += 2*dalpha*set->dot_sum_w_psi_gt + set->psi_gt_sqr*SQR(dalpha);
    *sum_w += (*set->psi_gt*dalpha);
  }
  set->alpha = sum_alpha;
  sum_alpha_loss += D_i - set->D_i;
  sum_w_sqr += d_sum_w_sqr;
  sum_dual = -sum_w_sqr/(2*sum_w_scale) + sum_alpha_loss;
  regularization_error = (sum_w_sqr/SQR(sum_w_scale))*lambda/2;
  set->D_i = D_i;


  set->slack_after = set->num_samples ? (sum_w->dot(*set->samples[0].psi)-sum_w->dot(*set->psi_gt))/(sum_w_scale) + set->samples[0].loss : 0;
  
  /*
  fprintf(stderr, "%f: ", (float)dalpha);
  for(int i = 0; i < set->num_samples; i++)
    fprintf(stderr, " %f:%f->%f", (float)(alpha[i]/sum_w_scale), (float)set->samples[i].slack, (float)(sum_w->dot(*set->samples[i].psi-*set->psi_gt)/(sum_w_scale) + set->samples[i].loss));
  */

#ifdef DONT_STORE_SAMPLES
  if(!alphas[set->i]) {
    alphas[set->i] = new double[num_classes];
    for(int i = 0; i < num_classes; i++)
      alphas[set->i][i] = 0;
  }
  for(int i = 0; i < set->num_samples; i++) {
    alphas[set->i][((MulticlassStructuredLabel*)set->samples[i].ybar)->class_id-1] = set->samples[i].alpha;
    free_SVM_cached_sample(&set->samples[i]);
  }
  if(set->samples) free(set->samples);
  set->samples = NULL;
  set->num_samples = 0;
#endif
  delete [] alpha;

  //fprintf(stderr, "sum dual is %lg, dual_change=%lg, D=%lg\n", sum_dual, -d_sum_w_sqr/(2*sum_w_scale) + set->D_i - D_i_orig, set->D_i);
}


SparseVector MulticlassStructuredSVM::Psi(StructuredData *x, StructuredLabel *y) {
  // The dimensionality of Psi(x,y) is num_featuresXnum_classes, by concatenating
  // num_features features for each class. The entries for Psi are equal to x->psi for 
  // the true class y and 0 for all classes other than y.  
  MulticlassStructuredData *m_x = (MulticlassStructuredData*)x;
  MulticlassStructuredLabel *m_y = (MulticlassStructuredLabel*)y;
  return m_x->psi->shift((m_y->class_id-1)*num_features);
}

double MulticlassStructuredSVM::Loss(StructuredLabel *y_gt, StructuredLabel *y_pred) {
  // Computes the loss of prediction y_pred against the correct label y_gt. 
  MulticlassStructuredLabel *m_y_gt = (MulticlassStructuredLabel*)y_gt;
  MulticlassStructuredLabel *m_y_pred = (MulticlassStructuredLabel*)y_pred;
  if(classConfusionCosts)
    return classConfusionCosts[m_y_gt->class_id][m_y_pred->class_id];
  else 
    return m_y_gt->class_id == m_y_pred->class_id ? 0 : 1;
}

Json::Value MulticlassStructuredSVM::Save() {
  Json::Value root;
  
  root["version"] = VERSION;
  root["Num Classes"] = num_classes;
  root["Num Features"] = num_features;

  //if(classConfusionCosts) {
    Json::Value c;
    int n = 0;
    for(int i = 1; i <= num_classes; i++) {
      for(int j = 1; j <= num_classes; j++) {
	if((!classConfusionCosts && i != j) || (classConfusionCosts && classConfusionCosts[i][j])) {
	  Json::Value o;
	  o["c_gt"] = i;
	  o["c_pred"] = j;
	  o["loss"] = classConfusionCosts ? classConfusionCosts[i][j] : 1;
	  c[n++] = o;
	}
      }
    }
    root["Class Confusion Costs"] = c;
  //}

    char fname[1000];
    sprintf(fname, "%s.export", modelfile);
    ExportModel(fname);

  return root;
}


//#define EXPORT_WEIGHTS_SVM_LIGHT_FORMAT
void MulticlassStructuredSVM::ExportModel(const char *fname) {
  FILE *fout = fopen(fname, "w"); 
  if(!fout) { fprintf(stderr, "Couldn't open %s for writing\n", fname); return; }

  SparseVector *w = GetCurrentWeights(false);
  w->make_non_sparse(false);
  char *str = new char[3000000];
  for(int i = 0; i < num_classes; i++) {
    SparseVector class_weights = w->extract_subset(i*num_features, (i+1)*num_features);
#ifdef EXPORT_WEIGHTS_SVM_LIGHT_FORMAT
    char *str = class_weights.to_string();
    fprintf(fout, "%d %s\n", i+1, str);
    free(str);
#else
    double *w_c =class_weights.get_non_sparse<double>(num_features);
    assert(!w_c[0]);
    for(int j = 1; j < num_features; j++)
      fprintf(fout, "%lg\n", w_c[j]);
#endif
  }

  fclose(fout);
  delete [] str;
  delete w;
}

bool MulticlassStructuredSVM::Load(const Json::Value &root) {
  fprintf(stdout, "loading parameters\n");
  if(strcmp(root.get("version", "").asString().c_str(), VERSION)) {
    fprintf(stderr, "Version of parameter file does not match version of the software"); 
    return false;
  }
  num_classes = root.get("Num Classes",0).asInt();
  num_features = root.get("Num Features",0).asInt();
  
  sizePsi = num_features*num_classes;

  if(root.isMember("Class Confusion Costs") && root["Class Confusion Costs"].isArray()) {
    classConfusionCosts = (double**)malloc((num_classes+1)*(sizeof(double*)+(num_classes+1)*sizeof(double)));
    double *ptr = (double*)(classConfusionCosts+(num_classes+1));
    for(int i = 0; i <= num_classes; i++, ptr += (num_classes+1)) {
      classConfusionCosts[i] = ptr;
      for(int j = 0; j <= num_classes; j++)
	classConfusionCosts[i][j] = 0;
    }
    Json::Value a = root["Class Confusion Costs"];
    for(int i = 0; i < (int)a.size(); i++) {
      int c_gt = a[i].get("c_gt",-1).asInt();
      int c_pred = a[i].get("c_pred",-1).asInt();
      double l = a[i].get("loss",0).asDouble();
      if(c_gt > num_classes || c_pred > num_classes || c_gt <= 0 || c_pred <= 0) {
	fprintf(stderr, "Error reading Class Confusion Costs\n");
	return false;
      }
      classConfusionCosts[c_gt][c_pred] = l;
    }
  }

  return true;
}



// Ordinarily, one need not override this function and should use the default StructuredSVM::LoadDataset() function
// instead.  We override because we want to import SVM^light format files
StructuredDataset *MulticlassStructuredSVM::LoadDataset(const char *fname) {
  if(debugLevel > 0) fprintf(stderr, "Reading dataset %s...", fname);
  
  Lock();

  bool detectNumFeatures = 1;//num_features == 0;
  bool detectNumClasses = num_classes == 0;

  FILE *fin = fopen(fname, "r");
  if(!fin) {
    fprintf(stderr, "Couldn't open dataset file %s\n", fname);
    Unlock();
    return NULL;
  }

  StructuredDataset *d = new StructuredDataset();
  char *line = new char[1000000];
  
  while(fgets(line, 99999, fin) && strlen(line) > 1) {
    chomp(line);
    StructuredExample *ex = new StructuredExample;
    ex->x = NewStructuredData();
    ex->y = NewStructuredLabel(ex->x);
    SparseVector *psi = ((MulticlassStructuredData*)ex->x)->psi = new SparseVector();

    // Assume each example is a line containing a class id followed by %d:%f pairs 
    if(!sscanf(line, "%d", &((MulticlassStructuredLabel*)ex->y)->class_id) || !psi->from_string(strstr(line, " ")+1)) { 
      fprintf(stderr, "Error parsing dataset example %s\n", line);
      delete ex;
      delete d;
      fclose(fin);
      return false;
    }

    if(detectNumFeatures) num_features = my_max(num_features, psi->Length());
    if(detectNumClasses) num_classes = my_max(num_classes, ((MulticlassStructuredLabel*)ex->y)->class_id);

    d->AddExample(ex);
  }
  delete [] line;
  fclose(fin);
  Unlock();

  if(detectNumFeatures || detectNumClasses)
    sizePsi = num_features*num_classes;

  if(debugLevel > 0) fprintf(stderr, "done\n");

  return d;
}

// Ordinarily, one need not override this function and should use the default StructuredSVM::SaveDataset() function
// instead.  We override because we want to import SVM^light format files
bool MulticlassStructuredSVM::SaveDataset(StructuredDataset *d, const char *fname, int start_from) {
  if(debugLevel > 0 && start_from == 0) fprintf(stderr, "Saving dataset %s...", fname);

  Lock();

  FILE *fout = fopen(fname, start_from>0 ? "a" : "w");
  if(!fout) {
    fprintf(stderr, "Couldn't open dataset file %s for writing\n", fname);
    Unlock();
    return false;
  }

  for(int i = start_from; i < d->num_examples; i++) {
    char *data = ((MulticlassStructuredData*)d->examples[i]->x)->psi->to_string();
    fprintf(fout, "%d %s\n", ((MulticlassStructuredLabel*)d->examples[i]->y)->class_id, data);
    free(data);
  }
  fclose(fout);
  Unlock();

  if(debugLevel > 0 && start_from == 0) fprintf(stderr, "done\n");

  return true;
}


#ifndef NO_SERVER

#include "online_interactive_server.h"

/** 
 * @brief Run the structured learning server
 *
 * To run as a standalone structured learning algorithm with a fixed dataset file, run something like
 *   Train: ./online_interactive_server -p data/params.txt -d data/train.dat -o data/learned_model.txt
 *   Test:  ./online_interactive_server -i data/learned_model.txt -t data/test.dat data/predictions.dat
 *
 * where 
 *   data/params.txt is in the format of StructuredSVM::Load()
 *   data/train.dat is a training set in the format of  StructuredSVM::LoadDataset() (the default implementation 
 *       reads each training example as a line "y x", where y is in the format of StructuredLabel::read()
 *       and x is in the format of StructuredData::read()
 *   data/learned_model.txt is the file where the learned model is written
 *   data/test.dat is a testset in the format of  StructuredSVM::LoadDataset() 
 *   data/predictions.txt is the file where predictions for all labels are written, where each line
 *      corresponds to a test example in the format 
 *          "y_predicted y_ground_truth loss score_prediction score_ground_truth"
 *
 *
 *
 * To run as a server that trains in online fashion, allowing a client to interactively classify examples
 * and add new training examples, run something like
 *   ./online_interactive_server -P 8086 -p data/params.txt -d data/initial_train.txt
 *
 * where data/initial_train.dat is optional and 8086 is the port in which the serve listens on.
 *
 * See StructuredLearnerRpc for info on the network protocol
 *
 **/ 
int main(int argc, const char **argv) {
  StructuredLearnerRpc v(new MulticlassStructuredSVM);
  v.main(argc, argv);
}

#endif
