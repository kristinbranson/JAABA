#include "structured_svm.h"

#ifdef WIN32
//#include <windows.h>
//void Sleep(unsigned long  dwMilliseconds);
void usleep(int us) { /*Sleep(us*1000);*/ }
#else
#include <unistd.h>
#endif

#ifdef HAVE_SVM_STRUCT
extern StructuredSVM *g_learner;
int main_train (int argc, char* argv[]);
#endif


bool ReadString(char *str, FILE *fin) {
  int len;
  return fread(&len, sizeof(int), 1, fin) && fread(str, sizeof(char), len, fin);
}

bool WriteString(char *str, FILE *fout) {
  int len=strlen(str)+1;
  return fwrite(&len, sizeof(int), 1, fout) && fwrite(str, sizeof(char), len, fout);
}

void StructuredSVM::Train(const char *modelout, bool saveFull, const char *initial_sample_set)  {
  runForever = saveFull && method != SPO_MAP_TO_BINARY_MINE_HARD_NEGATIVES;
  if(!trainset) trainset = new StructuredDataset;

  finished = false;
  lambda = 1.0/C;

  // If useFixedSampleSet=true, never call Inference() or ImportanceSample() to concentrate on 
  // different samples.  Instead, use a fixed set of of feature vectors
  useFixedSampleSet = method == SPO_MAP_TO_BINARY || 
                      method == SPO_MAP_TO_BINARY_MINE_HARD_NEGATIVES || 
                      method == SPO_FIXED_SAMPLE_SET;

  // If cache_old_examples=true, then cache features for extracted samples, such that we can
  // optimize with respect to extracted samples in future iterations before calling 
  // Inference() or ImportanceSample() again
  cache_old_examples = method == SPO_DUAL_UPDATE_WITH_CACHE ||
                       method == SPO_DUAL_MULTI_SAMPLE_UPDATE_WITH_CACHE || 
                       useFixedSampleSet;

  // If isMultiSample=true, then call ImportanceSample() instead of Inference(), such that
  // we extract a set of samples in each iteration instead of just one
  isMultiSample = method == SPO_DUAL_MULTI_SAMPLE_UPDATE ||
                  method == SPO_DUAL_MULTI_SAMPLE_UPDATE_WITH_CACHE ||
                  method == SPO_FIXED_SAMPLE_SET;


  // Run SVM^struct instead of our own optimizer
  if(method == SPO_CUTTING_PLANE || method == SPO_CUTTING_PLANE_SMALL_BATCH_SIZE) {
#ifdef HAVE_SVM_STRUCT
    g_learner = this;
    char C_str[1000], eps_str[1000], modelfile[1000];
    sprintf(C_str, "%lf", C);
    sprintf(eps_str, "%lf", eps);
    if(modelout) strcpy(modelfile, modelout);
    char* argv[11] = {(char*)"./svm_struct_learn.out", (char*)"-c", C_str, (char*)"-e", eps_str, (char*)"-M", modelfile, (char*)"-b", (char*)"1", (char*)"-w", (char*)"4" };
    start_time = get_runtime();
    main_train(method == SPO_CUTTING_PLANE_SMALL_BATCH_SIZE ? 11 : (modelout ? 7 : 5), argv);
#else
    fprintf(stderr, "ERROR: to train using SVM^struct, you must define -DHAVE_SVM_STRUCT in the Makefile and link the SVM^struct library\n");
#endif
    return;
  }

  // Begin with a pre-extracted set of samples (e.g., cached features from a previous training session)
  if(initial_sample_set) {
    LoadCachedExamples(initial_sample_set);
    SetSumWScale(lambda*n);
  }

  start_time = get_runtime();
  if(method == SPO_MAP_TO_BINARY || method == SPO_MAP_TO_BINARY_MINE_HARD_NEGATIVES) {
    // Convert structured training set to a binary training set by randomly extracting samples
    // or by mining hard negatives
    TrainBinary(modelout, saveFull, initial_sample_set);
  } else {
    // Do structured learning
    TrainMain(modelout, saveFull, initial_sample_set);
  }

  fprintf(stderr, "Training finished in %lf seconds\n", GetElapsedTime()); 

  if(modelout)
    Save(modelout, saveFull);
}


void StructuredSVM::TrainBinary(const char *modelout, bool saveFull, const char *initial_sample_set)  {
  int num_iter = method == SPO_MAP_TO_BINARY_MINE_HARD_NEGATIVES ? numMineHardNegativesRound : 1;
 
  if(debugLevel > 2) debugLevel--;

  // For some methods (mostly those based on binary classification), pre-extract a training set of samples 
  for(int iter = 0; iter < num_iter; iter++) {
    long old_n = n;
    StructuredDataset *old_trainset = trainset;
    char sample_set_name[1000];
    if(method == SPO_MAP_TO_BINARY_MINE_HARD_NEGATIVES) {
      sprintf(sample_set_name, "%s.samples.%d", modelout, iter);
      maxIters = my_min(maxIters,1000000);
    } else
      sprintf(sample_set_name, "%s.samples", modelout);

    // Extract a binary sample set, or load a cached sample set from disk if it exists
    double elapsedTime = GetElapsedTime();
    if(!FileExists(sample_set_name)) {
      ExtractSampleSet(numHardNegativesPerExample, iter > 0);
      SaveCachedExamples(sample_set_name);
    } else if(method == SPO_MAP_TO_BINARY_MINE_HARD_NEGATIVES) {
      int i = iter;
      while(FileExists(sample_set_name)) 
	sprintf(sample_set_name, "%s.samples.%d", modelout, ++i);
      iter = mineHardNegativesRound = i-1;
      sprintf(sample_set_name, "%s.samples.%d", modelout, mineHardNegativesRound);
      fprintf(stderr, "Loading %s...\n", sample_set_name);
      LoadCachedExamples(sample_set_name);
    }
    double elapsedTime2 = GetElapsedTime();
    if(mineHardNegativesRound) 
      start_time += elapsedTime2-elapsedTime;
    old_n = n;
    old_trainset = trainset;
    if(method == SPO_MAP_TO_BINARY || method == SPO_MAP_TO_BINARY_MINE_HARD_NEGATIVES)
      ConvertCachedExamplesToBinaryTrainingSet();

    // Train a binary classifier
    TrainMain(modelout, saveFull, initial_sample_set);
    if(modelout)
      DumpModelIfNecessary(modelout, true);

    // Restore the original structured training set
    for(int i = 0; i < n; i++) {
      if(trainset->examples[i]->set) {
	trainset->examples[i]->set->samples = NULL;
	trainset->examples[i]->set->num_samples = 0;
	trainset->examples[i]->set->psi_gt = NULL;
	free_SVM_cached_sample_set(trainset->examples[i]->set);
	trainset->examples[i]->set = NULL;
      }
    }
    trainset = old_trainset;
    n = old_n;

    mineHardNegativesRound++;
  }
}


void StructuredSVM::TrainMain(const char *modelout, bool saveFull, const char *initial_sample_set) {
  int num = 0;

  int numThreads = num_thr;
  if(!runMultiThreaded) numThreads = 1;
  else if(runMultiThreaded > 1) numThreads = runMultiThreaded;
  //if(numThreads < 2 && cache_old_examples) numThreads = 2;

  SetTrainset(trainset);
  
  if(!sum_w) {
    sum_w = new SparseVector;
    sum_w->make_non_sparse(true, sizePsi);
  }

  if(cache_old_examples) {
    n = trainset->num_examples;
    SetSumWScale(lambda*n);
  }
  //window = n;

  if(!examples_by_iteration_number) {
    examples_by_iteration_number = (int*)realloc(examples_by_iteration_number, sizeof(int)*(M+2));
    for(int i = 0; i <= M+1; i++) examples_by_iteration_number[i] = -1;
  }

  /* some training information */
  if(debugLevel > 0) {
    char mstr[1000];  OptimizationMethodToString(method, mstr);
    printf("Number of threads=%d\n", numThreads);
    printf("Optimization method=%s\n", mstr);
    printf("Regularization constant (lambda): %.8lg\n", lambda);
    printf("Approximation factor (epsilon): %.8lg\n", eps);
    printf("Number of training examples (n): %ld\n", (long)trainset->num_examples);
    printf("Feature space dimension (sizePsi): %d\n", sizePsi); fflush(stdout);
  }

  #pragma omp parallel num_threads(numThreads)
  {
    int tid = omp_get_thread_num();
    if(tid > 0 || !cache_old_examples || numThreads<=1 || !updateFromCacheThread) {
      // Worker threads, continuously call ybar=find_most_violated_constraint and then use ybar to update the current weights
      int i = -1;
      while(!finished) {
        // Choose a training example 'i' to process
        Lock();
        double w_scale = numThreads==1 && canScaleW ? (t ? 1.0/sum_w_scale : 1) : 1;
        i = ChooseNextExample();
        if(i < 0) {
          Unlock();
          usleep(100000);
          continue;
        }
        StructuredExample *ex = trainset->examples[i];

        SVM_cached_sample_set *set = trainset->examples[i]->set;
        if(!set) {
          set = new_SVM_cached_sample_set(i, Psi(ex->x, ex->y).ptr());
          set->psi_gt_sqr = set->psi_gt->dot(*set->psi_gt);
        } 
        set->lock = true;
        if(cache_old_examples && !useFixedSampleSet && numCacheUpdatesPerIteration && t) {
          // Optimize weights with respect to previously cached samples for this example, to ensure returned samples are as independent as possible
          if(ex_num_iters[i])  
            UpdateFromCache(false, &num, i);
          if(!updateFromCacheThread)   // Optimize weights with respect to cached samples that are randomly selected
            for(int k = 0; k < numCacheUpdatesPerIteration; k++)
              UpdateFromCache(false, &num, rand()%my_min(t,n));
        }
        SparseVector *w = numThreads==1 && canScaleW ? sum_w : GetCurrentWeights(false);
        set->score_gt = w->dot(*set->psi_gt)*w_scale;   // <w_t,psi(x_i,y_i)>
        Unlock();

        VFLOAT score = 0;
        double score_loss = 0;

        if(!useFixedSampleSet) {
          if(!isMultiSample) {
            // Find the most violated label ybar = max_y <w_t,psi(x_i,y)>+loss(y_i,y)
            StructuredLabel *ybar = NewStructuredLabel(ex->x);
            score_loss = Inference(ex->x, ybar, w, NULL, ex->y, w_scale);
            SVM_cached_sample *sample = SVM_cached_sample_set_add_sample(set, ybar);
            sample->slack = score_loss - set->score_gt;
          } else {
	    // Extract a set of samples with non-zero slack.  Should include the most violated label ybar
            score_loss = ImportanceSample(ex->x, w, ex->y, set, w_scale);
          }

          SVM_cached_sample_set_compute_features(set, ex);
          set->loss = set->num_samples ? set->samples[0].loss : 0;
          if(set->ybar) delete set->ybar;
          set->ybar = NULL;
          if(set->num_samples) {
            set->ybar = NewStructuredLabel(ex->x);
            Json::Value yy = set->samples[0].ybar->save(this);
            set->ybar->load(yy, this);
          }
          score = score_loss-set->loss;      // <w_t,psi(x_i,y)>
          set->slack_before = score_loss - set->score_gt;

          OnFinishedIteration(ex->x, ex->y);  // called because the API user might want to free temporary memory caches
        } else {
	  // Optimize over pre-extracted sample set, instead of dynamically choosing new samples using Inference() or ImportanceSample()
          set->slack_before = -1000000;
          for(int j = 0; j < set->num_samples; j++) {
            set->samples[j].dot_w = sum_w->dot(*set->samples[j].psi) - set->score_gt*sum_w_scale;
            set->samples[j].slack = set->samples[j].dot_w/sum_w_scale + set->samples[j].loss;
            if(set->samples[j].slack > set->slack_before) {
              set->slack_before = set->samples[j].slack;
              set->loss = set->samples[j].loss;
              score_loss = set->slack_before + set->score_gt;
              score = score_loss-set->loss;
            }
          }
        }
        set->numIters++;
        set->sumSlack += set->slack_before;

        Lock();
        UpdateWeights(set, -1);       // Use the sample(s) for training example i to update the model weights
        UpdateStatistics(set, t-1);   // Update estimates of training and model error
        UpdateExampleIterationQueue(set, t-1);  // Maintain a queue used to select which example to process next
        CheckConvergence();           // Check if optimization is finished

        if(cache_old_examples) 
	  trainset->examples[i]->set = set;

        if(debugLevel > 2) {
          printf("Example %d: m=%d slack=%f->%f score=%f score_gt=%f loss=%f alpha=%f num_samples=%d\n",
                 i, ex_num_iters[i], (VFLOAT)set->slack_before, (VFLOAT)set->slack_after, (VFLOAT)(score), 
                 (VFLOAT)(set->score_gt), (VFLOAT)set->loss, (VFLOAT)set->alpha, set->num_samples);
        }
        set->lock = false;

        // Cleanup
        if(numThreads!=1 || !canScaleW) delete w;
        if(!cache_old_examples) 
          free_SVM_cached_sample_set(set);
        else if(!useFixedSampleSet && mergeSamples) {
          CondenseSamples(set);
        }

	if(t % n == n-1) {
	  OnFinishedPassThroughTrainset();   
	  if(trainLatent) {
	    // If training with respect to latent variables, infer new latent variables every time we pass through the
	    // dataset once
	    InferLatentValues(trainset); 
	  }  
	}

        Unlock();

        if(modelout)
          DumpModelIfNecessary(modelout);
      }
    } else  {
      // If statement occurs only if updateFromCacheThread==true
      // Optimization thread, each iteration randomly selects a previous cached example ybar_i (earlier calls
      // to finding the most violated constraint) and then makes an update to w with respect to that sample.  Can be useful
      // if the operation of finding the most violated constraint is much more expensive than the operation of updating w
      while(!finished) {
        //usleep(1);
        while(t <= numThreads) // wait for at least one first constraint
          usleep(100000);
	UpdateFromCache(true, &num);
      }
    }
  }
}


void StructuredSVM::UpdateFromCache(bool lock, int *num, int i) {
  // Choose a label from a random iteration (a label ybar from a previous worker thread iteration), and optimize its
  // dual parameters
  if(!n) return;
  if(lock) Lock();
  if(i < 0) i = rand()%(n);
  if(!trainset->examples[i]->set || trainset->examples[i]->set->lock) {
    if(lock) Unlock();
    return;
  }
  trainset->examples[i]->set->score_gt = sum_w->dot(*trainset->examples[i]->set->psi_gt)/sum_w_scale;   // <w_t,psi(x_i,y_i)>
  for(int j = 0; j < trainset->examples[i]->set->num_samples; j++) 
    trainset->examples[i]->set->samples[j].dot_w = sum_w ? sum_w->dot(*trainset->examples[i]->set->samples[j].psi) - 
      trainset->examples[i]->set->score_gt*sum_w_scale : NULL;
  UpdateWeights(trainset->examples[i]->set, i);

  if(num) {
    *num++;
    if(*num >= t*50) {
      // Sanity check: occasionally recompute weights from dual parameters.  Could help avoid drifting due to numerical precision errors
      RecomputeWeights(true);
      *num = 0;
    }
  }
  if(lock) Unlock();
}


/*
 * The rules for choosing which example to process next are (from highest to lowest precedence):
 *   1) Never allow different threads to process the same example at the same time
 *   2) Prefer processing an example that has been iterated over at least once but less than minItersBeforeNewExample
 *      iterations before going onto a new example (allows better estimation of generalization and optimization errors)
 *   3) Otherwise, prefer processing an example that has been iterated over the fewest times
 */
int StructuredSVM::ChooseNextExample() {
  if(!trainset->num_examples) return -1;

  // Choose the examples in order, selecting the one that has been processed the fewest number of times
  int it = currMinIterByExample;
  while(it <= M && examples_by_iteration_number[it] < 0)
    it++;
  if(hasConverged && it == M)
    return -1;
  if(it > M || examples_by_iteration_number[it] < 0)
    return -1;

  if(it < minItersBeforeNewExample) {
    // For better estimation of generalization error, we can process each example minItersBeforeNewExample times before moving onto the next example
    int it2 = 1;
    while(it2 <= M && it2 < minItersBeforeNewExample && examples_by_iteration_number[it2] < 0)
      it2++;
    if(it2 < minItersBeforeNewExample && examples_by_iteration_number[it2] >= 0 && it2 <= M)
      it = it2;
  }

  // Remove (temporarily) the selected example, such that no other threads can process it at the same time
  int retval = examples_by_iteration_number[it];
  if(examples_by_iteration_next_ind[retval] == retval) {
    examples_by_iteration_number[it] = -1;
  } else {
    examples_by_iteration_number[it] = examples_by_iteration_next_ind[retval];
    examples_by_iteration_prev_ind[examples_by_iteration_next_ind[retval]] = examples_by_iteration_prev_ind[retval];
    examples_by_iteration_next_ind[examples_by_iteration_prev_ind[retval]] = examples_by_iteration_next_ind[retval];
  }
  examples_by_iteration_next_ind[retval] = examples_by_iteration_prev_ind[retval] = -1;

  return retval;
}


VFLOAT StructuredSVM::Test(const char *testfile, const char *predictionsFile, const char *htmlDir, double *svm_err) {
  StructuredDataset *testset = LoadDataset(testfile);
  VFLOAT v = Test(testset, predictionsFile, htmlDir, svm_err, true);
  delete testset;
  return v;
}


void StructuredSVM::VisualizeDataset(StructuredDataset *dataset, const char *htmlDir, int max_examples) {
  Lock();
  CreateDirectoryIfNecessary(htmlDir);
  int nc = NumHTMLColumns();
  char fname[1000];  sprintf(fname, "%s/index.html", htmlDir);
  FILE *fout = fopen(fname, "w");
  if(!fout) { fprintf(stderr, "Could not open %s for writing\n", fname); return; }
  fprintf(fout, "<html><table>\n");
  for(int i = 0; i < dataset->num_examples && (max_examples < 0 || i < max_examples); i++) {
    char *htmlStr = VisualizeExample(htmlDir, dataset->examples[i]);
    if(i%nc == 0) {
      if(i) fprintf(fout, "</tr>\n");
      fprintf(fout, "<tr>\n");
    }
    fprintf(fout, "<td>%s</td>\n", htmlStr);
    free(htmlStr);
  }
  fprintf(fout, "</tr></table></html>\n");
  fclose(fout);
  Unlock();
}

VFLOAT StructuredSVM::Test(StructuredDataset *testset, const char *predictionsFile, const char *htmlDir, double *svm_err, bool getLock) {
  SparseVector *w = GetCurrentWeights(getLock);
  double sum_slack = 0;

  if(getLock) Lock();
  int nc = NumHTMLColumns();
  double sum_los = 0;
  char **strs;
  int num = 0;
  omp_lock_t l_lock;
  omp_init_lock(&l_lock);

  FILE *htmlOut = NULL;
  if(htmlDir) {
    char fname[1000];  sprintf(fname, "%s/index.html", htmlDir);
    htmlOut = fopen(fname, "w");
    fprintf(htmlOut, "<html><table>\n");
  }

  if(debugLevel > 0) fprintf(stderr, "Evaluating testset...\n");
  if(predictionsFile)
    strs = (char**)malloc(sizeof(char*)*testset->num_examples);

  // Getting strange problems in Visual C++ Release mode.  num and sum_los
  // aren't getting shared in the main loop below unless I run this loop first
  #pragma omp parallel 
  for(int i = 0; i < 100; i++) {
    num++;
    sum_los++;
    sum_slack++;
  }
  num = 0;
  sum_los = sum_slack = 0;

#pragma omp parallel for
  for(int i = 0; i < testset->num_examples; i++) {
    StructuredLabel *y = NewStructuredLabel(testset->examples[i]->x);
    double score = Inference(testset->examples[i]->x, y, w);
    double los = Loss(testset->examples[i]->y, y);
    double score_gt = w->dot(Psi(testset->examples[i]->x, testset->examples[i]->y)); 
    double slack = score-score_gt+los;

    if(predictionsFile) {
      Json::Value o, pred, gt;
      pred["y"] = y->save(this);
      pred["score"] = score;
      gt["y"] = testset->examples[i]->y->save(this);
      gt["score"] = score_gt;
      o["predicted"] = y->save(this);
      o["ground_truth"] = gt;
      o["loss"] = los;
      Json::FastWriter writer;
      char tmp[100000];
      strcpy(tmp, writer.write(o).c_str());
      //fprintf(stderr, "%s\n", tmp);
      strs[i] = StringCopy(tmp);
    }
    delete y;

    omp_set_lock(&l_lock);
    sum_los += los;
    sum_slack += slack;
    num++;
    if(htmlOut) {
      char *htmlStr = VisualizeExample(htmlDir, testset->examples[i]); 
      if(i%nc == 0) {
	if(i) fprintf(htmlOut, "</tr>\n");
	fprintf(htmlOut, "<tr>\n");
      }
      fprintf(htmlOut, "<td>%s</td>\n", htmlStr);
      free(htmlStr);
    }
    fprintf(stderr, "After %d examples: ave_loss=%f, ave_slack=%f\n",  num, (float)(sum_los/num), (float)(sum_slack/num));
    omp_unset_lock(&l_lock);

    OnFinishedIteration(testset->examples[i]->x, testset->examples[i]->y);
  }
  if(htmlOut) {
    fprintf(htmlOut, "</tr></table></html>\n");
    fclose(htmlOut);
  }
  double svm_error = (float)(sum_slack/testset->num_examples) + w->dot(*w)*lambda/2;
  printf("Average loss was %f, slack=%f, svm_err=%f\n", (float)(sum_los/testset->num_examples), 
	 (float)(sum_slack/testset->num_examples), svm_error);

  omp_destroy_lock(&l_lock);
  if(predictionsFile) {
    FILE *fout = fopen(predictionsFile, "w");
    for(int i = 0; i < testset->num_examples; i++) {
      fprintf(fout, "%s", strs[i]);
      free(strs[i]);
    }
    fclose(fout);
    free(strs);
  }

  if(getLock) Unlock();

  if(svm_err) *svm_err = svm_error;

  return (sum_los/testset->num_examples);
}

void StructuredSVM::SetSumWScale(double sum_w_scale_new) {
  sum_dual += .5*((sum_w_scale ? 1/sum_w_scale : 0) - 1/sum_w_scale_new)*sum_w_sqr;
  sum_w_scale = sum_w_scale_new;
  regularization_error = (sum_w_sqr/SQR(sum_w_scale))*lambda/2;
}

void StructuredSVM::UpdateWeights(SVM_cached_sample_set *ex, int iterInd) {
  bool bound_w = method != SPO_SGD; 
  
  for(int j = 0; j < ex->num_samples; j++)
    if(ex->samples[j].loss > maxLoss) 
      maxLoss = ex->samples[j].loss;

  if(ex_num_iters[ex->i] == 0) {
    ex_first_iter[ex->i] = t;
    last_example = ex->i;
    if(n < trainset->num_examples) {
      n++;
      if(cache_old_examples) 
        SetSumWScale(lambda*n);
    }
  }

  if(iterInd == -1) {// When we increase t, w^{t+1} = t/(t+1) * w^{t}, causing a change in the regularization error and dual
    t++;
    if(!cache_old_examples) 
      SetSumWScale(lambda*t);
  }

 
  if(!isMultiSample && !ex->u_i) {
    // Update the model by taking a step in the direction of the subgradient v=ex->psi_gt-ex->sample[0].psi
    SingleSampleUpdate(ex, method != SPO_SGD && method != SPO_SGD_PEGASOS);
  } else {
    // Take a more complicated update step with respect to multiple samples, instead of just the
    // sub-gradient
    MultiSampleUpdate(ex, trainset->examples[ex->i], iterInd >= 0 ? 1 : numMultiSampleIterations);
  }

  if(bound_w) {
    // Project sum_w onto the L2 ball, such that ||w||^2 <= 1/lambda
    // This is the projection step used by a Pegasos-like update
    if(regularization_error > .5*maxLoss) {
      // scale w by s = (1/sqrt(lambda))/|w| = sqrt(1/lambda/w^2) = sqrt(1/(2*regularization_error))
      double s = sqrt(maxLoss / (2*regularization_error));
      *sum_w *= s;
      sum_w_sqr = 1.0/lambda*SQR(sum_w_scale);
      sum_alpha_loss *= s;
      sum_dual = -sum_w_sqr/(2*sum_w_scale) + sum_alpha_loss;
      //sum_dual = s*sum_dual + (1-s)*t/(2*s);
      regularization_error = .5*maxLoss;
    }
  }

  if(t > 1000 && ex->i % n == n-1)  // avoid numerical drifting problems
    RecomputeWeights(false);
}


// Book-keeping stuff, for estimating generalization error, optimization error, and regret when
// a new example ex->i is processed in iteration t=iter
void StructuredSVM::UpdateStatistics(SVM_cached_sample_set *ex, int iter) {
  double e = ex->slack_before;  // the slack at iteration iter 
  double elapsedTime = GetElapsedTime();

  if(ex_num_iters[ex->i] == 0) {
    // If this is the first time ex was processed, update estimate of the generalization error
    sum_generalization_error += my_max(e,0); // Stores the sum of the (online estimate of) test error
    generalization_errors_by_n[ex->i] = my_max(e,0)+regularization_error;  
    if(ex->i >= window) sum_generalization_error_window -= generalization_errors_by_n[ex->i-window];
    sum_generalization_error_window += generalization_errors_by_n[ex->i];
  }

  // Allocate memory buffers
  if(t+1 > alloc_t) {
    alloc_t = (int)(alloc_t*1.1)+10;
    iter_examples = (long*)realloc(iter_examples, sizeof(long)*alloc_t);
    iter_errors_by_t = (double*)realloc(iter_errors_by_t, sizeof(double)*alloc_t);
    sum_dual_by_t = (double*)realloc(sum_dual_by_t, sizeof(double)*alloc_t);
    generalization_errors_by_t = (double*)realloc(generalization_errors_by_t, sizeof(double)*alloc_t);
    regularization_errors_by_t = (double*)realloc(regularization_errors_by_t, sizeof(double)*alloc_t);
    losses_by_t = (double*)realloc(losses_by_t, sizeof(double)*alloc_t);
    elapsed_time_by_t = (double*)realloc(elapsed_time_by_t, sizeof(double)*alloc_t);
  }

  sum_iter_error += my_max(e,0);   // Stores the sum of the (online estimate of) training error
  iter_errors_by_t[iter] = my_max(e,0);  // Stores the slack (training error) in each iteration
  sum_dual_by_t[iter] = sum_dual;  // Stores the value of the dual objective, which lowerbounds the minimum achievable model error
  generalization_errors_by_t[iter] = ex_num_iters[ex->i] == 0 ? generalization_errors_by_n[ex->i] : sum_generalization_error_window/my_min(n,window);  
  regularization_errors_by_t[iter] = regularization_error;
  losses_by_t[iter] = ex->num_samples ? ex->samples[0].loss : 0;  // Stores the loss (not including the slack violation)
  elapsed_time_by_t[iter] = (double)GetElapsedTime();

  // Error measured over last set of examples/iterations of size window
  int curr_window_t = my_min(t,window), curr_window_n = my_min(n,window);
  if(iter >= window) 
    sum_iter_error_window += iter_errors_by_t[iter] - iter_errors_by_t[iter-window];
  else
    sum_iter_error_window += iter_errors_by_t[iter];

  double nn = (double)(cache_old_examples ? n : t);
  if(debugLevel > 2 || (debugLevel > 1 && iter%10000==9999))
    printf("t=%d, n=%d, time=%f: Average Training Error=%f (Model error=%f, Optimization error=%f, Regularization error=%f), Out of sample error=%f\n",
	   (int)t, (int)n, (float)elapsedTime, (VFLOAT)(sum_iter_error/t)+(VFLOAT)regularization_error, (VFLOAT)(sum_dual/nn), (VFLOAT)((sum_iter_error)/t-sum_dual/nn+regularization_error), (VFLOAT)regularization_error, (VFLOAT)(sum_generalization_error/n));


    
  if(debugLevel > 2 || (debugLevel > 1 && iter%10000==9999))
    printf("Last %d iters: Average Training Error=%f (Model error=%f, Optimization error=%f, Regularization error=%f), Out of sample error=%f\n",
	   (int)curr_window_t, (float)(sum_iter_error_window/curr_window_t)+(float)regularization_error, (float)(sum_dual/nn), (float)((sum_iter_error_window)/curr_window_t)-sum_dual/nn+regularization_error, (float)regularization_error, (float)(sum_generalization_error_window/curr_window_n));
}


void StructuredSVM::DumpModelIfNecessary(const char *modelfile, bool force) {
  double elapsedTime = GetElapsedTime();
  if((dumpModelStartTime && elapsedTime >= dumpModelStartTime*pow(dumpModelFactor,numModelDumps)) || force) {
    double tm_beg = get_runtime();
    char modelfile_times[1000], modelfile_dump[1000];
    sprintf(modelfile_times, "%s.times", modelfile);
    sprintf(modelfile_dump, "%s.%d", modelfile, numModelDumps);
    Save(modelfile_dump, false);
    double loss = 0, svm_err = 0;
    if(validationfile)
      loss = Test(validationfile, NULL, NULL, &svm_err);

    FILE *fout = fopen(modelfile_times, numModelDumps ? "a" : "w");
    fprintf(fout, "%d %lf %d %lf %lf %lf\n", numModelDumps, elapsedTime, (int)t, sum_dual/(cache_old_examples ? n : t), loss, svm_err);
    fclose(fout);
    while(elapsedTime >= dumpModelStartTime*pow(dumpModelFactor,numModelDumps)) 
      numModelDumps++;
    double tm_end = get_runtime();
    start_time += tm_end-tm_beg;
  }
}

void StructuredSVM::CheckConvergence() {
  int curr_window_t = my_min(t,window), curr_window_n = my_min(n,window);
  double nn = (double)(cache_old_examples ? n : t);
  double eps_empirical_measured = (sum_iter_error_window) / curr_window_t - sum_dual/nn + regularization_error;
  double eps_generalization_measured = sum_generalization_error_window / curr_window_n - sum_dual/nn;

  if(!hasConverged && ((eps && eps_empirical_measured < eps && !finished && (!runForever || !cache_old_examples)) || 
		       t > maxIters)) {
    if(t > window) {
      if(!runForever)
	finished = true;
      if(debugLevel > 0 && !hasConverged) {
	printf("%s at t=%d: epsilon_measured=%f\n", runForever ? "Convergence of empirical error detected" : "Finishing", (int)t, (float)eps_empirical_measured);
	printf("Last %d iters: Average Training Error=%f (Model error=%f, Optimization error=%f, Regularization error=%f), Out of sample error=%f\n",
	       (int)curr_window_t, (float)(sum_iter_error_window/curr_window_t)+(float)regularization_error, (float)(sum_dual/nn), (float)((sum_iter_error_window)/curr_window_t)-sum_dual/nn+regularization_error, (float)regularization_error, (float)(sum_generalization_error_window/curr_window_n));
      }
      if(!runForever && M) {
	finished = true;
      }
    }
  } else if(0 && !finished && n > window && eps_generalization_measured < eps && !runForever) {
    printf("Convergence of generalization error detected at t=%d: epsilon_measured=%f\n", (int)n, (float)eps_generalization_measured);
    finished = true;
  }
} 

void StructuredSVM::UpdateExampleIterationQueue(SVM_cached_sample_set *ex, int iter) {
  iter_examples[iter] = ex->i;  
  ex_num_iters[ex->i]++;

  if(ex_num_iters[ex->i] > M) {
    assert(M == ex_num_iters[ex->i]-1);
    M++;
    double nn = (double)(cache_old_examples ? n : t);
    printf("M=%d, t=%d, n=%d, dual=%f\n",
	   (int)M-1, (int)t, (int)n, (float)(sum_dual/nn));
    
    examples_by_iteration_number = (int*)realloc(examples_by_iteration_number, sizeof(int)*(M+2));
    examples_by_iteration_number[M] = -1;
  }
  
  // Add this example back to the queue of examples to be processed
  int r = examples_by_iteration_number[ex_num_iters[ex->i]];
  if(r >= 0) {
    examples_by_iteration_prev_ind[ex->i] = examples_by_iteration_prev_ind[r];
    examples_by_iteration_next_ind[ex->i] = r;
    examples_by_iteration_next_ind[examples_by_iteration_prev_ind[r]] = ex->i;
    examples_by_iteration_prev_ind[r] = ex->i;
  } else {
    examples_by_iteration_prev_ind[ex->i] = examples_by_iteration_next_ind[ex->i] = ex->i;
    examples_by_iteration_number[ex_num_iters[ex->i]] = ex->i;
  }
}


StructuredExample *StructuredSVM::CopyExample(StructuredData *x, StructuredLabel *y, StructuredLabel *y_latent) {
  StructuredExample *copy = new StructuredExample();
  copy->x = NewStructuredData();
  copy->y = NewStructuredLabel(copy->x);
  Json::Value xx = x->save(this);
  bool b = copy->x->load(xx, this);
  assert(b);
  Json::Value yy = y->save(this);
  b = copy->y->load(yy, this);
  assert(b);
  if(y_latent) {
    Json::Value yy_latent = y_latent->save(this);
    b = copy->y_latent->load(yy_latent, this);
    assert(b);
  }

  return copy;
}
int StructuredSVM::AddExample(StructuredData *x, StructuredLabel *y) {
  int retval = -1;
  Lock();
  hasConverged = false;

  // Add copy of ex to list of examples
  assert(trainset);
  retval = trainset->num_examples;
  trainset->AddExample(CopyExample(x,y));
  CreateTrainingExampleQueues(retval);
  Unlock();

  return retval;
}


void StructuredSVM::CreateTrainingExampleQueues(int ind) {
  if(trainset->num_examples+1 > alloc_n) {
    alloc_n = (int)(alloc_n*1.1 + 10);
    ex_num_iters = (int*)realloc(ex_num_iters, sizeof(int)*alloc_n);
    ex_first_iter = (int*)realloc(ex_first_iter, sizeof(int)*alloc_n);
    generalization_errors_by_n = (double*)realloc(generalization_errors_by_n, sizeof(double)*alloc_n);
    examples_by_iteration_next_ind = (int*)realloc(examples_by_iteration_next_ind, sizeof(int)*alloc_n);
    examples_by_iteration_prev_ind = (int*)realloc(examples_by_iteration_prev_ind, sizeof(int)*alloc_n);
  }
  ex_num_iters[ind] = 0;
  ex_first_iter[ind] = -1;
  generalization_errors_by_n[ind] = 0;
  

  // Add this example to the front of the list of examples that have been iterated over 0 times
  if(examples_by_iteration_number[0] >= 0) {
    examples_by_iteration_prev_ind[ind] = examples_by_iteration_prev_ind[examples_by_iteration_number[0]];
    examples_by_iteration_next_ind[ind] = examples_by_iteration_number[0];
    examples_by_iteration_next_ind[examples_by_iteration_prev_ind[examples_by_iteration_number[0]]] = ind;
    examples_by_iteration_prev_ind[examples_by_iteration_number[0]] = ind;
  } else {
    examples_by_iteration_prev_ind[ind] = examples_by_iteration_next_ind[ind] = ind;
    examples_by_iteration_number[0] = ind;
  }
  this->currMinIterByExample = 0;
}

bool StructuredSVM::SaveTrainingSet(int start_from) {
  if(!trainfile && modelfile) {
    char tmp[1000];
    sprintf(tmp, "%s.train", modelfile);
    trainfile = StringCopy(tmp);
  }
  if(trainfile) return SaveDataset(trainset, trainfile, start_from);

  return false;
}

// Sanity check: recompute weights from dual parameters.  Could avoid drifting due to numerical precision errors
void StructuredSVM::RecomputeWeights(bool full) {
  int i;
  if(full) {
    if(cache_old_examples || useFixedSampleSet) {
      SparseVector *sum_w_new = new SparseVector;
      sum_w_new->make_non_sparse(true, sizePsi);
      sum_alpha_loss = 0;
      for(i = 0; i < n; i++) {
	if(trainset->examples[i]->set) {
	  if(trainset->examples[i]->set->alpha) 
	    *sum_w_new += (*trainset->examples[i]->set->psi_gt*trainset->examples[i]->set->alpha);
	  if(trainset->examples[i]->set->u_i) {
	    *sum_w_new -= *trainset->examples[i]->set->u_i;
	    sum_alpha_loss += trainset->examples[i]->set->D_i;
	  } else {
	    for(int j = 0; j < trainset->examples[i]->set->num_samples; j++) {
	      if(trainset->examples[i]->set->samples[j].alpha)
		*sum_w_new -= (*trainset->examples[i]->set->samples[j].psi * trainset->examples[i]->set->samples[j].alpha);
	      sum_alpha_loss += trainset->examples[i]->set->samples[j].alpha * trainset->examples[i]->set->samples[j].loss;
	    }
	  }
	}
      }
      //assert(sum_sqr_diff_ss(sum_w, sum_w_new) < .001*sum_w_scale);
      delete sum_w;
      sum_w = sum_w_new;
      sum_iter_error = 0;
      for(i = 0; i < t; i++) 
	sum_iter_error += iter_errors_by_t[i];
    }
  }
  sum_iter_error_window = 0;
  for(i = my_max(0,t-window); i < t; i++) 
    sum_iter_error_window += iter_errors_by_t[i];
  sum_w_sqr = sum_w->dot(*sum_w, regularize);
  regularization_error = sum_w_sqr/SQR(sum_w_scale)*lambda/2;
  sum_dual = -sum_w_sqr/(2*sum_w_scale) + sum_alpha_loss;
}

void StructuredSVM::OptimizeAllConstraints(int num_iter) {
  for(int i = 0; i < num_iter; i++) {
    for(int j = 0; j < n; j++)
      UpdateWeights(trainset->examples[i]->set, j);
    RecomputeWeights();
  }
}
void StructuredSVM::SetLambda(double l, int num_iter) {
  Lock();
  lambda = l;
  C = (VFLOAT)(1.0/lambda);
  if(n) SetSumWScale(cache_old_examples ? lambda*n :lambda*t);
  if(sizePsi && t && num_iter) {
    RecomputeWeights();
    OptimizeAllConstraints(num_iter);
  }
  Unlock();
}

double *ComputeWindowAverageArray(double *a, int n, int window, double *a2=NULL, float s2=0, double *a3=NULL, float s3=0) {
  double cumSum = 0;
  int w = 0, i;
  double *retval = (double*)malloc(sizeof(double)*n);

  for(i = 0; i < n; i++) {
    if(w < window) w++;
    else cumSum -= a[i-window];
    cumSum += a[i];
    retval[i] = cumSum/w;
  }
  if(a2) {
    cumSum = 0;
    w = 0;
    for(i = 0; i < n; i++) {
      if(w < window) w++;
      else cumSum -= a2[i-window];
      cumSum += a2[i];
      retval[i] += cumSum/w*s2;
    }
  }
  if(a3) {
    cumSum = 0;
    w = 0;
    for(i = 0; i < n; i++) {
      if(w < window) w++;
      else cumSum -= a3[i-window];
      cumSum += a3[i];
      retval[i] += cumSum/w*s3;
    }
  }
  return retval;
}



void StructuredSVM::GetStatisticsByExample(int ave, long *nn, double **gen_err_buff, double **emp_err_buff, double **model_err_buff, double **reg_err_buff, double **loss_err_buff) {
  Lock();
  double *gen_errors_by_n = (double*)malloc(sizeof(double)*n*5);
  double *emp_errors_by_n = gen_errors_by_n+n;
  double *mod_errors_by_n = emp_errors_by_n+n;
  double *reg_errors_by_n = mod_errors_by_n+n;
  double *loss_by_n = reg_errors_by_n+n;
  int num = 0;
  double tt = cache_old_examples ? n : t;
  for(int i = 0; i < t; i++) {
    if(ex_first_iter[iter_examples[i]] == i) {
      gen_errors_by_n[num] = generalization_errors_by_t[i];
      emp_errors_by_n[num] = iter_errors_by_t[i]+regularization_errors_by_t[i];
      mod_errors_by_n[num] = sum_dual_by_t[i]/my_min(i+1,tt);
      reg_errors_by_n[num] = regularization_errors_by_t[i];
      loss_by_n[num] = losses_by_t[i]/t;
      num++;
    }
  }
  *nn = num;
  if(gen_err_buff) *gen_err_buff = ComputeWindowAverageArray(gen_errors_by_n, num, ave);
  if(emp_err_buff) *emp_err_buff = ComputeWindowAverageArray(emp_errors_by_n, num, ave); 
  if(model_err_buff) *model_err_buff = ComputeWindowAverageArray(mod_errors_by_n, num, ave);
  if(reg_err_buff) *reg_err_buff = ComputeWindowAverageArray(reg_errors_by_n, num, ave);
  if(loss_err_buff) *loss_err_buff = ComputeWindowAverageArray(loss_by_n, num, ave);
  free(gen_errors_by_n);
  Unlock();
}

void StructuredSVM::GetStatisticsByIteration(int ave, long *tt, double **gen_err_buff, double **emp_err_buff, double **model_err_buff,
					     double **reg_err_buff, double **loss_err_buff, double **time_buff) {
  Lock();
  *tt = t;
  double nn = cache_old_examples ? n : t;
  double *ave_dual_by_t = (double*)malloc(sizeof(double)*t);
  for(int i = 0; i < t; i++) ave_dual_by_t[i] = sum_dual_by_t[i]/my_min(i+1,nn);
  if(gen_err_buff) *gen_err_buff = ComputeWindowAverageArray(generalization_errors_by_t, t, ave);
  if(emp_err_buff) *emp_err_buff = ComputeWindowAverageArray(iter_errors_by_t, t, ave, regularization_errors_by_t, 1); 
  if(model_err_buff) *model_err_buff = ComputeWindowAverageArray(ave_dual_by_t, t, ave); 
  if(reg_err_buff) *reg_err_buff = ComputeWindowAverageArray(regularization_errors_by_t, t, ave);
  if(loss_err_buff) *loss_err_buff = ComputeWindowAverageArray(losses_by_t, t, ave);
  if(time_buff) { *time_buff = (double*)malloc(sizeof(double)*(t+1)); memcpy(*time_buff, elapsed_time_by_t, sizeof(double)*t); }
  free(ave_dual_by_t);
  Unlock();
}

SparseVector *StructuredSVM::GetCurrentWeights(bool lock) {
  if(lock) Lock();
  SparseVector *retval = sum_w->mult_scalar(sum_w_scale ? 1.0/(sum_w_scale) : 0, regularize);
  if(lock) Unlock();
  return retval;
}


void free_SVM_cached_sample(SVM_cached_sample *s) {
  if(s->ybar) delete s->ybar;
  if(s->psi) delete s->psi;
  s->ybar = NULL;
  s->psi = NULL;
}

void read_SVM_cached_sample(SVM_cached_sample *s, FILE *fin, StructuredSVM *svm, StructuredData *x, bool readFull) {
  char str[100000];
  ReadString(str, fin);
  s->ybar = svm->NewStructuredLabel(x);
  Json::Reader reader;
  Json::Value v;
  bool b = reader.parse(str, v);
  assert(b);
  s->ybar->load(v, svm);

  s->psi = NULL;
  clear_SVM_cached_sample(s);
  if(readFull) {
    bool b = (fread(&s->loss, sizeof(VFLOAT), 1, fin) && fread(&s->alpha, sizeof(VFLOAT), 1, fin) && fread(&s->sqr, sizeof(VFLOAT), 1, fin) &&
	      fread(&s->slack, sizeof(VFLOAT), 1, fin) && fread(&s->dot_psi_gt_psi, sizeof(VFLOAT), 1, fin) && 
	      fread(&s->dot_w, sizeof(VFLOAT), 1, fin));
    assert(b);
    s->psi = new SparseVector;
    s->psi->read(fin);
  }
}


void write_SVM_cached_sample(SVM_cached_sample *s, FILE *fout, StructuredSVM *svm, bool saveFull) {
  char str[100000];
  Json::FastWriter writer;
  Json::Value v = s->ybar->save(svm);
  strcpy(str, writer.write(v).c_str());
  WriteString(str, fout);
  
  if(saveFull) {
    bool b = (fwrite(&s->loss, sizeof(VFLOAT), 1, fout) && fwrite(&s->alpha, sizeof(VFLOAT), 1, fout) && fwrite(&s->sqr, sizeof(VFLOAT), 1, fout) &&
	      fwrite(&s->slack, sizeof(VFLOAT), 1, fout) && fwrite(&s->dot_psi_gt_psi, sizeof(VFLOAT), 1, fout) && 
	      fwrite(&s->dot_w, sizeof(VFLOAT), 1, fout));
    s->psi->write(fout);
    assert(b);
  }
}



void clear_SVM_cached_sample_set(SVM_cached_sample_set *s) {
  if(s->u_i)
    delete s->u_i;
  s->u_i = NULL;
  if(s->ybar)
    delete s->ybar;
  s->ybar = NULL;
  s->num_samples = 0;
  s->alpha = 0;
  s->loss = 0;
  s->slack_before = 0;
  s->slack_after = 0;
  s->psi_gt_sqr = 0;
  s->D_i = 0;
  s->dot_u_psi_gt = 0;
  s->u_i_sqr = 0;
  s->drift_bits = 0;
  s->lock = false;
  s->numIters = 0;
  s->sumSlack = 0;
}

SVM_cached_sample_set *new_SVM_cached_sample_set(int i, SparseVector *psi_gt) {
  SVM_cached_sample_set *retval = (SVM_cached_sample_set*)malloc(sizeof(SVM_cached_sample_set));
  retval->samples = NULL;
  retval->u_i = NULL;
  retval->ybar = NULL;
  retval->psi_gt = psi_gt;
  retval->i = i;
  clear_SVM_cached_sample_set(retval);
  return retval;
}


void free_SVM_cached_sample_set(SVM_cached_sample_set *s) {
  for(int i = 0; i < s->num_samples; i++)
    free_SVM_cached_sample(&s->samples[i]);
  if(s->samples) free(s->samples);
  if(s->psi_gt) delete s->psi_gt;
  if(s->u_i) delete s->u_i;
  if(s->ybar) delete s->ybar;
  free(s);
}

SVM_cached_sample_set *read_SVM_cached_sample_set(FILE *fin, StructuredSVM *svm, StructuredData *x, bool readFull) {
  int i;
  bool hasFull;
  bool b = (fread(&hasFull, sizeof(bool), 1, fin)); assert(b);
  b = (fread(&i, sizeof(int), 1, fin)); assert(b);
  SVM_cached_sample_set *s = new_SVM_cached_sample_set(i);
  b = (fread(&s->num_samples, sizeof(int), 1, fin)); assert(b);
  s->samples = (SVM_cached_sample*)realloc(s->samples, sizeof(SVM_cached_sample)*(s->num_samples+1));
  for(int i = 0; i < s->num_samples; i++) {
    read_SVM_cached_sample(&s->samples[i], fin, svm, x, hasFull);
    if(!readFull) 
      clear_SVM_cached_sample(&s->samples[i]);
  }
  s->psi_gt = new SparseVector;  
  s->psi_gt->read(fin);

  if(hasFull) {
    b = (fread(&s->alpha, sizeof(double), 1, fin));  assert(b);
    b = (fread(&s->loss, sizeof(double), 1, fin));  assert(b);
    b = (fread(&s->slack_before, sizeof(double), 1, fin));  assert(b);
    b = (fread(&s->slack_after, sizeof(double), 1, fin));  assert(b);
    b = (fread(&s->psi_gt_sqr, sizeof(double), 1, fin));  assert(b);

    b = (fread(&s->dot_sum_w_psi_gt, sizeof(double), 1, fin));  assert(b);
    b = (fread(&s->D_i, sizeof(VFLOAT), 1, fin));  assert(b);
    b = (fread(&s->dot_u_psi_gt, sizeof(VFLOAT), 1, fin));  assert(b);
    b = (fread(&s->u_i_sqr, sizeof(double), 1, fin));  assert(b);
    b = (fread(&s->drift_bits, sizeof(double), 1, fin));  assert(b);
    b = (fread(&s->sumSlack, sizeof(double), 1, fin));  assert(b);
    b = (fread(&s->numIters, sizeof(int), 1, fin));  assert(b);
    bool has_u_i;
    b = (fread(&has_u_i, sizeof(bool), 1, fin));  assert(b);
    if(has_u_i) {
      s->u_i = new SparseVector;  
      s->u_i->read(fin);
    }
    if(!readFull)
      clear_SVM_cached_sample_set(s);
    s->psi_gt_sqr = s->psi_gt->dot(*s->psi_gt);
  }
  return s;
}

void write_SVM_cached_sample_set(SVM_cached_sample_set *s, FILE *fout, StructuredSVM *svm, bool fullWrite) {
  bool b = (fwrite(&fullWrite, sizeof(bool), 1, fout));  assert(b);
  b = (fwrite(&s->i, sizeof(int), 1, fout)); assert(b);
  b = (fwrite(&s->num_samples, sizeof(int), 1, fout)); assert(b);
  for(int i = 0; i < s->num_samples; i++)
    write_SVM_cached_sample(&s->samples[i], fout, svm, fullWrite);
  s->psi_gt->write(fout);

  if(fullWrite) {
    b = (fwrite(&s->alpha, sizeof(double), 1, fout));  assert(b);
    b = (fwrite(&s->loss, sizeof(double), 1, fout));  assert(b);
    b = (fwrite(&s->slack_before, sizeof(double), 1, fout));  assert(b);
    b = (fwrite(&s->slack_after, sizeof(double), 1, fout));  assert(b);
    b = (fwrite(&s->psi_gt_sqr, sizeof(double), 1, fout));  assert(b);

    b = (fwrite(&s->dot_sum_w_psi_gt, sizeof(double), 1, fout));  assert(b);
    b = (fwrite(&s->D_i, sizeof(VFLOAT), 1, fout));  assert(b);
    b = (fwrite(&s->dot_u_psi_gt, sizeof(VFLOAT), 1, fout));  assert(b);
    b = (fwrite(&s->u_i_sqr, sizeof(double), 1, fout));  assert(b);
    b = (fwrite(&s->drift_bits, sizeof(double), 1, fout));  assert(b);
    b = (fwrite(&s->sumSlack, sizeof(double), 1, fout));  assert(b);
    b = (fwrite(&s->numIters, sizeof(int), 1, fout));  assert(b);
    bool has_u_i = s->u_i != NULL;
    b = (fwrite(&has_u_i, sizeof(bool), 1, fout));  assert(b);
    if(has_u_i) s->u_i->write(fout);
  }
}


void clear_SVM_cached_sample(SVM_cached_sample *s) {
  if(s->psi) delete s->psi;
  s->psi = NULL;
  s->loss = 0;
  s->alpha = 0;
  s->slack = 0;
  s->sqr = 0;
  s->dot_w = 0;
  s->dot_psi_gt_psi = 0;
}

SVM_cached_sample *SVM_cached_sample_set_add_sample(SVM_cached_sample_set *s, StructuredLabel *ybar) {
  s->samples = (SVM_cached_sample*)realloc(s->samples, sizeof(SVM_cached_sample)*(s->num_samples+1));
  SVM_cached_sample *retval = s->samples+s->num_samples;
  retval->ybar = ybar;
  retval->psi = NULL;
  clear_SVM_cached_sample(retval);
  s->num_samples++;

  return retval;
}


void StructuredSVM::SVM_cached_sample_set_compute_features(SVM_cached_sample_set *set, StructuredExample *ex) {
  for(int j = 0; j < set->num_samples; j++) {
    if(!set->samples[j].psi) {
      set->samples[j].psi = Psi(ex->x, set->samples[j].ybar).ptr();
      set->samples[j].loss = Loss(ex->y, set->samples[j].ybar);
      set->samples[j].dot_psi_gt_psi = set->psi_gt->dot(*set->samples[j].psi);
      set->samples[j].sqr = set->psi_gt_sqr - 2*set->samples[j].dot_psi_gt_psi + set->samples[j].psi->dot(*set->samples[j].psi);
      set->samples[j].dot_psi_gt_psi = set->psi_gt->dot(*set->samples[j].psi);
      set->samples[j].dot_w = sum_w ? sum_w->dot(*set->samples[j].psi) - set->score_gt*sum_w_scale : NULL;
      assert(!isnan(set->samples[j].sqr));
      //set->samples[j].slack =
    }
  }
}

void StructuredSVM::SetTrainset(StructuredDataset *t) { 
  M = 0;
  trainset=t; 
  int i;
  examples_by_iteration_number = (int*)realloc(examples_by_iteration_number, sizeof(int)*(M+2));
  for(i = 0; i <= M+1; i++)
    examples_by_iteration_number[i] = -1;
  for(int i = 0; i < t->num_examples; i++) {
    CreateTrainingExampleQueues(i); 
  }
  regularization_error = 0;
  sum_dual = 0;
  sum_alpha_loss = 0;
  sum_w_sqr = 0;
  this->t = 0;
  if(sum_w) delete sum_w;
  sum_w = NULL;
}

void StructuredSVM::SingleSampleUpdate(SVM_cached_sample_set *set, bool useSmartStepSize) {
  assert(set->num_samples == 1);
  if(runMultiThreaded)
    set->samples[0].dot_w = sum_w->dot(*set->samples[0].psi) - set->score_gt*sum_w_scale;

  // SPO_SGD and SPO_SGD_PEGASOS: take a step of size -1/(lambda*t) in the direction of the
  // sub-gradient psi(ybar,x)-psi(y_i,x), which corresponds to setting dalpha=1
  double dalpha = 1;
  SVM_cached_sample *s = &set->samples[0];
  SparseVector dpsi = *s->psi - *set->psi_gt;
  double dot = s->dot_w, d_sum_w_sqr;

  // SPO_DUAL_UPDATE and SPO_DUAL_UPDATE_WITH_CACHE: take a step in the direction of the 
  // sub-gradient psi(ybar,x)-psi(y_i,x), where the chosen step size maximizes the dual objective
  if(useSmartStepSize) 
    dalpha = (dot + s->loss*(sum_w_scale)) / my_max(s->sqr,.0000000001);
    
  // SPO_DUAL_UPDATE and SPO_DUAL_UPDATE_WITH_CACHE: ensure 0 <= alpha <= 1
  dalpha = my_min(1-s->alpha, my_max(-s->alpha,dalpha));

  if(dalpha != 0) {
    // Take a step of size dalpha the direction of the sub-gradient psi(ybar,x)-psi(y_i,x)
    // (lambda*t)*w_t = (lambda*(t-1))*w_{t-1} - dalpha*(psi(ybar,x)-psi(y_i,x))
    *sum_w -= (dpsi * (dalpha));
      
    // Keep track of the change in the dual objective, regularization_error, and w^2
    s->alpha += dalpha;
    d_sum_w_sqr = -2*dalpha*dot + SQR(dalpha)*(s->sqr);
    sum_dual += -d_sum_w_sqr/(2*sum_w_scale) + dalpha*s->loss;
    sum_alpha_loss += dalpha*s->loss;
    sum_w_sqr += d_sum_w_sqr;
    regularization_error = (sum_w_sqr/SQR(sum_w_scale))*lambda/2;
    set->alpha = s->alpha;
    set->loss = s->loss;
    set->slack_after = dot/(sum_w_scale)+s->loss - dalpha*(s->sqr)/(sum_w_scale);
    assert(!isnan(sum_w_sqr) && !isnan(sum_dual));
  }
}



//#define DEBUG_MULTI_SAMPLE_UPDATE


// Iteratively optimize alpha parameters for a set of samples 's'.  Implements 'Multi-Sample Dual Update Step' of writeup
void StructuredSVM::MultiSampleUpdate(SVM_cached_sample_set *set, StructuredExample *ex, int R) {
  VFLOAT dot_w_u = 0;   // dot_w_u = <w,u>
  VFLOAT L_i = 0;       // L_i = <w,-u_i/(lambda*t)> + \sum_ybar alpha_{i,ybar} loss(y_i,ybar)
  VFLOAT s_u = 1;
  int j, r;

  //if(!set->num_samples && !set->alpha) return;

  if(!set->u_i) {
    if(!u_i_buff) {
      u_i_buff = (double*)malloc(sizeof(double)*sizePsi);
      for(int i = 0; i < sizePsi; i++) u_i_buff[i] = 0;
    }
    set->u_i = new SparseVector;
  } else {
    dot_w_u = (sum_w->dot(*set->u_i)-sum_w->dot(*set->psi_gt)*set->alpha)/(sum_w_scale);
    L_i = set->D_i+dot_w_u;
  }
  double dot_w_u_orig = dot_w_u, u_i_sqr_orig = set->u_i_sqr;
  
  // For computational reasons, we want each iteration of the main loop over set->samples to run
  // in time proportional to the number of non-zero entries in v=set->samples[j]->psi as
  // opposed to sizePsi, since v will often be very sparse
  // (e.g., for multiclass classification and mixture models it is very sparse).  Our vector
  // algebra operations have the following properties:
  //   1) Scalar multiplication takes time proportional to the number of non-zero entries
  //   2) Addition, subtraction, and dot products of 2 sparse vectors takes time proportional
  //      to the sum of the number of non-zero entries in both vectors
  //   3) Addition, subtraction, and dot products of a sparse vector and a non-sparse vector
  //      are faster, taking time proportional to the number of non-zero entries in the sparse vector
  // Since sum_w and u_i will not be sparse, we want to avoid doing any scalar multiplies (1) of sum_w
  // and u_i, and we want to avoid any operations between sum_w and u_i (2).  We therefore decompose 
  // w_sum and u_i as follows:
  //    u_i = s_u*u_i_scaled
  //    sum_w = sum_w_without_u_i - u_i
  //    u_i_sqr = <u_i,u_i>
  // We can maintain computation of u_i_scaled, s_u, sum_w_minus_u_i, u_i_sqr, L_i, and dot_w_u
  // using only operation (3) in each iteration.  In each iteration, it requires computing only one 
  // sparse dot product:
  //   <u_i_scaled,v>
  // and one sparse vector addition and scalar multiply:
  //   u_i_scaled = u_i_scaled + v*(step_size/s_u)
  // All other operations are arithmetic operations between scalars.
  SparseVector *sum_w_without_u_i = sum_w;
  SparseVector *u_i_scaled = set->u_i;
  if(set->alpha || runMultiThreaded) {
    if(set->alpha) {
      *sum_w_without_u_i += *set->u_i;
      *sum_w_without_u_i -= (*set->psi_gt*set->alpha);
    }
    double score_gt_without_u = sum_w_without_u_i->dot(*set->psi_gt);
    for(j = 0; j < set->num_samples; j++) 
      set->samples[j].dot_w = sum_w_without_u_i->dot(*set->samples[j].psi)-score_gt_without_u;
  }
  VFLOAT D_i_orig = set->D_i;
  VFLOAT dot_u_w_without_u = dot_w_u + set->u_i_sqr/(sum_w_scale);

  if(set->alpha && set->u_i_sqr) {
    // Can compute the optimal scale on u_i inexpensively
    s_u = 1 + (L_i*sum_w_scale) / set->u_i_sqr;
    s_u = my_min(1.0/set->alpha, my_max(0, s_u));
    
    set->alpha *= s_u;
    dot_u_w_without_u *= s_u;
    dot_w_u = s_u*dot_w_u - (s_u*(s_u-1)*set->u_i_sqr) / (sum_w_scale);
    set->u_i_sqr *= SQR(s_u);
    set->dot_u_psi_gt *= s_u;
    set->D_i *= s_u;
    L_i = dot_w_u + set->D_i;

    set->drift_bits += s_u <= 0 ? MAX_DRIFT_BITS : my_abs(LOG2(s_u));
    if(set->drift_bits > MAX_DRIFT_BITS) {
      // To avoid numerical precision issues, recompute set->u_i_sqr, set->dot_u_psi_gt
      set->drift_bits = 0;
      *set->u_i *= s_u;
      double d_u_gt = set->u_i->dot(*set->psi_gt);
      set->dot_u_psi_gt = d_u_gt - set->alpha*set->psi_gt_sqr;
      set->u_i_sqr = set->u_i->dot(*set->u_i) - 2*set->alpha*d_u_gt + SQR(set->alpha)*set->psi_gt_sqr;
      for(int jj = 0; jj < set->num_samples; jj++)
        set->samples[jj].alpha *= s_u;
      s_u = 1;
    }
  }
  
  u_i_scaled->make_non_sparse(true, sizePsi, true, u_i_buff);  // convert to a non-sparse vector, to accelerate -= and dot() operations
  for(r = 0; r < R; r++) {
    for(j = 0; j < set->num_samples; j++) {
      SVM_cached_sample *s = &set->samples[j];

      VFLOAT dot_u_v = u_i_scaled->dot(*s->psi)*s_u - set->dot_u_psi_gt - set->alpha*s->dot_psi_gt_psi;  // <u_i,v>
      VFLOAT dot = s->dot_w - dot_u_v;   // (lambda*t)<w,v>

#ifdef DEBUG_MULTI_SAMPLE_UPDATE
      SparseVector w_sum_new = *sum_w + (*set->psi_gt*set->alpha) - (*u_i_scaled*s_u);
      SparseVector u_new = *u_i_scaled*s_u - (*set->psi_gt*set->alpha);
      double dot_u_v_real = u_new.dot(*s->psi - *set->psi_gt);
      double dot_real = w_sum_new.dot(*s->psi - *set->psi_gt);
      fprintf(stderr, "t=%d, i=%d, j=%d, dot_u_v=%lg:%lg, dot=%lg:%lg\n", (int)t, set->i, j, dot_u_v_real, dot_u_v, dot_real, dot);
      assert(!dot_u_v || (dot_u_v_real/dot_u_v > .999999999 && dot_u_v_real/dot_u_v < 1.00000001));
      assert(!dot || (dot_real/dot > .999999999 && dot_real/dot < 1.00000001));
#endif

      VFLOAT scale=1, new_alpha;
      VFLOAT dalpha = (dot + s->loss*(sum_w_scale)) / my_max(s->sqr,.0000000001);
      if(set->alpha < 1 || dalpha < 0) {
	// alpha expand: solve for the optimal amount 'dalpha' to increase s->alpha
	// and then scale all set->samples[:].alpha (the value of 'dalpha' and 'scale' 
	// that maximizes the increase in the dual objective). 
	dalpha = my_min(1-set->alpha, my_max(-s->alpha*s_u,dalpha));
	if(set->u_i_sqr && set->alpha) {
	  scale = 1 + (L_i*sum_w_scale - dalpha*dot_u_v) / set->u_i_sqr;
	  scale = my_min((1-dalpha)/set->alpha, my_max(0, scale));
	  if(s->alpha*s_u*scale+dalpha < 0) dalpha = -s->alpha*s_u*scale;
	}
	new_alpha = set->alpha*scale + dalpha;
      } else {
	// alpha swap: solve for the optimal amount 'dalpha' to increase s->alpha
	// while scaling down all set->samples[:].alpha, such that we preserve sum_k{s->samples[:].alpha}=1
	// (choose the value of 'dalpha' that maximizes the increase in the dual objective)
	VFLOAT e = dot/(sum_w_scale) + s->loss;
	VFLOAT sqr = set->u_i_sqr + 2*dot_u_v + s->sqr;
	dalpha = (e-L_i)*(sum_w_scale) / my_max(sqr,.00000000001);
	dalpha = my_min(1-s->alpha*s_u, my_max(-s->alpha*s_u,dalpha));
	scale = 1-dalpha;
	//assert(scale > 0 && scale <= 1);
	new_alpha = 1;
      }
      if(!(scale >= 0 && new_alpha >= 0 && new_alpha <= 1.000000001)) {
	fprintf(stderr, "scale=%f, new_alpha=%f, dalpha=%f, set->alpha=%f, u=%f\n", (float)scale, (float)new_alpha, (float)dalpha, (float)set->alpha, (float)set->u_i_sqr);
      }
      assert(scale >= 0 && new_alpha >= -0.000000001 && new_alpha <= 1.000000001);
      new_alpha = my_min(1,my_max(new_alpha,0));

      if(dalpha != 0 || scale != 1) {
	s_u *= scale;
	set->drift_bits += s_u <= 0 ? MAX_DRIFT_BITS : my_abs(LOG2(scale));
	if(set->drift_bits >= MAX_DRIFT_BITS) {
	  // To avoid numerical precision problems, set s_u=1 and recompute some values manually
	  set->drift_bits = 0;
	  u_i_scaled->make_non_sparse(false, -1, false, u_i_buff);
	  *u_i_scaled *= s_u;
	  double d_u_gt = u_i_scaled->dot(*set->psi_gt), alpha = set->alpha*scale;
	  set->dot_u_psi_gt = d_u_gt - alpha*set->psi_gt_sqr;
	  set->u_i_sqr = u_i_scaled->dot(*u_i_scaled) - 2*alpha*d_u_gt + SQR(alpha)*set->psi_gt_sqr;
	  set->D_i = scale*set->D_i;
	  dot_u_w_without_u = sum_w_without_u_i->dot(*u_i_scaled-(*set->psi_gt*alpha))/(sum_w_scale);
	  dot_w_u = dot_u_w_without_u - set->u_i_sqr/(sum_w_scale);
	  dot_u_v = u_i_scaled->dot(*s->psi) - set->dot_u_psi_gt - alpha*s->dot_psi_gt_psi; 
	  dot = s->dot_w - dot_u_v;
	  for(int jj = 0; jj < set->num_samples; jj++)
	    set->samples[jj].alpha *= s_u;
	  scale = 1;
	  s_u = 1;
	  u_i_scaled->make_non_sparse(true, sizePsi, true, u_i_buff);
	}

	*u_i_scaled += (*s->psi * (dalpha/s_u));
	s->alpha += dalpha/s_u;

	// Keep track of L_i, D_i, u_i_sqr, dot_w_u, dot_u_psi_gt using inexpensive online updates
	set->alpha = new_alpha;
	dot_u_w_without_u = scale*dot_u_w_without_u + dalpha*s->dot_w/(sum_w_scale);
        dot_w_u = scale*dot_w_u + (dalpha*dot - scale*(scale-1)*set->u_i_sqr - (2*scale*dalpha-dalpha)*dot_u_v - s->sqr*SQR(dalpha)) / (sum_w_scale);
	set->u_i_sqr = SQR(scale)*set->u_i_sqr + 2*scale*dalpha*dot_u_v + s->sqr*SQR(dalpha);
	set->dot_u_psi_gt = scale*set->dot_u_psi_gt + dalpha*(s->dot_psi_gt_psi - set->psi_gt_sqr);
	set->D_i = scale*set->D_i + dalpha*s->loss;
	L_i = dot_w_u + set->D_i;
	assert(!isnan(L_i));
#ifdef DEBUG_MULTI_SAMPLE_UPDATE
	SparseVector w_sum_new = *sum_w + (*set->psi_gt*set->alpha) - (*u_i_scaled*s_u);
	SparseVector u_new = *u_i_scaled*s_u - (*set->psi_gt*set->alpha);
	double dot_w_u_real = w_sum_new.dot(u_new)/(sum_w_scale);
	double dot_u_psi_gt_real = u_new.dot(*set->psi_gt);
	double u_i_sqr_real = u_new.dot(u_new);
	double dot_u_w_without_u_real = sum_w->dot(u_new)/(sum_w_scale);
	fprintf(stderr, "t=%d, i=%d, j=%d, scale=%f, dalpha=%f, s_u=%f, dot_w_u=%lg:%lg, dot_u_psi_gt=%lg:%lg, u_i_sqr=%lg:%lg\n", (int)t, set->i, j, (float)scale, (float)dalpha, (float)s_u, dot_w_u_real, dot_w_u,  dot_u_psi_gt_real, set->dot_u_psi_gt,  u_i_sqr_real, set->u_i_sqr);
	assert(dot_w_u_real/dot_w_u > .999999999 && dot_w_u_real/dot_w_u < 1.00000001);
	assert(dot_u_psi_gt_real/set->dot_u_psi_gt > .999999999 && dot_u_psi_gt_real/set->dot_u_psi_gt < 1.00000001);
	assert(u_i_sqr_real/set->u_i_sqr > .999999999 && u_i_sqr_real/set->u_i_sqr < 1.00000001);
	if(dot_u_w_without_u) assert(dot_u_w_without_u_real/dot_u_w_without_u > .999999999 && dot_u_w_without_u_real/dot_u_w_without_u < 1.00000001);
#endif
      }
    }
  }

  
  // Update sum_w, sum_dual, sum_w_sqr, and regularization_error, taking into account the new value of u_i
  set->u_i->make_non_sparse(false, -1, false, u_i_buff);
  *set->u_i *= s_u;
  *sum_w -= *set->u_i;   // Add u_i back into (lambda*t)w
  *sum_w += (*set->psi_gt*set->alpha);
  double d_sum_w_sqr = 2*(dot_w_u_orig-dot_u_w_without_u)*sum_w_scale + set->u_i_sqr + u_i_sqr_orig;
  sum_dual += -d_sum_w_sqr/(2*sum_w_scale) + set->D_i - D_i_orig;
  sum_alpha_loss += set->D_i - D_i_orig;
  sum_w_sqr += d_sum_w_sqr;
  regularization_error = (sum_w_sqr/SQR(sum_w_scale))*lambda/2;
  for(j = 0; j < set->num_samples; j++)
    set->samples[j].alpha *= s_u;
  set->slack_after = set->num_samples ? sum_w->dot(*set->samples[0].psi-*set->psi_gt)/(sum_w_scale) + set->samples[0].loss : 0;
  

#ifdef DEBUG_MULTI_SAMPLE_UPDATE
  double sum_w_sqr_real = sum_w->dot(*sum_w);
  fprintf(stderr, "t=%d, i=%d, sum_w_sqr=%lg:%lg\n", (int)t, set->i, sum_w_sqr_real, sum_w_sqr);
  assert(sum_w_sqr_real/sum_w_sqr > .99999 && sum_w_sqr_real/sum_w_sqr < 1.0001);
#endif

  //fprintf(stderr, "sum dual is %lg, dual_change=%lg, D=%lg\n", sum_dual, -d_sum_w_sqr/(2*sum_w_scale) + set->D_i - D_i_orig, set->D_i);
}




bool StructuredSVM::SaveOnlineData(const char *fname) {
  FILE *fout = fopen(fname, "wb");
  if(!fout) return false;

  long tm = GetElapsedTime();
  bool b = (fwrite(&curr, sizeof(long), 1, fout) &&
         fwrite(&minItersBeforeNewExample, sizeof(long), 1, fout) &&
         fwrite(&M, sizeof(int), 1, fout) &&
         fwrite(&last_example, sizeof(int), 1, fout) &&
         fwrite(&sum_iter_error, sizeof(double), 1, fout) &&
         fwrite(&sum_iter_error_window, sizeof(double), 1, fout) &&
         fwrite(&sum_generalization_error, sizeof(double), 1, fout) &&
         fwrite(&sum_generalization_error_window, sizeof(double), 1, fout) &&
         fwrite(&regularization_error, sizeof(double), 1, fout) &&
         fwrite(&tm, 1, sizeof(long), fout) &&
         fwrite(&n, sizeof(int), 1, fout) &&
         fwrite(&sum_dual, sizeof(double), 1, fout) &&
         fwrite(&sum_alpha_loss, sizeof(double), 1, fout) &&
         fwrite(&sum_w_sqr, sizeof(double), 1, fout) &&
         fwrite(&sum_w_scale, sizeof(double), 1, fout) &&
         fwrite(ex_num_iters, sizeof(int), trainset->num_examples, fout) &&
         fwrite(ex_first_iter, sizeof(int), trainset->num_examples, fout));
  assert(b);

  if(generalization_errors_by_n) (fwrite(generalization_errors_by_n, sizeof(double), n, fout));
  if(generalization_errors_by_t) (fwrite(generalization_errors_by_t, sizeof(double), t, fout));
  if(iter_errors_by_t) (fwrite(iter_errors_by_t, sizeof(double), t, fout));
  if(sum_dual_by_t) (fwrite(sum_dual_by_t, sizeof(double), t, fout));
  if(regularization_errors_by_t) (fwrite(regularization_errors_by_t, sizeof(double), t, fout));
  if(losses_by_t) (fwrite(losses_by_t, sizeof(double), t, fout));
  if(elapsed_time_by_t) (fwrite(elapsed_time_by_t, sizeof(double), t, fout));

  if(iter_examples) (fwrite(iter_examples, sizeof(long), t, fout));
  if(examples_by_iteration_number) (fwrite(examples_by_iteration_number, sizeof(int), M+1, fout));
  if(examples_by_iteration_next_ind) (fwrite(examples_by_iteration_next_ind, sizeof(int), trainset->num_examples, fout));
  if(examples_by_iteration_prev_ind) (fwrite(examples_by_iteration_prev_ind, sizeof(int), trainset->num_examples, fout));

  char cache_name[1000];  
  sprintf(cache_name, "%s.cache", fname);
  SaveCachedExamples(cache_name);

  fclose(fout);
  return true;
}

bool StructuredSVM::LoadOnlineData(const char *fname) {
  FILE *fin = fopen(fname, "rb");
  if(!fin) return false;

  if(trainfile)
    assert((trainset=LoadDataset(trainfile)) != NULL);

  alloc_n = (trainset->num_examples+1);
  ex_num_iters = (int*)realloc(ex_num_iters, sizeof(int)*alloc_n);
  bool b = (fread(&curr, sizeof(long), 1, fin) &&
         fread(&minItersBeforeNewExample, sizeof(long), 1, fin) &&
         fread(&M, sizeof(int), 1, fin) &&
         fread(&last_example, sizeof(int), 1, fin) &&
         fread(&sum_iter_error, sizeof(double), 1, fin) &&
         fread(&sum_iter_error_window, sizeof(double), 1, fin) &&
         fread(&sum_generalization_error, sizeof(double), 1, fin) &&
         fread(&sum_generalization_error_window, sizeof(double), 1, fin) &&
         fread(&regularization_error, sizeof(double), 1, fin) &&
         fread(&base_time, 1, sizeof(double), fin) &&
         fread(&sum_dual, sizeof(double), 1, fin) &&
         fread(&sum_alpha_loss, sizeof(double), 1, fin) &&
         fread(&sum_w_sqr, sizeof(double), 1, fin) &&
         fread(&sum_w_scale, sizeof(double), 1, fin) &&
         fread(&n, sizeof(int), 1, fin) &&
         fread(ex_num_iters, sizeof(int), trainset->num_examples, fin) &&
         fread(ex_first_iter, sizeof(int), trainset->num_examples, fin));
  assert(b);

  alloc_t = t;
  generalization_errors_by_n = (double*)realloc(generalization_errors_by_n, sizeof(double)*alloc_n);
  generalization_errors_by_t = (double*)realloc(generalization_errors_by_t, sizeof(double)*alloc_t);
  iter_errors_by_t = (double*)realloc(iter_errors_by_t, sizeof(double)*alloc_t);
  sum_dual_by_t = (double*)realloc(sum_dual_by_t, sizeof(double)*alloc_t);
  regularization_errors_by_t = (double*)realloc(regularization_errors_by_t, sizeof(double)*alloc_t);
  losses_by_t = (double*)realloc(losses_by_t, sizeof(double)*alloc_t);
  elapsed_time_by_t = (double*)realloc(elapsed_time_by_t, sizeof(double)*alloc_t);
  iter_examples = (long*)realloc(iter_examples, sizeof(long)*alloc_t);
  examples_by_iteration_number = (int*)realloc(examples_by_iteration_number, sizeof(int)*(M+2));
  examples_by_iteration_next_ind = (int*)realloc(examples_by_iteration_next_ind, sizeof(int)*alloc_n);
  examples_by_iteration_prev_ind = (int*)realloc(examples_by_iteration_prev_ind, sizeof(int)*alloc_n);
  examples_by_iteration_number[M+1] = -1;

  b = (fread(generalization_errors_by_n, sizeof(double), n, fin)); assert(b);
  b = (fread(generalization_errors_by_t, sizeof(double), t, fin)); assert(b);
  b = (fread(iter_errors_by_t, sizeof(double), t, fin)); assert(b);
  b = (fread(sum_dual_by_t, sizeof(double), t, fin)); assert(b);
  b = (fread(regularization_errors_by_t, sizeof(double), t, fin)); assert(b);
  b = (fread(losses_by_t, sizeof(double), t, fin)); assert(b);
  b = (fread(elapsed_time_by_t, sizeof(double), t, fin)); assert(b);
  b = (fread(iter_examples, sizeof(long), t, fin)); assert(b);
  b = (fread(examples_by_iteration_number, sizeof(int), M+1, fin)); assert(b);
  b = (fread(examples_by_iteration_next_ind, sizeof(int), trainset->num_examples, fin)); assert(b);
  b = (fread(examples_by_iteration_prev_ind, sizeof(int), trainset->num_examples, fin)); assert(b);

  char cache_name[1000];  
  sprintf(cache_name, "%s.cache", fname);
  LoadCachedExamples(cache_name);

  fclose(fin);
  return true;
}

void StructuredSVM::LoadTrainset(const char *fname) {
  if(trainfile) free(trainfile);
  trainfile = StringCopy(fname);

  trainset = LoadDataset(fname);
  assert(trainset != NULL);

  int i;
  examples_by_iteration_number = (int*)realloc(examples_by_iteration_number, sizeof(int)*(M+2));
  for(i = 0; i <= M+1; i++)
    examples_by_iteration_number[i] = -1;
  for(i = 0; i < trainset->num_examples; i++) {
    CreateTrainingExampleQueues(i);
  }
}

void StructuredSVM::ExtractSampleSet(int num_per_negative, bool augment) {
  StructuredDataset *train = GetTrainset(); 

  max_samples = num_per_negative;

  Lock();
  if(augment) {
    for(int i = 0; i < n; i++) 
      if(trainset->examples[i]->set) 
	CondenseSamples(trainset->examples[i]->set);
  }
  Unlock();

  SparseVector *w = sum_w ? GetCurrentWeights(false) : NULL;
  int num = 0;

#pragma omp parallel for
  for(int i = 0; i < train->num_examples; i++) {
    if(!augment) {
      // Sample negative examples randomly
      assert(!trainset->examples[i]->set);
      trainset->examples[i]->set = new_SVM_cached_sample_set(i, Psi(train->examples[i]->x, train->examples[i]->y).ptr());
      trainset->examples[i]->set->psi_gt_sqr = trainset->examples[i]->set->psi_gt->dot(*trainset->examples[i]->set->psi_gt);
    }
    num++;
    assert(trainset->examples[i]->set);
    trainset->examples[i]->set->score_gt = sum_w ? sum_w->dot(*trainset->examples[i]->set->psi_gt)/sum_w_scale : 0; 
    fprintf(stderr, "Extracting sample %d of %d...\n", num, train->num_examples);

    // Mine hard negative examples.  If the model is uninitialized, this instead draws random
    // negative examples.  Features are extracted for each sample
    ImportanceSample(train->examples[i]->x, w, train->examples[i]->y, trainset->examples[i]->set, 1);
    SVM_cached_sample_set_compute_features(trainset->examples[i]->set, train->examples[i]);
    OnFinishedIteration(train->examples[i]->x, train->examples[i]->y);
  }

  n = train->num_examples;
  SetSumWScale(lambda*n);
  delete w;
}

void StructuredSVM::SaveCachedExamples(const char *output_name, bool saveFull) {
  FILE *fout = fopen(output_name, "wb");
  fwrite(&n, sizeof(long), 1, fout);
  for(int i = 0; i < n; i++) {
    bool b = trainset->examples[i]->set ? true : false;
    fwrite(&b, sizeof(bool), 1, fout);
    if(b)
      write_SVM_cached_sample_set(trainset->examples[i]->set, fout, this, saveFull);
  }
  fclose(fout);
}


void StructuredSVM::LoadCachedExamples(const char *fname, bool loadFull) {
  int i;
  for(i = 0; i < n; i++)
    if(trainset->examples[i]->set)
      free_SVM_cached_sample_set(trainset->examples[i]->set);
  trainset->examples[i]->set = NULL;
  FILE *fin = fopen(fname, "rb");
  if(fin) {
    long n2;
    int b = fread(&n2, sizeof(long), 1, fin); assert(b);
    assert(n2 == n);
    for(i = 0; i < n; i++) {
      bool b;
      fread(&b, sizeof(bool), 1, fin);
      if(b)
	trainset->examples[i]->set = read_SVM_cached_sample_set(fin, this, trainset->examples[i]->x, loadFull);
    }
    fclose(fin);
  }
}

void StructuredSVM::CondenseSamples(SVM_cached_sample_set *set) {
  if(!set->u_i && !set->alpha && set->num_samples && method != SPO_MAP_TO_BINARY && method != SPO_MAP_TO_BINARY_MINE_HARD_NEGATIVES) {
    set->alpha = set->samples[0].alpha;	
    set->u_i = set->samples[0].psi;
    set->samples[0].psi = NULL;
    *set->u_i *= set->samples[0].alpha;	
    set->D_i = set->samples[0].alpha*set->samples[0].loss;
    set->dot_u_psi_gt = set->samples[0].alpha*set->samples[0].dot_psi_gt_psi;
    set->u_i_sqr = SQR(set->samples[0].alpha)*set->samples[0].sqr;
  }

  if(set->num_samples > maxCachedSamplesPerExample) {
    qsort(set->samples, set->num_samples, sizeof(SVM_cached_sample), SVM_cached_sample_cmp);

    for(int i = maxCachedSamplesPerExample; i < set->num_samples; i++)
      free_SVM_cached_sample(&set->samples[i]);
    set->num_samples = maxCachedSamplesPerExample;
  }
}	  

void StructuredSVM::ConvertCachedExamplesToBinaryTrainingSet() {
  long num = 0;
  StructuredDataset *binary_trainset = new StructuredDataset;

  for(int i = 0; i < n; i++) {
    for(int j = 0; j < trainset->examples[i]->set->num_samples; j++) {
      binary_trainset->examples[num]->set = new_SVM_cached_sample_set(num, trainset->examples[i]->set->psi_gt);
      *binary_trainset->examples[num]->set = *trainset->examples[i]->set;
      binary_trainset->examples[num]->set->i = num;
      binary_trainset->examples[num]->set->num_samples = 1;
      binary_trainset->examples[num++]->set->samples = trainset->examples[i]->set->samples+j;
      binary_trainset->AddExample(CopyExample(trainset->examples[i]->x, trainset->examples[i]->y, trainset->examples[i]->y_latent));
    }
  }
  
  trainset = binary_trainset;
  n = num;
  SetSumWScale(lambda*n);
}

double StructuredSVM::ImportanceSample(StructuredData *x, SparseVector *w, StructuredLabel *y_gt, struct _SVM_cached_sample_set *set, double w_scale) {
  //char str[1000]; OptimizationMethodToString(method, str);
  //fprintf(stderr, "ERROR: method %s is not supported, because ImportanceSample() is not implemented\n", str);
  //exit(0);
    
  // Usually, the user of this API should define a custom importance sampling routine, which adds a set of predicted labels with non-zero slack.  
  // If unimplemented, we simply add a single most violated label in each iteration, which is appended to the set of most violated labels
  // from previous iterations.
  StructuredLabel *ybar = NewStructuredLabel(x);
  double retval = Inference(x, ybar, w, NULL, y_gt, w_scale);
  SVM_cached_sample_set_add_sample(set, ybar);
  SVM_cached_sample s = set->samples[0];
  set->samples[0] = set->samples[set->num_samples-1];
  set->samples[set->num_samples-1] = s;

  return retval;
}

int int_compare (const void * a, const void * b)
{
  return ( *(int*)a - *(int*)b );
}

int SVM_cached_sample_cmp(const void * a, const void * b) {
  double d = ((SVM_cached_sample*)b)->slack - ((SVM_cached_sample*)a)->slack;
  return d < 0 ? -1 : (d > 0 ? 1 : 0);
}


int SVM_cached_sample_set_ave_slack_cmp(const void *a, const void *b) {
  SVM_cached_sample_set *set1 = *((SVM_cached_sample_set**)a);
  SVM_cached_sample_set *set2 = *((SVM_cached_sample_set**)b);
  double v1 = set1 && set1->numIters ? set1->sumSlack/set1->numIters : -10000000;
  double v2 = set2 && set2->numIters ? set2->sumSlack/set2->numIters : -10000000;
  double d = v2-v1;
  return d < 0 ? -1 : (d > 0 ? 1 : 0);
}

int SVM_cached_sample_set_alpha_cmp(const void *a, const void *b) {
  SVM_cached_sample_set *set1 = *((SVM_cached_sample_set**)a);
  SVM_cached_sample_set *set2 = *((SVM_cached_sample_set**)b);
  double v1 = set1 ? set1->alpha : 0;
  double v2 = set2 ? set2->alpha : 0;
  double d = v2-v1;
  return d < 0 ? -1 : (d > 0 ? 1 : SVM_cached_sample_set_ave_slack_cmp(a,b));
}
