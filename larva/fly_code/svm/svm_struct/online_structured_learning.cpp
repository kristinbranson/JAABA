#include "online_structured_learning.h"
#include "../../blob.h"

#ifdef WIN32
#include <windows.h>
//void Sleep(unsigned long  dwMilliseconds);
void usleep(double us) { Sleep(us*1000); }
#else
#include <unistd.h>
#endif

//#define MAX_THREADS 2

/*
 * The rules for choosing which example to process next are (from highest to lowest precedence):
 *   1) Never allow different threads to process the same example at the same time
 *   2) Prefer processing an example that has been iterated over at least once but less than minItersBeforeNewExample
 *      iterations before going onto a new example (allows better estimation of generalization and optimization errors)
 *   3) Otherwise, prefer processing an example that has been iterated over the fewest times
 */
int StructuredSVMOnlineLearner::ChooseNextExample() {
  // Choose the examples in order, selecting the one that has been processed the fewest number of times
  int it = currMinIterByExample;
  while(it <= M && examples_by_iteration_number[it] < 0)
    it++;
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

void StructuredSVMOnlineLearner::Train() {
  cache_old_examples = sparm->method == SPO_DUAL_UPDATE_WITH_CACHE ||
                            sparm->method == SPO_DUAL_MULTI_SAMPLE_UPDATE_WITH_CACHE; 
  //lambda = .3;
  minItersBeforeNewExample = 5;

  /* some training information */
  fprintf(stderr, "Number of threads=%d\n", num_thr);
  printf("Regularization constant (lambda): %.8g\n", lambda);
  printf("Approximation factor (epsilon): %.8g\n", learn_parm->eps);
  printf("Number of training examples (n): %ld\n", (long)sample.n); 
  printf("Feature space dimension (sizePsi): %ld\n", sm.sizePsi); fflush(stdout);

  start_time = time(NULL);

  #pragma omp parallel num_threads(num_thr) 
  {
    int tid = omp_get_thread_num();
    if(tid == 0 && cache_old_examples) {
      // If statement occurs only if sparm->method == SPO_DUAL_UPDATE_WITH_CACHE || sparm->method == SPO_MULTI_SAMPLE_DUAL_UPDATE_WITH_CACHE
      // Optimization thread, continuously optimize dual parameters for labels for previous cached examples (earlier calls
      // to find_most_violated_constraint)
      int num = 0;
      while(!finished) {
	usleep(1);
	while(t <= num_thr) // wait for at least one first constraint
	  usleep(100000);

	// Choose a label from a random iteration (a label ybar from a previous worker thread iteration), and optimize its
	// dual parameters
	omp_set_lock(&my_lock);
	int i = rand()%t;
	if(!cached_examples[i]) {
	  omp_unset_lock(&my_lock);
	  continue;
	}
	UpdateWeights(cached_examples[i], i);
	num++;
	omp_unset_lock(&my_lock);

	if(num >= t*50) { 
	  // Sanity check: occasionally recompute weights from dual parameters.  Could help avoid drifting due to numerical precision errors
	  RecomputeWeights(); 
	  num = 0; 
	}  
      }
    } else {
      // Worker threads, continuously call ybar=find_most_violated_constraint and then use ybar to update the current weights
      int i = -1;
      while(!finished) {
	STRUCTMODEL sm2 = sm;

	// Choose a training example 'i' to process
	omp_set_lock(&my_lock);
	SVECTOR *w = smult_s(sum_w, 1.0/(lambda*t));
	sm2.w = create_nvector(sm.sizePsi, w);
	i = ChooseNextExample();
	if(i < 0) {
	  usleep(100000);
	  omp_unset_lock(&my_lock);
	  continue;
	}
    LABEL y = sample.examples[i].y;
    EXAMPLE ex_curr = sample.examples[i];
	omp_unset_lock(&my_lock);
    
    
	// Find the most violated label ybar = max_y <w_t,psi(x_i,y)>+l(y_i,y) 
	SVM_cached_sample_set *set = new_SVM_cached_sample_set(i);
	double score_loss;
	LABEL ybar = m->find_most_violated_constraint_marginrescaling(&ex_curr.x, &ex_curr.y, &sm2, sparm, &score_loss);
	SVECTOR *fybar = m->psi(&ex_curr.x,&ybar,&sm2,sparm);
	double lossval = m->loss(ex_curr.y,ybar,sparm);
	SVECTOR *y_psi = m->psi(&ex_curr.x, &ex_curr.y, &sm, sparm);
	SVECTOR *dpsi_s = multadd_ss(fybar, y_psi, 1, -1);  // gradient <w_t,psi(x_i,y)>-<w_t,psi(x_i,y_i)> (not including regularization)
	double score = score_loss-lossval;             // <w_t,psi(x_i,y)>
	double score_gt = sprod_ns(sm2.w, y_psi);   //score_loss-slack;            // <w_t,psi(x_i,y_i)>
	assert(lossval >= 0);

	//release_memory_caches(&sample.examples[i].x, &sample.examples[i].y, sparm);

	// Use ybar to update the weights
	omp_set_lock(&my_lock);
    sample.examples[i].x = ex_curr.x; sample.examples[i].y = ex_curr.y;
	SVM_cached_sample_set_add_sample(set, ybar, dpsi_s, lossval);
	UpdateWeights(set, -1);
	omp_unset_lock(&my_lock);

	fprintf(stderr, "Example %d: m=%d slack=%f score=%f score_gt=%f loss=%f alpha=%f\n", 
		i, ex_num_iters[i], (float)(score_loss-score_gt), (float)(score), (float)(score_gt), (float)lossval, (float)set->samples[0].alpha);

	char ename[1000]; sprintf(ename, "%d", i);
	m->on_finished_find_most_violated_constraint(&ybar, &sample.examples[i].y, ex_num_iters[i], sparm, ename);

	// Cleanup
	free_svector(y_psi);
	free_svector(fybar);
	free(sm2.w);
	free_svector(w);
	if(!cache_old_examples) 
	  free_SVM_cached_sample_set(set, m);
      }
    }
  }

  SaveModel(modelfile);
}

void StructuredSVMOnlineLearner::UpdateWeights(SVM_cached_sample_set *ex, int iterInd) {
  long tt = t;
  if(iterInd == -1) {
    iter_examples = (long*)realloc(iter_examples, sizeof(long)*(t+1));
    iter_examples[t] = ex->i;

    // Book-keeping stuff, for estimating generalization error, optimization error, and regret
    double e = sprod_ss(sum_w, ex->samples[0].dpsi)/(lambda*(t+1)) + ex->samples[0].loss;
    if(ex_num_iters[ex->i] == 0) {
      sum_generalization_error += my_max(e,0);
      generalization_errors_by_n[ex->i] = optimization_errors_by_n[ex->i] = model_errors_by_n[ex->i] = my_max(e,0);
      losses_by_n[ex->i] = ex->samples[0].loss;
      n++;
    } 
    iter_errors_by_t = (double*)realloc(iter_errors_by_t, sizeof(double)*(t+1));
    generalization_errors_by_t = (double*)realloc(generalization_errors_by_t, sizeof(double)*(t+1));
    constraint_errors_by_t = (double*)realloc(constraint_errors_by_t, sizeof(double)*(t+1));
    model_errors_by_t = (double*)realloc(model_errors_by_t, sizeof(double)*(t+1));
    regularization_errors_by_t = (double*)realloc(regularization_errors_by_t, sizeof(double)*(t+1));
    losses_by_t = (double*)realloc(losses_by_t, sizeof(double)*(t+1));
    elapsed_time_by_t = (double*)realloc(elapsed_time_by_t, sizeof(double)*(t+1));

    sum_iter_error += my_max(e,0);
    iter_errors_by_t[t] = constraint_errors_by_t[t] = optimization_errors_by_n[ex->i] = my_max(e,0);
    generalization_errors_by_t[t] = generalization_errors_by_n[ex->i];
    losses_by_t[t] = ex->samples[0].loss;
    elapsed_time_by_t[t] = (double)GetElapsedTime();

    ex_num_iters[ex->i]++;
    if(ex_num_iters[ex->i] > M) {
      assert(M == ex_num_iters[ex->i]-1);
      M++;
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
    

    tt = t++;
    if(cache_old_examples) {
      cached_examples = (SVM_cached_sample_set**)realloc(cached_examples, (t)*sizeof(SVM_cached_sample_set*));
      cached_examples[tt] = ex;
    }
  }
  switch(sparm->method) {
    case SPO_SGD:
    case SPO_SGD_PEGASOS:
      // Stochastic gradient descent: take a step of size -1/(lambda*t) in the direction of psi(ybar,x)-psi(y_i,x)
      assert(iterInd == -1 && ex->num_samples == 1);
      sum_w = SVM_cached_sample_set_update_parameters(&ex->samples[0], sum_w, 1, 1, lambda, t);

      if(sparm->method == SPO_SGD_PEGASOS) {
	// Project sum_w onto the L2 ball, such that ||w||^2 < 1/lambda
	double m = sprod_ss(sum_w,sum_w)/SQR(lambda*t);
	if(m > 1/lambda) {
	  SVECTOR *sum_w_new = smult_s(sum_w, 1/lambda / m);
	  free_svector(sum_w);
	  sum_w = sum_w_new;
	}
      }
      break;

    case SPO_DUAL_UPDATE:
    case SPO_DUAL_UPDATE_WITH_CACHE:
      // Take a step in the direction of psi(ybar,x)-psi(y_i,x), where the chosen step size maximizes the dual objective
      assert(ex->num_samples == 1);
      sum_w = SVM_cached_sample_optimize_dual(&ex->samples[0], sum_w, lambda, t);
      break;

    case SPO_DUAL_MULTI_SAMPLE_UPDATE:
    case SPO_DUAL_MULTI_SAMPLE_UPDATE_WITH_CACHE:
      // Jointly optimize the dual parameters for a set of different labels within this example
      sum_w = SVM_cached_sample_set_optimize_dual(ex, sum_w, lambda, t, R);
      break;

    default:
      fprintf(stderr, "svm_learn_main not implemented: %d\n", sparm->method);
      assert(0);  
  }

  double e = sprod_ss(sum_w, ex->samples[0].dpsi)/(lambda*t) + ex->samples[0].loss;
  if(iterInd == -1) {
    // Book-keeping stuff, for estimating generalization error, optimization error, regret, and convergence

    // Error measured over entire training history
    constraint_errors_by_t[tt] = my_max(e,0);
    sum_model_error += my_max(e,0);
    model_errors_by_n[ex->i] = model_errors_by_t[tt] = sum_model_error/t;
    regularization_errors_by_t[tt] = regularization_errors_by_n[ex->i] = 
      regularization_error = sprod_ss(sum_w,sum_w)/SQR(lambda*t)*lambda/2;
    double eps_measured = (sum_iter_error-sum_model_error) / t;
    fprintf(stderr, "t=%d, n=%d: Average Training Error=%f (Model error=%f, Optimization error=%f, Regularization error=%f), Out of sample error=%f\n", 
	    (int)t, (int)n, (float)(sum_iter_error/t)+(float)regularization_error, (float)(sum_model_error/t), (float)((sum_iter_error-sum_model_error)/t), (float)regularization_error, (float)(sum_generalization_error/n));
    

    // Error measured over last 100 examples
    double sum_iter_error2=0, sum_model_error2=0, sum_generalization_error2=0;
    int ttt = my_min(t, 100), nnn = my_min(ex->i, 100);
    for(int i = 0; i < ttt; i++) {
      sum_iter_error2 += iter_errors_by_t[tt-i];
      sum_model_error2 += constraint_errors_by_t[tt-i];
    }
    for(int j = 0; j < nnn; j++) 
      sum_generalization_error2 += generalization_errors_by_n[ex->i-j];
    fprintf(stderr, "Last %d iters: Average Training Error=%f (Model error=%f, Optimization error=%f, Regularization error=%f), Out of sample error=%f\n", 
	    (int)ttt, (float)(sum_iter_error2/ttt)+(float)regularization_error, (float)(sum_model_error2/ttt), (float)((sum_iter_error2-sum_model_error2)/ttt), (float)regularization_error, (float)(sum_generalization_error2/nnn));
    
    //if(eps_measured < eps && t > 10000)
      //finished = true;
  } else {
    sum_model_error += my_max(e,0) - constraint_errors_by_t[iterInd];
    constraint_errors_by_t[iterInd] = my_max(e,0);
    losses_by_t[iterInd] = ex->samples[0].loss;
    if(ex_num_iters[ex->i] == M) {
      model_errors_by_n[ex->i] = my_max(e,0);
      losses_by_n[ex->i] = ex->samples[0].loss;
    }
  }
}

void StructuredSVMOnlineLearner::Init() {
  omp_init_lock(&my_lock);
  this->R = 1;
  base_time = 0;

  this->currMinIterByExample = 0;
  this->curr = 0;
  this->t = 0;
  this->n = 0;
  this->M = 0;
  this->finished = false;
  minItersBeforeNewExample = 0;
  memset(&sm, 0, sizeof(STRUCTMODEL));

  this->sum_generalization_error = this->sum_model_error = this->sum_iter_error = 0;

  this->ex_num_iters = NULL;
  this->cached_examples = NULL;
  this->sum_w = NULL;

  this->trainfile = NULL;
  this->modelfile = NULL;
  this->learn_parm = NULL;
  this->kernel_parm = NULL;
  this->sparm = NULL;

  generalization_errors_by_n = optimization_errors_by_n = model_errors_by_n = regularization_errors_by_n = losses_by_n = NULL;
  generalization_errors_by_t = iter_errors_by_t = model_errors_by_t = constraint_errors_by_t = regularization_errors_by_t = losses_by_t = elapsed_time_by_t = NULL;
  iter_examples = NULL;

  examples_by_iteration_number = NULL;
  examples_by_iteration_next_ind = NULL;
  examples_by_iteration_prev_ind = NULL;

  num_thr = omp_get_num_procs();
#ifdef MAX_THREADS
  if(num_thr > MAX_THREADS) num_thr = MAX_THREADS;
#endif
  if(num_thr < 2) num_thr = 2;
}

StructuredSVMOnlineLearner::StructuredSVMOnlineLearner(SVMStructMethod *m, const char *fname) {
  this->m = m;
  Init();
  ReadModel(fname);
}

StructuredSVMOnlineLearner::StructuredSVMOnlineLearner(SVMStructMethod *m, LEARN_PARM *learn_parm, KERNEL_PARM *kernel_parm, STRUCT_LEARN_PARM *sparm, SAMPLE *samp, const char *trainfile, const char *modelfile) {
  this->m = m;
  Init();

  this->lambda = (double)(1.0/sparm->C);
  this->eps = learn_parm->eps;

  //sample = m->read_struct_examples(trainfile, sparm); 
  this->sample = *samp;
  ex_num_iters = (int*)malloc(sizeof(int)*sample.n);
  memset(ex_num_iters, 0, sizeof(int)*sample.n);
  m->init_struct_model(sample, &sm, sparm, learn_parm, kernel_parm);
  if(!sm.w) {
    sm.w = create_nvector(sm.sizePsi+1);
  }
  clear_nvector(sm.w, sm.sizePsi);
  this->sum_w = create_svector_n(sm.w, sm.sizePsi, NULL, 1);

  generalization_errors_by_n = (double*)realloc(generalization_errors_by_n, sizeof(double)*(sample.n+1));
  optimization_errors_by_n = (double*)realloc(optimization_errors_by_n, sizeof(double)*(sample.n+1));
  model_errors_by_n = (double*)realloc(model_errors_by_n, sizeof(double)*(sample.n+1));
  regularization_errors_by_n = (double*)realloc(regularization_errors_by_n, sizeof(double)*(sample.n+1));
  losses_by_n = (double*)realloc(losses_by_n, sizeof(double)*(sample.n+1));
  
  // Build a queue defining the order examples will be processed
  int i;
  examples_by_iteration_number = (int*)realloc(examples_by_iteration_number, sizeof(int)*(M+2));
  for(i = 0; i <= M+1; i++)
    examples_by_iteration_number[i] = -1;
  if(sample.n) {
    examples_by_iteration_number[0] = 0;
    examples_by_iteration_next_ind = (int*)realloc(examples_by_iteration_next_ind, sizeof(int)*(sample.n+1));
    examples_by_iteration_prev_ind = (int*)realloc(examples_by_iteration_prev_ind, sizeof(int)*(sample.n+1));
    for(i = 0; i < sample.n; i++) {
      examples_by_iteration_prev_ind[i] = i-1;
      examples_by_iteration_next_ind[i] = i+1;
    }
    examples_by_iteration_next_ind[sample.n-1] = 0;
    examples_by_iteration_prev_ind[0] = sample.n-1;
  }

  this->trainfile = StringCopy(trainfile);
  this->modelfile = StringCopy(modelfile);
  this->learn_parm = learn_parm;
  this->kernel_parm = kernel_parm;
  this->sparm = sparm;
}

StructuredSVMOnlineLearner::~StructuredSVMOnlineLearner() {
  if(cached_examples) {
    for(int i = 0; i < t; i++) 
      if(cached_examples[i])
	free_SVM_cached_sample_set(cached_examples[i], m);
    free(cached_examples);
  }

  m->free_struct_sample(sample);
  m->free_struct_model(sm);
  free(modelfile);
  free(trainfile);

  if(examples_by_iteration_number) free(examples_by_iteration_number);
  if(examples_by_iteration_next_ind) free(examples_by_iteration_next_ind);
  if(examples_by_iteration_prev_ind) free(examples_by_iteration_prev_ind);
}




void StructuredSVMOnlineLearner::AddExample(const char *fname) {
  omp_set_lock(&my_lock);
  sample.examples = (EXAMPLE*)realloc(sample.examples, sizeof(EXAMPLE)*(sample.n+1));
  sample.examples[sample.n] = m->read_struct_example(fname, sparm);
  ex_num_iters = (int*)realloc(ex_num_iters, sizeof(int)*(sample.n+1));
  ex_num_iters[sample.n] = 0;

  generalization_errors_by_n = (double*)realloc(generalization_errors_by_n, sizeof(double)*(sample.n+1));
  optimization_errors_by_n = (double*)realloc(optimization_errors_by_n, sizeof(double)*(sample.n+1));
  model_errors_by_n = (double*)realloc(model_errors_by_n, sizeof(double)*(sample.n+1));
  regularization_errors_by_n = (double*)realloc(regularization_errors_by_n, sizeof(double)*(sample.n+1));
  losses_by_n = (double*)realloc(losses_by_n, sizeof(double)*(sample.n+1));
  generalization_errors_by_n[sample.n] = optimization_errors_by_n[sample.n] = losses_by_n[sample.n] =
    model_errors_by_n[sample.n] = regularization_errors_by_n[sample.n] = 0;
  sample.n++;

  
  // Add this example to the front of the list of examples that have been iterated over 0 times
  examples_by_iteration_next_ind = (int*)realloc(examples_by_iteration_next_ind, sizeof(int)*(sample.n+1));
  examples_by_iteration_prev_ind = (int*)realloc(examples_by_iteration_prev_ind, sizeof(int)*(sample.n+1));
  if(examples_by_iteration_number[0] >= 0) {
    examples_by_iteration_prev_ind[sample.n-1] = examples_by_iteration_prev_ind[examples_by_iteration_number[0]];
    examples_by_iteration_next_ind[sample.n-1] = examples_by_iteration_number[0];
    examples_by_iteration_next_ind[examples_by_iteration_prev_ind[examples_by_iteration_number[0]]] = sample.n-1;
    examples_by_iteration_prev_ind[examples_by_iteration_number[0]] = sample.n-1;
  } else {
    examples_by_iteration_prev_ind[sample.n-1] = examples_by_iteration_next_ind[sample.n-1] = sample.n-1;
  }
  examples_by_iteration_number[0] = sample.n-1;
  this->currMinIterByExample = 0;

  omp_unset_lock(&my_lock);
}

const char *StructuredSVMOnlineLearner::SaveModel(const char *fname) {
  omp_set_lock(&my_lock);
  const char *modelfile = fname ? fname : this->modelfile;
  if(fname) { 
    char tmp[1000];
    sprintf(tmp, "%s.trainlist%d", trainfile, sample.n);
    free(trainfile);
    trainfile = StringCopy(tmp);
  }

  /* write structural model */
  char modelfile2[1000];
  free(sm.w);
  SVECTOR *w = smult_s(sum_w, 1.0/(lambda*t));
  sm.w = create_nvector(sm.sizePsi, w);
  bool freeModel = false;
  if(!sm.svm_model) {
    sm.svm_model = (MODEL*)malloc(sizeof(MODEL));
	memset(sm.svm_model, 0, sizeof(MODEL));
    sm.svm_model->kernel_parm = *kernel_parm;
    sm.svm_model->totwords = sm.sizePsi;
    sm.svm_model->totdoc = sample.n;
    freeModel = true;
  }
  fprintf(stderr, "write model %s\n", modelfile);
  m->write_struct_model(modelfile, &sm, sparm);
  free_svector(w);

  sprintf(modelfile2, "%s.learner", modelfile);
  FILE *fout = fopen(modelfile2, "wb");
  if(!fout) { fprintf(stderr, "Failed to open %s for writing\n", modelfile2); return NULL; }
  long tm = GetElapsedTime();
  assert(fwrite(&lambda, sizeof(double), 1, fout) &&
	 fwrite(&eps, sizeof(double), 1, fout) &&
	 fwrite(&R, sizeof(int), 1, fout) &&
	 fwrite(&curr, sizeof(long), 1, fout) &&
	 fwrite(&t, sizeof(long), 1, fout) &&
	 fwrite(&n, sizeof(long), 1, fout) &&
	 fwrite(&minItersBeforeNewExample, sizeof(long), 1, fout) &&
	 fwrite(&M, sizeof(int), 1, fout) &&
	 fwrite(&sum_generalization_error, sizeof(double), 1, fout) &&
	 fwrite(&sum_iter_error, sizeof(double), 1, fout) &&
	 fwrite(&sum_model_error, sizeof(double), 1, fout) &&
	 fwrite(&regularization_error, sizeof(double), 1, fout) &&
	 fwrite(learn_parm, 1, sizeof(LEARN_PARM), fout) &&
	 fwrite(kernel_parm, 1, sizeof(KERNEL_PARM), fout) &&
	 fwrite(&tm, 1, sizeof(long), fout) );

  int len = strlen(modelfile)+1;  
  assert(fwrite(&len, sizeof(int), 1, fout) && fwrite(modelfile, sizeof(char), len, fout));

  len = strlen(trainfile)+1;  
  assert(fwrite(&len, sizeof(int), 1, fout) && fwrite(trainfile, sizeof(char), len, fout));

  write_svector(sum_w, fout);

  assert(fwrite(ex_num_iters, sizeof(int), sample.n, fout));
  len = cached_examples ? t : 0;
  assert(fwrite(&len, sizeof(int), 1, fout));
  if(cached_examples) {
    for(int i = 0; i < t; i++)
      write_SVM_cached_sample_set(cached_examples[i], fout, sparm);
  }
  
  
  if(generalization_errors_by_n) assert(fwrite(generalization_errors_by_n, sizeof(double), n, fout));
  if(optimization_errors_by_n) assert(fwrite(optimization_errors_by_n, sizeof(double), n, fout));
  if(model_errors_by_n) assert(fwrite(model_errors_by_n, sizeof(double), n, fout));
  if(regularization_errors_by_n) assert(fwrite(regularization_errors_by_n, sizeof(double), n, fout));
  if(losses_by_n) assert(fwrite(losses_by_n, sizeof(double), n, fout));

  if(generalization_errors_by_t) assert(fwrite(generalization_errors_by_t, sizeof(double), t, fout));
  if(iter_errors_by_t) assert(fwrite(iter_errors_by_t, sizeof(double), t, fout));
  if(model_errors_by_t) assert(fwrite(model_errors_by_t, sizeof(double), t, fout));
  if(constraint_errors_by_t) assert(fwrite(constraint_errors_by_t, sizeof(double), t, fout));
  if(regularization_errors_by_t) assert(fwrite(regularization_errors_by_t, sizeof(double), t, fout));
  if(losses_by_t) assert(fwrite(losses_by_t, sizeof(double), t, fout));
  if(elapsed_time_by_t) assert(fwrite(elapsed_time_by_t, sizeof(double), t, fout));

  if(iter_examples) assert(fwrite(iter_examples, sizeof(long), t, fout));

  fclose(fout);

  m->save_examples(trainfile, sample);

  if(freeModel) { free(sm.svm_model); sm.svm_model = NULL; }

  omp_unset_lock(&my_lock);

  return modelfile;
}

void StructuredSVMOnlineLearner::ReadModel(const char *fname) {
  const char *modelfile = fname ? fname : this->modelfile;

  /* read structural model */
  char modelfile2[1000];  

  sparm = (STRUCT_LEARN_PARM*)malloc(sizeof(STRUCT_LEARN_PARM));
  memset(sparm, 0, sizeof(STRUCT_LEARN_PARM));
  sm = m->read_struct_model(modelfile, sparm);

  sprintf(modelfile2, "%s.learner", modelfile);
  FILE *fin = fopen(modelfile2, "rb");
  assert(fin);
  learn_parm = (LEARN_PARM*)malloc(sizeof(LEARN_PARM));
  kernel_parm = (KERNEL_PARM*)malloc(sizeof(KERNEL_PARM));
  assert(fread(&lambda, sizeof(double), 1, fin) &&
	 fread(&eps, sizeof(double), 1, fin) &&
	 fread(&R, sizeof(int), 1, fin) &&
	 fread(&curr, sizeof(long), 1, fin) &&
	 fread(&t, sizeof(long), 1, fin) &&
	 fread(&n, sizeof(long), 1, fin) &&
	 fread(&minItersBeforeNewExample, sizeof(long), 1, fin) &&
	 fread(&M, sizeof(int), 1, fin) &&
	 fread(&sum_generalization_error, sizeof(double), 1, fin) &&
	 fread(&sum_iter_error, sizeof(double), 1, fin) &&
	 fread(&sum_model_error, sizeof(double), 1, fin) &&
	 fread(&regularization_error, sizeof(double), 1, fin) &&
	 fread(learn_parm, 1, sizeof(LEARN_PARM), fin) &&
	 fread(kernel_parm, 1, sizeof(KERNEL_PARM), fin) &&
	 fread(&base_time, 1, sizeof(long), fin) );

  kernel_parm->gram_matrix = NULL;
  learn_parm->svm_cost = NULL;

  int len;
  assert(fread(&len, sizeof(int), 1, fin));
  this->modelfile = (char*)malloc(len*sizeof(char)); 
  assert(fread(this->modelfile, sizeof(char), len, fin));

  assert(fread(&len, sizeof(int), 1, fin));
  trainfile = (char*)malloc(len*sizeof(char)); 
  assert(fread(trainfile, sizeof(char), len, fin));

  sample = m->read_struct_examples(trainfile, sparm); 
  m->init_struct_model(sample, &sm, sparm, learn_parm, kernel_parm);

  sum_w = read_svector(fin);

  ex_num_iters = (int*)malloc(sample.n*sizeof(int)); 
  assert(fread(ex_num_iters, sizeof(int), sample.n, fin));
  assert(fread(&len, sizeof(int), 1, fin));
  if(len) {
    assert(len == t);
    cached_examples = (SVM_cached_sample_set**)malloc(t*sizeof(SVM_cached_sample_set*)); 
    for(int i = 0; i < t; i++)
      cached_examples[i] = read_SVM_cached_sample_set(fin, sparm);
  }

  generalization_errors_by_n = (double*)realloc(generalization_errors_by_n, sizeof(double)*(sample.n+1));
  optimization_errors_by_n = (double*)realloc(optimization_errors_by_n, sizeof(double)*(sample.n+1));
  model_errors_by_n = (double*)realloc(model_errors_by_n, sizeof(double)*(sample.n+1));
  regularization_errors_by_n = (double*)realloc(regularization_errors_by_n, sizeof(double)*(sample.n+1));
  losses_by_n = (double*)realloc(losses_by_n, sizeof(double)*(sample.n+1));
  generalization_errors_by_t = (double*)realloc(generalization_errors_by_t, sizeof(double)*t);
  iter_errors_by_t = (double*)realloc(iter_errors_by_t, sizeof(double)*t);
  model_errors_by_t = (double*)realloc(model_errors_by_t, sizeof(double)*t);
  constraint_errors_by_t = (double*)realloc(constraint_errors_by_t, sizeof(double)*t);
  regularization_errors_by_t = (double*)realloc(regularization_errors_by_t, sizeof(double)*t);
  losses_by_t = (double*)realloc(losses_by_t, sizeof(double)*t);
  elapsed_time_by_t = (double*)realloc(elapsed_time_by_t, sizeof(double)*t);
  iter_examples = (long*)realloc(iter_examples, sizeof(long)*t);
  assert(fread(generalization_errors_by_n, sizeof(double), n, fin));
  assert(fread(optimization_errors_by_n, sizeof(double), n, fin));
  assert(fread(model_errors_by_n, sizeof(double), n, fin));
  assert(fread(regularization_errors_by_n, sizeof(double), n, fin));
  assert(fread(losses_by_n, sizeof(double), n, fin));
  assert(fread(generalization_errors_by_t, sizeof(double), t, fin));
  assert(fread(iter_errors_by_t, sizeof(double), t, fin));
  assert(fread(model_errors_by_t, sizeof(double), t, fin));
  assert(fread(constraint_errors_by_t, sizeof(double), t, fin));
  assert(fread(regularization_errors_by_t, sizeof(double), t, fin));
  assert(fread(losses_by_t, sizeof(double), t, fin));
  assert(fread(elapsed_time_by_t, sizeof(double), t, fin));
  assert(fread(iter_examples, sizeof(long), t, fin));

  // Build a queue defining the order examples will be processed
  int i;
  examples_by_iteration_number = (int*)realloc(examples_by_iteration_number, sizeof(int)*(M+2));
  examples_by_iteration_next_ind = (int*)realloc(examples_by_iteration_next_ind, sizeof(int)*(sample.n+1));
  examples_by_iteration_prev_ind = (int*)realloc(examples_by_iteration_prev_ind, sizeof(int)*(sample.n+1));
  for(i = 0; i <= M+1; i++)
    examples_by_iteration_number[i] = -1;
  for(i = sample.n-1; i >= 0; i--) {
    assert(ex_num_iters[i] <= M && ex_num_iters[i] >= 0);
    if(examples_by_iteration_number[ex_num_iters[i]] >= 0) {
      examples_by_iteration_prev_ind[i] = examples_by_iteration_prev_ind[examples_by_iteration_number[ex_num_iters[i]]];
      examples_by_iteration_next_ind[i] = examples_by_iteration_number[ex_num_iters[i]];
      examples_by_iteration_next_ind[examples_by_iteration_prev_ind[examples_by_iteration_number[ex_num_iters[i]]]] = i;
      examples_by_iteration_prev_ind[examples_by_iteration_number[ex_num_iters[i]]] = i;
    } else {
      examples_by_iteration_prev_ind[i] = examples_by_iteration_next_ind[i] = i;
    }
    examples_by_iteration_number[ex_num_iters[i]] = i;
  }

  fclose(fin);
}

// Sanity check: recompute weights from dual parameters.  Could avoid drifting due to numerical precision errors
void StructuredSVMOnlineLearner::RecomputeWeights() {
  omp_set_lock(&my_lock);
  SVECTOR *sum_w_new = create_svector_zero();
  sum_model_error = 0;
  int i, j;
  for(i = 0; i < t; i++) {
    for(j = 0; cached_examples[i] && j < cached_examples[i]->num_samples; j++) {
      SVECTOR *sum_w_new2 = multadd_ss(sum_w_new, cached_examples[i]->samples[j].dpsi, 1, -cached_examples[i]->samples[j].alpha);
      free_svector(sum_w_new);
      sum_w_new = sum_w_new2;
    }
    sum_model_error += constraint_errors_by_t[i];
  }
  assert(sum_sqr_diff_ss(sum_w, sum_w_new) < .001*lambda*t);
  free_svector(sum_w);
  sum_w = sum_w_new;
  omp_unset_lock(&my_lock);
}

double *ComputeWindowAverageArray(double *a, int n, int window, double *a2, float s, int window2=-1) {
  double cumSum = 0;
  int w = 0, i;
  double *retval = (double*)malloc(sizeof(double)*n);
  if(window2 < 0) window2 = window;

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
      if(w < window2) w++;
      else cumSum -= a2[i-window2];
      cumSum += a2[i];
      retval[i] += cumSum/w*s;
    }
  }
  return retval;
} 



void StructuredSVMOnlineLearner::GetStatisticsByExample(int ave, long *nn, double **gen_err_buff, double **opt_err_buff, double **model_err_buff, 
							double **reg_err_buff, double **train_err_buff, double **test_err_buff) {
  omp_set_lock(&my_lock);
  *nn = n;
  if(gen_err_buff) *gen_err_buff = ComputeWindowAverageArray(generalization_errors_by_n, n, ave, regularization_errors_by_n, 1); //optimization_errors_by_n, -1);
  if(opt_err_buff) *opt_err_buff = ComputeWindowAverageArray(optimization_errors_by_n, n, ave, regularization_errors_by_n, 1); //, model_errors_by_n, -1);
  if(model_err_buff) *model_err_buff = ComputeWindowAverageArray(model_errors_by_n, n, ave, regularization_errors_by_n, 1);
  if(reg_err_buff) *reg_err_buff = ComputeWindowAverageArray(regularization_errors_by_n, n, ave, NULL, 0);
  if(test_err_buff) *test_err_buff = ComputeWindowAverageArray(generalization_errors_by_n, n, ave, regularization_errors_by_n, 1);
  if(train_err_buff) *train_err_buff = ComputeWindowAverageArray(model_errors_by_n, n, ave, regularization_errors_by_n, 1);
  omp_unset_lock(&my_lock);
}

void StructuredSVMOnlineLearner::GetStatisticsByIteration(int ave, long *tt, long *tm, double **gen_err_buff, double **opt_err_buff, double **model_err_buff, 
							double **reg_err_buff, double **train_err_buff, double **test_err_buff, double **time_buff) {
  omp_set_lock(&my_lock);
  *tt = t; 
  *tm = GetElapsedTime();
  if(gen_err_buff) *gen_err_buff = ComputeWindowAverageArray(generalization_errors_by_t, t, ave, regularization_errors_by_t, 1); //iter_errors_by_t, -1);
  if(opt_err_buff) *opt_err_buff = ComputeWindowAverageArray(iter_errors_by_t, t, ave, regularization_errors_by_t, 1); //model_errors_by_t, -1);
  if(model_err_buff) *model_err_buff = ComputeWindowAverageArray(model_errors_by_t, t, ave, regularization_errors_by_t, 1); //NULL, 0);
  if(reg_err_buff) *reg_err_buff = ComputeWindowAverageArray(regularization_errors_by_t, t, ave, NULL, 0);
  if(test_err_buff) *test_err_buff = ComputeWindowAverageArray(generalization_errors_by_t, t, ave, regularization_errors_by_t, 1);
  if(train_err_buff) *train_err_buff = ComputeWindowAverageArray(model_errors_by_t, t, ave, regularization_errors_by_t, 1);
  if(time_buff) { *time_buff = (double*)malloc(sizeof(double)*(t+1)); memcpy(*time_buff, elapsed_time_by_t, sizeof(double)*t); }
  omp_unset_lock(&my_lock);
}

double *StructuredSVMOnlineLearner::GetCurrentWeights() {
  omp_set_lock(&my_lock);
  SVECTOR *w = smult_s(sum_w, 1.0/(lambda*t));
  double *ww = create_nvector(sm.sizePsi, w);
  free_svector(w);
  omp_unset_lock(&my_lock);
  return ww;
}

void free_SVM_cached_sample(SVM_cached_sample *s, SVMStructMethod *m) {
  m->free_label(s->ybar);
  free_svector(s->dpsi);
  s->dpsi = NULL;
  memset(&s->ybar, 0, sizeof(LABEL));
}

void read_SVM_cached_sample(SVM_cached_sample *s, FILE *fin, STRUCT_LEARN_PARM *sparm) {
  assert(fread(&s->loss, sizeof(double), 1, fin) && fread(&s->alpha, sizeof(double), 1, fin) && fread(&s->sqr, sizeof(double), 1, fin));
  //m->read_label(&s->ybar, fin, sparm);
  s->dpsi = read_svector(fin);
}

void write_SVM_cached_sample(SVM_cached_sample *s, FILE *fout, STRUCT_LEARN_PARM *sparm) {
  assert(fwrite(&s->loss, sizeof(double), 1, fout) && fwrite(&s->alpha, sizeof(double), 1, fout) && fwrite(&s->sqr, sizeof(double), 1, fout));
  //m->write_label(&s->ybar, fout, sparm);
  write_svector(s->dpsi, fout);
}



SVM_cached_sample_set *new_SVM_cached_sample_set(int i) {
  SVM_cached_sample_set *retval = (SVM_cached_sample_set*)malloc(sizeof(SVM_cached_sample_set));
  retval->samples = NULL;
  retval->num_samples = 0;
  retval->i = i;
  return retval;
}

void free_SVM_cached_sample_set(SVM_cached_sample_set *s, SVMStructMethod *m) {
  for(int i = 0; i < s->num_samples; i++)
    free_SVM_cached_sample(&s->samples[i], m);
  if(s->samples) free(s->samples);
  free(s);
}

SVM_cached_sample_set *read_SVM_cached_sample_set(FILE *fin, STRUCT_LEARN_PARM *sparm) {
  int i;
  assert(fread(&i, sizeof(int), 1, fin));
  SVM_cached_sample_set *s = new_SVM_cached_sample_set(i);
  assert(fread(&s->num_samples, sizeof(int), 1, fin));
  s->samples = (SVM_cached_sample*)realloc(s->samples, sizeof(SVM_cached_sample)*(s->num_samples+1));
  for(int i = 0; i < s->num_samples; i++)
    read_SVM_cached_sample(&s->samples[i], fin, sparm);
  return s;
}

void write_SVM_cached_sample_set(SVM_cached_sample_set *s, FILE *fout, STRUCT_LEARN_PARM *sparm) {
  assert(fwrite(&s->i, sizeof(int), 1, fout));
  assert(fwrite(&s->num_samples, sizeof(int), 1, fout));
  for(int i = 0; i < s->num_samples; i++)
    write_SVM_cached_sample(&s->samples[i], fout, sparm);
}

void SVM_cached_sample_set_add_sample(SVM_cached_sample_set *s, LABEL ybar, SVECTOR *dpsi, double l) {
  s->samples = (SVM_cached_sample*)realloc(s->samples, sizeof(SVM_cached_sample)*(s->num_samples+1));
  s->samples[s->num_samples].ybar = ybar;
  s->samples[s->num_samples].dpsi = dpsi;
  s->samples[s->num_samples].loss = l;
  s->samples[s->num_samples].alpha = 0;
  s->samples[s->num_samples].sqr = sprod_ss(dpsi, dpsi);
  s->num_samples++;
}


// Update weights sum_w when the alpha parameters for sample set 'set' get scaled by 's', and the alpha parameter for sample 'c' gets 
// incremented by 'd'
SVECTOR *SVM_cached_sample_set_update_parameters(SVM_cached_sample *c, SVECTOR *sum_w, double s, double d, double lambda, long t, 
						 SVM_cached_sample_set *set, SVECTOR **w_i, double *sum_alpha, double *L_i) {
  if(s == 1 && d == 0)  // for speed, we don't need to compute anything when nothing changes
    return sum_w;

  // Update alpha
  if(s != 1) 
    for(int j = 0; j < set->num_samples; j++) 
      set->samples[j].alpha *= s;
  c->alpha += d;
  assert(!isnan(c->alpha));

  // Update sum_w, w_i and L_i
  if(w_i) {
    // w_i^q = s*w_i^{q-1} - d*s->dpsi
    // sum_w^q = sum_w^{q-1} + w_i^q - w_i^{q-1}
    SVECTOR *tmp = smult_s(*w_i, s);
    SVECTOR *w_i_new = multadd_ss(tmp, c->dpsi, 1, -d);
    SVECTOR *tmp2 = add_ss(sum_w, w_i_new);
    SVECTOR *sum_w_new = sub_ss(tmp2, *w_i);

    if(L_i) 
      *L_i = s*(*L_i + sprod_ss(sum_w, *w_i)/(lambda*t)) - sprod_ss(sum_w_new, w_i_new)/(lambda*t) + d*c->loss;

    free_svector(tmp2);
    free_svector(tmp);
    free_svector(sum_w);
    free_svector(*w_i);
    *w_i = w_i_new;
    sum_w = sum_w_new;
  } else {
    assert(s == 1 && !L_i);
    SVECTOR *sum_w_new = multadd_ss(sum_w, c->dpsi, 1, -d);
    free_svector(sum_w);
    sum_w = sum_w_new;
  }

  return sum_w;
}

// Optimize alpha parameter for a single sample 's'.  Implements 'Online Dual Update Step' of writeup
SVECTOR *SVM_cached_sample_optimize_dual(SVM_cached_sample *s, SVECTOR *sum_w, double lambda, long t, 
					 SVM_cached_sample_set *set, SVECTOR **w_i, double *sum_alpha, double *L_i) {
  // When there is just a single label, the alpha parameter which maximizes the dual objective can be solved in closed form
  double d = (sprod_ss(sum_w, s->dpsi) + s->loss*(lambda*t)) / my_max(s->sqr,.0000000001);
  d = my_min(1-s->alpha, my_max(-s->alpha,d));
  return SVM_cached_sample_set_update_parameters(s, sum_w, 1, d, lambda, t, set, w_i, sum_alpha, L_i);
}

// Iteratively optimize alpha parameters for a set of samples 's'.  Implements 'Multi-Sample Dual Update Step' of writeup
SVECTOR *SVM_cached_sample_set_optimize_dual(SVM_cached_sample_set *s, SVECTOR *sum_w, double lambda, long t, int R) {
  if(s->num_samples == 1) {  
    // alpha expand
    return SVM_cached_sample_optimize_dual(&s->samples[0], sum_w, lambda, t);
  } else {
    double sum_alpha = 0;  // \sum_ybar \alpha_{i,ybar}
    double L_i = 0;        // L_i = <w,-w_i> + \sum_ybar alpha_{i,ybar} loss(y_i,ybar)
    SVECTOR *w_i;         // w_i = \sum_y alpha_{i,y} (psi(x_i,y_i)-psi(x_i,y))
    int j, r;

    // Initialize w_i, L_i, sum_alpha
    w_i = create_svector_zero();
    for(j = 0; j < s->num_samples; j++) {
      sum_alpha += s->samples[j].alpha;
      SVECTOR *w_i_new = multadd_ss(w_i, s->samples[j].dpsi, 1, -s->samples[j].alpha);
      free_svector(w_i);
      w_i = w_i_new;
      L_i += s->samples[j].alpha * s->samples[j].loss;
    }
    L_i -= sprod_ss(sum_w, w_i)/(lambda*t);
    assert(sum_alpha >= 0 && sum_alpha <= 1);

    for(r = 0; r < R; r++) {
      for(j = 0; j < s->num_samples; j++) {
	if(sum_alpha < 1) {  
	  // alpha expand: solve for the optimal amount 'd' to increase s->samples[j].alpha 
	  // (the value of 'd' that maximizes the increase in the dual objective)
	  sum_w = SVM_cached_sample_optimize_dual(&s->samples[j], sum_w, lambda, t, s, &w_i, &sum_alpha, &L_i);
	} else { 
	  // alpha swap: solve for the optimal amount 'd' to increase s->samples[j].alpha 
	  // while scaling down all s->samples[:].alpha, such that we preserve sum_k{s->samples[:].alpha}=1 
	  // (chose the value of 'd' that maximizes the increase in the dual objective)
	  SVECTOR *v = add_ss(w_i, s->samples[j].dpsi);
	  double e = sprod_ss(sum_w, s->samples[j].dpsi)/(lambda*t) + s->samples[j].loss;
	  double sqr = sprod_ss(v,v);
	  double d = (e-L_i)*(lambda*t) / my_max(sqr,.00000000001);
	  d = my_min(1-s->samples[j].alpha, my_max(-s->samples[j].alpha,s->samples[j].alpha+d));
	  sum_w = SVM_cached_sample_set_update_parameters(&s->samples[j], sum_w, 1-d, d, lambda, t, s, &w_i, &sum_alpha, &L_i);
	  free_svector(v);
	}
      }
    }
  }

  return sum_w;
}