#include "structured_svm.h"

#ifdef WIN32
//#include <windows.h>
//void Sleep(unsigned long  dwMilliseconds);
void usleep(int us) { /*Sleep(us*1000);*/ }
#else
#include <unistd.h>
#endif

bool ReadString(char *str, FILE *fin) {
  int len;
  return fread(&len, sizeof(int), 1, fin) && fread(str, sizeof(char), len, fin);
}

bool WriteString(char *str, FILE *fout) {
  int len=strlen(str)+1;
  return fwrite(&len, sizeof(int), 1, fout) && fwrite(str, sizeof(char), len, fout);
}


void StructuredSVM::Train(const char *modelout, bool saveFull) {
  runForever = saveFull;
  if(!trainset) trainset = new StructuredDataset;

  //C = 100.0; // EYRUN: try a smaller lambda
  lambda = 1.0/C;
  cache_old_examples = method == SPO_DUAL_UPDATE_WITH_CACHE ||
                       method == SPO_DUAL_MULTI_SAMPLE_UPDATE_WITH_CACHE;
  minItersBeforeNewExample = 1;

  if(!examples_by_iteration_number) {
    examples_by_iteration_number = (int*)realloc(examples_by_iteration_number, sizeof(int)*(M+2));
    for(int i = 0; i <= M+1; i++) examples_by_iteration_number[i] = -1;
  }

  if(!sum_w)
    sum_w = new SparseVector;

  //window = trainset->num_examples;  //EYRUN: try to use different window size
  //if (num_thr > trainset->num_examples)
  //	num_thr = 2;
  //	//num_thr = trainset->num_examples;

  /* some training information */
  if(debugLevel > 0) {
    printf("Number of threads=%d\n", num_thr);
    printf("Regularization constant (lambda): %.8lg\n", lambda);
    printf("Approximation factor (epsilon): %.8lg\n", eps);
    printf("Number of training examples (n): %ld\n", (long)trainset->num_examples);
    printf("Feature space dimension (sizePsi): %d\n", sizePsi); fflush(stdout);
  }
  start_time = (long)time(NULL);

  #pragma omp parallel num_threads(num_thr)
  {
    int tid = omp_get_thread_num();
    //printf("tid = %d\n", tid);
    if(tid > 0 || !cache_old_examples) {
      // Worker threads, continuously call ybar=find_most_violated_constraint and then use ybar to update the current weights
      int i = -1;
      while(!finished) {
        // Choose a training example 'i' to process
        Lock();
        SparseVector *w = GetCurrentWeights(false);
        i = ChooseNextExample();
        if(i < 0) {
          Unlock();
          usleep(100000);
          continue;
        }
        StructuredExample *ex = trainset->examples[i];
        Unlock();


        // Find the most violated label ybar = max_y <w_t,psi(x_i,y)>+loss(y_i,y)
        StructuredLabel *ybar = NewStructuredLabel(ex->x);
        VFLOAT score_loss = Inference(ex->x, ybar, w, NULL, ex->y)*featureScale;

	// Compute the gradient dpsi = psi(x_i,ybar)-psi(x_i,y_i) (not including regularization)
        SparseVector ybar_psi = Psi(ex->x, ybar);
        SparseVector y_psi = Psi(ex->x, ex->y);
        SparseVector dpsi = ybar_psi - y_psi;  

        VFLOAT lossval = Loss(ex->y, ybar);
        assert(lossval >= 0);

        OnFinishedIteration(ex->x, ex->y);  // called because the API user might want to free temporary memory caches

        // Use ybar to update the weights
        Lock();
        SVM_cached_sample_set *set = new_SVM_cached_sample_set(i);
        SVM_cached_sample_set_add_sample(set, ybar, dpsi.ptr(), lossval, regularize);
        double e = UpdateWeights(set, -1);
        Unlock();

        if(debugLevel > 2) {
	  //VFLOAT score = score_loss-lossval;      // <w_t,psi(x_i,y)>
	  VFLOAT score = w->dot(ybar_psi)*featureScale;
	  VFLOAT score_gt = w->dot(y_psi)*featureScale;   // <w_t,psi(x_i,y_i)>
	  VFLOAT loss_inference = (VFLOAT)(score_loss-score);
	  //if(abs(loss_inference - (VFLOAT)lossval) > 1)
	  //	printf("Score is off...!");
          printf("Example %d: m=%d slack=%f->%f score=%f score_gt=%f loss=%f loss_inference=%f alpha=%f\n",
                 i, ex_num_iters[i], (VFLOAT)(score_loss-score_gt), (VFLOAT)e, (VFLOAT)(score), (VFLOAT)(score_gt), (VFLOAT)lossval, loss_inference, (VFLOAT)set->samples[0].alpha);
	}

        // Cleanup
        delete w;
        if(!cache_old_examples) 
          free_SVM_cached_sample_set(set);
      }
    } else  {
      // If statement occurs only if method == SPO_DUAL_UPDATE_WITH_CACHE || method == SPO_MULTI_SAMPLE_DUAL_UPDATE_WITH_CACHE
      // Optimization thread, each iteration randomly selects a previous cached example ybar_i (earlier calls
      // to finding the most violated constraint) and then makes an update to w with respect to that sample.  Can be useful
      // if the operation of finding the most violated constraint is much more expensive than the operation of updating w
      int num = 0;
      while(!finished) {
        usleep(1);
        while(t <= num_thr) // wait for at least one first constraint
          usleep(100000);

        // Choose a label from a random iteration (a label ybar from a previous worker thread iteration), and optimize its
        // dual parameters
        Lock();
        int i = rand()%t;
        if(!cached_examples[i]) {
          Unlock();
          continue;
        }
        UpdateWeights(cached_examples[i], i);
        num++;
        Unlock();

        if(num >= t*50) {
          // Sanity check: occasionally recompute weights from dual parameters.  Could help avoid drifting due to numerical precision errors
          Lock();
          RecomputeWeights();
          Unlock();
          num = 0;
        }
      }
    } 
  }

  if(modelout)
    Save(modelout, saveFull);
}

/*
 * The rules for choosing which example to process next are (from highest to lowest precedence):
 *   1) Never allow different threads to process the same example at the same time
 *   2) Prefer processing an example that has been iterated over at least once but less than minItersBeforeNewExample
 *      iterations before going onto a new example (allows better estimation of generalization and optimization errors)
 *   3) Otherwise, prefer processing an example that has been iterated over the fewest times
 */
int StructuredSVM::ChooseNextExample() {
  if(!trainset->num_examples) {
	return -1;
  }

  // Choose the examples in order, selecting the one that has been processed the fewest number of times
  int it = currMinIterByExample;
  while(it <= M && examples_by_iteration_number[it] < 0)
    it++;
  if(hasConverged && it == M) {
    return -1;
  }
  if(it > M || examples_by_iteration_number[it] < 0) {
    return -1;
  }

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


VFLOAT StructuredSVM::Test(const char *testfile, const char *predictionsFile) {
  SparseVector *w = GetCurrentWeights();

  StructuredDataset *testset = LoadDataset(testfile);

  saveBoutFeatures(testset, "test_bout_feat.txt", true, false);

  Lock();
  double sum_los = 0;
  char **strs;
  omp_lock_t l_lock;
  omp_init_lock(&l_lock);

  if(debugLevel > 0) fprintf(stderr, "Evaluating testset...\n");
  if(predictionsFile)
    strs = (char**)malloc(sizeof(char*)*testset->num_examples);

#pragma omp parallel for
  for(int i = 0; i < testset->num_examples; i++) {
    StructuredLabel *y = NewStructuredLabel(testset->examples[i]->x);
    double score = Inference(testset->examples[i]->x, y, w); 
    double los = Loss(testset->examples[i]->y, y);

    if(predictionsFile) {
      double score_gt = w->dot(Psi(testset->examples[i]->x, testset->examples[i]->y)); 
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
    testset->examples[i]->y = y;
    //delete y;

    omp_set_lock(&l_lock);
    sum_los += los;
    omp_unset_lock(&l_lock);
  }
  printf("Average loss was %f\n", (VFLOAT)(sum_los/testset->num_examples));

  saveBoutFeatures(testset, "test_bout_feat_pred.txt", true, false);

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

  delete testset;
  Unlock();

  return (sum_los/testset->num_examples);
}


double StructuredSVM::UpdateWeights(SVM_cached_sample_set *ex, int iterInd) {
  long tt = t;
  double e = sum_w->dot(*ex->samples[0].dpsi)*featureScale/(lambda*(t+1)) + ex->samples[0].loss;
  //printf("<w,x> = %.2f,  loss = %.2f,  e = %.2f\n", sum_w->dot(*ex->samples[0].dpsi)*featureScale/(lambda*(t+1)), ex->samples[0].loss, e);

  if(iterInd == -1) 
    tt = UpdateWeightsAddStatisticsBefore(ex, iterInd, e);

  if(e > 0) {
    switch(method) {
      case SPO_SGD:
      case SPO_SGD_PEGASOS:
        // Stochastic gradient descent: take a step of size -1/(lambda*t) in the direction of psi(ybar,x)-psi(y_i,x)
        assert(iterInd == -1 && ex->num_samples == 1);
        *sum_w -= *ex->samples[0].dpsi * featureScale;  // sum_w = w*(lambda*t)

        if(method == SPO_SGD_PEGASOS) {
          // Project sum_w onto the L2 ball, such that ||w||^2 <= 1/lambda
          double m = sum_w->dot(*sum_w, regularize)/SQR(lambda*t);
          if(m > 1/lambda) 
	    *sum_w *= 1/lambda / m;
        }
        break;

      case SPO_DUAL_UPDATE:
      case SPO_DUAL_UPDATE_WITH_CACHE:
        // Take a step in the direction of psi(ybar,x)-psi(y_i,x), where the chosen step size maximizes the dual objective
        assert(ex->num_samples == 1);
        SVM_cached_sample_optimize_dual(&ex->samples[0], sum_w, lambda, t, featureScale);
	//printf("Updated w\n");
        break;

      case SPO_DUAL_MULTI_SAMPLE_UPDATE:
      case SPO_DUAL_MULTI_SAMPLE_UPDATE_WITH_CACHE:
        // Jointly optimize the dual parameters for a set of different labels within this example
        SVM_cached_sample_set_optimize_dual(ex, sum_w, lambda, t, featureScale, 1, regularize);
        break;

      default:
        fprintf(stderr, "svm_learn_main not implemented: %d\n", method);
        assert(0);
    }
  }

  e = sum_w->dot(*ex->samples[0].dpsi)*featureScale/(lambda*t) + ex->samples[0].loss;
  //printf("<w,x> = %.2f,  loss = %.2f,  e = %.2f,  hasConverged = %d\n", sum_w->dot(*ex->samples[0].dpsi)*featureScale/(lambda*t), ex->samples[0].loss, e, hasConverged);
  UpdateWeightsAddStatisticsAfter(ex, iterInd, e, tt);

  //return my_max(e, 0);
  return e;
}



long StructuredSVM::UpdateWeightsAddStatisticsBefore(SVM_cached_sample_set *ex, int iterInd, double e) {

  // Book-keeping stuff, for estimating generalization error, optimization error, and regret
  ex->samples[0].slack_orig = e;
  if(ex_num_iters[ex->i] == 0) {
    sum_generalization_error += my_max(e,0);
    generalization_errors_by_n[ex->i] = optimization_errors_by_n[ex->i] = model_errors_by_n[ex->i] = my_max(e,0);
    losses_by_n[ex->i] = ex->samples[0].loss;
    n++;
    if(ex->i >= window) sum_generalization_error_window -= generalization_errors_by_n[ex->i-window];
    sum_generalization_error_window += generalization_errors_by_n[ex->i];
  }
  if(t+1 > alloc_t) {
    alloc_t = (int)(alloc_t*1.1)+10;
    iter_examples = (long*)realloc(iter_examples, sizeof(long)*alloc_t);
    iter_errors_by_t = (double*)realloc(iter_errors_by_t, sizeof(double)*alloc_t);
    generalization_errors_by_t = (double*)realloc(generalization_errors_by_t, sizeof(double)*alloc_t);
    constraint_errors_by_t = (double*)realloc(constraint_errors_by_t, sizeof(double)*alloc_t);
    model_errors_by_t = (double*)realloc(model_errors_by_t, sizeof(double)*alloc_t);
    regularization_errors_by_t = (double*)realloc(regularization_errors_by_t, sizeof(double)*alloc_t);
    losses_by_t = (double*)realloc(losses_by_t, sizeof(double)*alloc_t);
    elapsed_time_by_t = (double*)realloc(elapsed_time_by_t, sizeof(double)*alloc_t);
  }
  iter_examples[t] = ex->i;

  sum_iter_error += my_max(e,0);
  iter_errors_by_t[t] = constraint_errors_by_t[t] = optimization_errors_by_n[ex->i] = my_max(e,0);
  generalization_errors_by_t[t] = generalization_errors_by_n[ex->i];
  losses_by_t[t] = ex->samples[0].loss;
  elapsed_time_by_t[t] = (double)GetElapsedTime();

  //printf("t=%d: iter_error_by_t = %.1f, constraint_errors_by_t = %.1f, e = %.1f\n", t, iter_errors_by_t[t], constraint_errors_by_t[t], e);

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

  long tt = t++;
  if(cache_old_examples) {
    cached_examples = (SVM_cached_sample_set**)realloc(cached_examples, (t)*sizeof(SVM_cached_sample_set*));
    cached_examples[tt] = ex;
  }

  return tt;
}

void StructuredSVM::UpdateWeightsAddStatisticsAfter(SVM_cached_sample_set *ex, int iterInd, double e, long tt) {
  ex->samples[0].slack = e;
  if(iterInd == -1) {
    // Book-keeping stuff, for estimating generalization error, optimization error, regret, and convergence

    // Error measured over entire training history
    model_errors_by_n[ex->i] = model_errors_by_t[tt] = constraint_errors_by_t[tt] = my_max(e,0);
    sum_model_error += my_max(e,0);
    regularization_errors_by_t[tt] = regularization_errors_by_n[ex->i] =
      regularization_error = sum_w->dot(*sum_w, regularize)/SQR(lambda*t)*lambda/2;
//    if(debugLevel > 2 || (debugLevel > 1 && tt%10000==9999))
//      printf("t=%d, n=%d: Average Training Error=%f (Model error=%f, Optimization error=%f, Regularization error=%f), Out of sample error=%f\n",
//             (int)t, (int)n, (VFLOAT)(sum_iter_error/t)+(VFLOAT)regularization_error, (VFLOAT)(sum_model_error/t), (VFLOAT)((sum_iter_error-sum_model_error)/t), (VFLOAT)regularization_error, (VFLOAT)(sum_generalization_error/n));


    // Error measured over last set of examples of size window
    int ttt = my_min(tt+1,window), nnn = my_min(n,window);
    //printf("tt=%d: iter_error_by_t = %.1f, constraint_errors_by_t = %.1f, e = %1.f\n", tt, iter_errors_by_t[tt], constraint_errors_by_t[tt], e);
    if(tt >= window) {
      sum_iter_error_window -= iter_errors_by_t[tt-window];
      sum_model_error_window -= constraint_errors_by_t[tt-window];
    } 
    sum_iter_error_window += iter_errors_by_t[tt];
    sum_model_error_window += constraint_errors_by_t[tt];

    double eps_empirical_measured = (sum_iter_error_window-sum_model_error_window) / ttt;
    double eps_generalization_measured = (sum_generalization_error_window) / nnn - sum_model_error_window / ttt;

    //printf("sum_iter_error_window = %.2f, sum_model_error_window = %.2f, eps_empirical_measured = %.2f, eps = %.2f\n", sum_iter_error_window, sum_model_error_window, eps_empirical_measured, eps);

    if(debugLevel > 2 || (debugLevel > 1 && tt%10000==9999))
//      printf("Last %d iters: Average Training Error=%f (Model error=%f, Optimization error=%f, Regularization error=%f), Out of sample error=%f\n",
//             (int)ttt, (float)(sum_iter_error_window/ttt)+(float)regularization_error, (float)(sum_model_error_window/ttt), (float)((sum_iter_error_window-sum_model_error_window)/ttt), (float)regularization_error, (float)(sum_generalization_error_window/nnn));

    if(!hasConverged && eps && eps_empirical_measured < eps && !finished) {
      if(t > window) {
        if(!runForever)
          finished = true;
        if(debugLevel > 0 && !hasConverged) {
          printf("%s at t=%d: epsilon_measured=%f\n", runForever ? "Convergence of empirical error detected" : "Finishing", (int)t, (float)eps_empirical_measured);
          printf("t=%d, n=%d: Average Training Error=%f (Model error=%f, Optimization error=%f, Regularization error=%f), Out of sample error=%f\n",
               (int)t, (int)n, (float)(sum_iter_error/t)+(float)regularization_error, (float)(sum_model_error/t), (float)((sum_iter_error-sum_model_error)/t), (float)regularization_error, (float)(sum_generalization_error/n));
        }
      }
      hasConverged = true;
    } else if(0 && !finished && n > window && eps_generalization_measured < eps && !runForever) {
      printf("Convergence of generalization error detected at t=%d: epsilon_measured=%f\n", (int)n, (float)eps_generalization_measured);
      finished = true;
    }
  } else {
    sum_model_error += my_max(e,0) - constraint_errors_by_t[iterInd];
    sum_model_error_window += my_max(e,0) - constraint_errors_by_t[iterInd];
    constraint_errors_by_t[iterInd] = my_max(e,0);
    losses_by_t[iterInd] = ex->samples[0].loss;
    if(ex_num_iters[ex->i] == M) {
      //model_errors_by_n[ex->i] = my_max(e,0);
      losses_by_n[ex->i] = ex->samples[0].loss;
    }
  }
}

StructuredExample *StructuredSVM::CopyExample(StructuredData *x, StructuredLabel *y) {
  StructuredExample *copy = new StructuredExample();
  copy->x = NewStructuredData();
  copy->y = NewStructuredLabel(copy->x);
  Json::Value xx = x->save(this);
  bool b = copy->x->load(xx, this);
  assert(b);
  Json::Value yy = y->save(this);
  b = copy->y->load(yy, this);
  assert(b);

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
    generalization_errors_by_n = (double*)realloc(generalization_errors_by_n, sizeof(double)*alloc_n);
    optimization_errors_by_n = (double*)realloc(optimization_errors_by_n, sizeof(double)*alloc_n);
    model_errors_by_n = (double*)realloc(model_errors_by_n, sizeof(double)*alloc_n);
    regularization_errors_by_n = (double*)realloc(regularization_errors_by_n, sizeof(double)*alloc_n);
    losses_by_n = (double*)realloc(losses_by_n, sizeof(double)*alloc_n);
    examples_by_iteration_next_ind = (int*)realloc(examples_by_iteration_next_ind, sizeof(int)*alloc_n);
    examples_by_iteration_prev_ind = (int*)realloc(examples_by_iteration_prev_ind, sizeof(int)*alloc_n);
  }
  ex_num_iters[ind] = 0;
  generalization_errors_by_n[ind] = optimization_errors_by_n[ind] = losses_by_n[ind] =
       model_errors_by_n[ind] = regularization_errors_by_n[ind] = 0;
  

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
void StructuredSVM::RecomputeWeights() {
  SparseVector *sum_w_new = new SparseVector;
  sum_iter_error = sum_model_error = sum_iter_error_window = sum_model_error_window = 0;
  int i, j;
  for(i = 0; i < t; i++) {
    for(j = 0; cached_examples[i] && j < cached_examples[i]->num_samples; j++) {
      *sum_w_new += (*cached_examples[i]->samples[j].dpsi*(-cached_examples[i]->samples[j].alpha*featureScale));
    }
    sum_model_error += constraint_errors_by_t[i];
    sum_iter_error += iter_errors_by_t[i];
    if(i >= t-window) {
      sum_iter_error_window += iter_errors_by_t[i];
      sum_model_error_window += constraint_errors_by_t[i];
    }
  }
  //assert(sum_sqr_diff_ss(sum_w, sum_w_new) < .001*lambda*t);
  delete sum_w;
  sum_w = sum_w_new;
}

void StructuredSVM::OptimizeAllConstraints(int num_iter) {
  for(int i = 0; i < num_iter; i++) {
    for(int j = 0; j < t; j++)
      UpdateWeights(cached_examples[j], j);
    RecomputeWeights();
  }
}
void StructuredSVM::SetLambda(double l, int num_iter) {
  Lock();
  lambda = l;
  C = (VFLOAT)(1.0/lambda);
  RecomputeWeights();
  OptimizeAllConstraints(num_iter);
  Unlock();
}

void StructuredSVM::SetFeatureScale(double fs, int num_iter) {
  Lock();
  featureScale = fs;
  RecomputeWeights();
  OptimizeAllConstraints(num_iter);
  Unlock();
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



void StructuredSVM::GetStatisticsByExample(int ave, long *nn, double **gen_err_buff, double **opt_err_buff, double **model_err_buff,
                                                        double **reg_err_buff, double **train_err_buff, double **test_err_buff) {
  Lock();
  *nn = n;
  if(gen_err_buff) *gen_err_buff = ComputeWindowAverageArray(generalization_errors_by_n, n, ave, regularization_errors_by_n, 1); //optimization_errors_by_n, -1);
  if(opt_err_buff) *opt_err_buff = ComputeWindowAverageArray(optimization_errors_by_n, n, ave, regularization_errors_by_n, 1); //, model_errors_by_n, -1);
  if(model_err_buff) *model_err_buff = ComputeWindowAverageArray(model_errors_by_n, n, ave, regularization_errors_by_n, 1);
  if(reg_err_buff) *reg_err_buff = ComputeWindowAverageArray(regularization_errors_by_n, n, ave, NULL, 0);
  if(test_err_buff) *test_err_buff = ComputeWindowAverageArray(generalization_errors_by_n, n, ave, regularization_errors_by_n, 1);
  if(train_err_buff) *train_err_buff = ComputeWindowAverageArray(model_errors_by_n, n, ave, regularization_errors_by_n, 1);
  Unlock();
}

void StructuredSVM::GetStatisticsByIteration(int ave, long *tt, long *tm, double **gen_err_buff, double **opt_err_buff, double **model_err_buff,
                                                        double **reg_err_buff, double **train_err_buff, double **test_err_buff, double **time_buff) {
  Lock();
  *tt = t;
  *tm = GetElapsedTime();
  if(gen_err_buff) *gen_err_buff = ComputeWindowAverageArray(generalization_errors_by_t, t, ave, regularization_errors_by_t, 1); //iter_errors_by_t, -1);
  if(opt_err_buff) *opt_err_buff = ComputeWindowAverageArray(iter_errors_by_t, t, ave, regularization_errors_by_t, 1); //model_errors_by_t, -1);
  if(model_err_buff) *model_err_buff = ComputeWindowAverageArray(model_errors_by_t, t, ave, regularization_errors_by_t, 1); //NULL, 0);
  if(reg_err_buff) *reg_err_buff = ComputeWindowAverageArray(regularization_errors_by_t, t, ave, NULL, 0);
  if(test_err_buff) *test_err_buff = ComputeWindowAverageArray(generalization_errors_by_t, t, ave, regularization_errors_by_t, 1);
  if(train_err_buff) *train_err_buff = ComputeWindowAverageArray(model_errors_by_t, t, ave, regularization_errors_by_t, 1);
  if(time_buff) { *time_buff = (double*)malloc(sizeof(double)*(t+1)); memcpy(*time_buff, elapsed_time_by_t, sizeof(double)*t); }
  Unlock();
}

SparseVector *StructuredSVM::GetCurrentWeights(bool lock) {
  if(lock) Lock();
  SparseVector *retval = sum_w->mult_scalar(1.0/(lambda*t), regularize);
  if(lock) Unlock();
  return retval;
}


void free_SVM_cached_sample(SVM_cached_sample *s) {
  delete s->ybar;
  delete s->dpsi;
  s->ybar = NULL;
  s->dpsi = NULL;
}

void read_SVM_cached_sample(SVM_cached_sample *s, FILE *fin, StructuredSVM *svm, StructuredData *x) {
  bool b = (fread(&s->loss, sizeof(VFLOAT), 1, fin) && fread(&s->alpha, sizeof(VFLOAT), 1, fin) && fread(&s->sqr, sizeof(VFLOAT), 1, fin));
  assert(b);
  char str[100000];
  ReadString(str, fin);
  s->ybar = svm->NewStructuredLabel(x);
  s->dpsi = new SparseVector;
  Json::Reader reader;
  Json::Value v;
  b = reader.parse(str, v);
  assert(b);
  s->ybar->load(v, svm);
  s->dpsi->read(fin);
}


void write_SVM_cached_sample(SVM_cached_sample *s, FILE *fout, StructuredSVM *svm) {
  bool b = (fwrite(&s->loss, sizeof(VFLOAT), 1, fout) && fwrite(&s->alpha, sizeof(VFLOAT), 1, fout) && fwrite(&s->sqr, sizeof(VFLOAT), 1, fout));
  assert(b);
  char str[100000];
  Json::FastWriter writer;
  Json::Value v = s->ybar->save(svm);
  strcpy(str, writer.write(v).c_str());
  WriteString(str, fout);
  s->dpsi->write(fout);
}



SVM_cached_sample_set *new_SVM_cached_sample_set(int i) {
  SVM_cached_sample_set *retval = (SVM_cached_sample_set*)malloc(sizeof(SVM_cached_sample_set));
  retval->samples = NULL;
  retval->num_samples = 0;
  retval->i = i;
  return retval;
}

void free_SVM_cached_sample_set(SVM_cached_sample_set *s) {
  for(int i = 0; i < s->num_samples; i++)
    free_SVM_cached_sample(&s->samples[i]);
  if(s->samples) free(s->samples);
  free(s);
}

SVM_cached_sample_set *read_SVM_cached_sample_set(FILE *fin, StructuredSVM *svm, StructuredData *x) {
  int i;
  bool b = (fread(&i, sizeof(int), 1, fin));
  assert(b);
  SVM_cached_sample_set *s = new_SVM_cached_sample_set(i);
  b = (fread(&s->num_samples, sizeof(int), 1, fin));
  assert(b);
  s->samples = (SVM_cached_sample*)realloc(s->samples, sizeof(SVM_cached_sample)*(s->num_samples+1));
  for(int i = 0; i < s->num_samples; i++)
    read_SVM_cached_sample(&s->samples[i], fin, svm, x);
  return s;
}

void write_SVM_cached_sample_set(SVM_cached_sample_set *s, FILE *fout, StructuredSVM *svm) {
  bool b = (fwrite(&s->i, sizeof(int), 1, fout));
  assert(b);
  b = (fwrite(&s->num_samples, sizeof(int), 1, fout));
  assert(b);
  for(int i = 0; i < s->num_samples; i++)
    write_SVM_cached_sample(&s->samples[i], fout, svm);
}

SVM_cached_sample *SVM_cached_sample_set_add_sample(SVM_cached_sample_set *s, StructuredLabel *ybar, SparseVector *dpsi, VFLOAT l, bool *regularize) {
  s->samples = (SVM_cached_sample*)realloc(s->samples, sizeof(SVM_cached_sample)*(s->num_samples+1));
  SVM_cached_sample *retval = s->samples+s->num_samples;
  retval->ybar = ybar;
  retval->dpsi = dpsi;
  retval->loss = l;
  retval->alpha = 0;
  retval->sqr = dpsi->dot(*dpsi);
  s->num_samples++;

  return retval;
}


// Update weights sum_w when the alpha parameters for sample set 'set' get scaled by 's', and the alpha parameter for sample 'c' gets
// incremented by 'd'
void SVM_cached_sample_set_update_parameters(SVM_cached_sample *c, SparseVector *sum_w, VFLOAT s, VFLOAT d, VFLOAT lambda, long t, VFLOAT featureScale,
					     SVM_cached_sample_set *set, SparseVector *w_i, VFLOAT *sum_alpha, VFLOAT *L_i) {
  if(s == 1 && d == 0)  // for speed, we don't need to compute anything when nothing changes
    return;

  // Update alpha
  if(s != 1)
    for(int j = 0; j < set->num_samples; j++)
      set->samples[j].alpha *= s;
  c->alpha += d;
  //assert(!isnan(c->alpha));

  // Update sum_w, w_i and L_i
  if(w_i) {
    double old_L = 0;
    if(L_i)
      old_L = *L_i + sum_w->dot(*w_i)/(lambda*t);

    // w_i^q = s*w_i^{q-1} - d*s->dpsi
    // sum_w^q = sum_w^{q-1} + w_i^q - w_i^{q-1}
    SparseVector wi_old = *w_i;
    if(s != 1) *w_i *= s;
    if(d) *w_i -= *c->dpsi * (d*featureScale);
    *sum_w += (*w_i) - wi_old;

    if(L_i)
      *L_i = s*old_L - (sum_w->dot(*w_i)/(lambda*t)) + d*c->loss;
  } else {
    assert(s == 1 && !L_i);
    *sum_w -= *c->dpsi * (d*featureScale);
  }
}

// Optimize alpha parameter for a single sample 's'.  Implements 'Online Dual Update Step' of writeup
void SVM_cached_sample_optimize_dual(SVM_cached_sample *s, SparseVector *sum_w, VFLOAT lambda, long t, VFLOAT featureScale,
				     SVM_cached_sample_set *set, SparseVector *w_i, VFLOAT *sum_alpha, VFLOAT *L_i) {
  // When there is just a single label, the alpha parameter which maximizes the dual objective can be solved in closed form
  VFLOAT d = (sum_w->dot(*s->dpsi)*featureScale + s->loss*(lambda*t)) / my_max(s->sqr*SQR(featureScale),.0000000001);
  d = my_min(1-s->alpha, my_max(-s->alpha,d));
  return SVM_cached_sample_set_update_parameters(s, sum_w, 1, d, lambda, t, featureScale, set, w_i, sum_alpha, L_i);
}

// Iteratively optimize alpha parameters for a set of samples 's'.  Implements 'Multi-Sample Dual Update Step' of writeup
void SVM_cached_sample_set_optimize_dual(SVM_cached_sample_set *s, SparseVector *sum_w, VFLOAT lambda, long t, VFLOAT featureScale, int R, bool *regularize) {
  if(s->num_samples == 1) {
    // alpha expand
    return SVM_cached_sample_optimize_dual(&s->samples[0], sum_w, lambda, t, featureScale);
  } else {
    VFLOAT sum_alpha = 0;  // \sum_ybar \alpha_{i,ybar}
    VFLOAT L_i = 0;        // L_i = <w,-w_i> + \sum_ybar alpha_{i,ybar} loss(y_i,ybar)
    SparseVector *w_i;         // w_i = \sum_y alpha_{i,y} (psi(x_i,y_i)-psi(x_i,y))
    int j, r;

    // Initialize w_i, L_i, sum_alpha
    w_i = new SparseVector;
    for(j = 0; j < s->num_samples; j++) {
      sum_alpha += s->samples[j].alpha;
      *w_i -= *s->samples[j].dpsi * (s->samples[j].alpha*featureScale);
      L_i += s->samples[j].alpha * s->samples[j].loss;
    }
    L_i -= sum_w->dot(*w_i)/(lambda*t);
    assert(sum_alpha >= 0 && sum_alpha <= 1);

    for(r = 0; r < R; r++) {
      for(j = 0; j < s->num_samples; j++) {
        if(sum_alpha < 1) {
          // alpha expand: solve for the optimal amount 'd' to increase s->samples[j].alpha
          // (the value of 'd' that maximizes the increase in the dual objective)
          SVM_cached_sample_optimize_dual(&s->samples[j], sum_w, lambda, t, featureScale, s, w_i, &sum_alpha, &L_i);
        } else {
          // alpha swap: solve for the optimal amount 'd' to increase s->samples[j].alpha
          // while scaling down all s->samples[:].alpha, such that we preserve sum_k{s->samples[:].alpha}=1
          // (chose the value of 'd' that maximizes the increase in the dual objective)
          SparseVector v = *w_i + (*s->samples[j].dpsi * featureScale);
          VFLOAT e = sum_w->dot(*s->samples[j].dpsi)*featureScale/(lambda*t) + s->samples[j].loss;
          VFLOAT sqr = v.dot(v, regularize);
          VFLOAT d = (e-L_i)*(lambda*t) / my_max(sqr,.00000000001);
          d = my_min(1-s->samples[j].alpha, my_max(-s->samples[j].alpha,s->samples[j].alpha+d));
          SVM_cached_sample_set_update_parameters(&s->samples[j], sum_w, 1-d, d, lambda, t, featureScale, s, w_i, &sum_alpha, &L_i);
        }
      }
    }
    delete w_i;
  }
}
  




bool StructuredSVM::SaveOnlineData(const char *fname) {
  FILE *fout = fopen(fname, "wb");
  if(!fout) return false;

  long tm = GetElapsedTime();
  bool b = (fwrite(&curr, sizeof(long), 1, fout) &&
         fwrite(&minItersBeforeNewExample, sizeof(long), 1, fout) &&
         fwrite(&M, sizeof(int), 1, fout) &&
         fwrite(&sum_generalization_error, sizeof(double), 1, fout) &&
         fwrite(&sum_iter_error, sizeof(double), 1, fout) &&
         fwrite(&sum_model_error, sizeof(double), 1, fout) &&
         fwrite(&sum_generalization_error_window, sizeof(double), 1, fout) &&
         fwrite(&sum_iter_error_window, sizeof(double), 1, fout) &&
         fwrite(&sum_model_error_window, sizeof(double), 1, fout) &&
         fwrite(&regularization_error, sizeof(double), 1, fout) &&
         fwrite(&tm, 1, sizeof(long), fout) &&
         fwrite(&n, sizeof(int), 1, fout) &&
         fwrite(ex_num_iters, sizeof(int), trainset->num_examples, fout));
  assert(b);

  if(generalization_errors_by_n) (fwrite(generalization_errors_by_n, sizeof(double), n, fout));
  if(optimization_errors_by_n) (fwrite(optimization_errors_by_n, sizeof(double), n, fout));
  if(model_errors_by_n) (fwrite(model_errors_by_n, sizeof(double), n, fout));
  if(regularization_errors_by_n) (fwrite(regularization_errors_by_n, sizeof(double), n, fout));
  if(losses_by_n) (fwrite(losses_by_n, sizeof(double), n, fout));

  if(generalization_errors_by_t) (fwrite(generalization_errors_by_t, sizeof(double), t, fout));
  if(iter_errors_by_t) (fwrite(iter_errors_by_t, sizeof(double), t, fout));
  if(model_errors_by_t) (fwrite(model_errors_by_t, sizeof(double), t, fout));
  if(constraint_errors_by_t) (fwrite(constraint_errors_by_t, sizeof(double), t, fout));
  if(regularization_errors_by_t) (fwrite(regularization_errors_by_t, sizeof(double), t, fout));
  if(losses_by_t) (fwrite(losses_by_t, sizeof(double), t, fout));
  if(elapsed_time_by_t) (fwrite(elapsed_time_by_t, sizeof(double), t, fout));

  if(iter_examples) (fwrite(iter_examples, sizeof(long), t, fout));
  if(examples_by_iteration_number) (fwrite(examples_by_iteration_number, sizeof(int), M+1, fout));
  if(examples_by_iteration_next_ind) (fwrite(examples_by_iteration_next_ind, sizeof(int), trainset->num_examples, fout));
  if(examples_by_iteration_prev_ind) (fwrite(examples_by_iteration_prev_ind, sizeof(int), trainset->num_examples, fout));

  int len = cached_examples ? t : 0;
  b = (fwrite(&len, sizeof(int), 1, fout));
  assert(b);
  if(cached_examples) {
    for(int i = 0; i < t; i++)
      write_SVM_cached_sample_set(cached_examples[i], fout, this);
  }

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
         fread(&sum_generalization_error, sizeof(double), 1, fin) &&
         fread(&sum_iter_error, sizeof(double), 1, fin) &&
         fread(&sum_model_error, sizeof(double), 1, fin) &&
         fread(&sum_generalization_error_window, sizeof(double), 1, fin) &&
         fread(&sum_iter_error_window, sizeof(double), 1, fin) &&
         fread(&sum_model_error_window, sizeof(double), 1, fin) &&
         fread(&regularization_error, sizeof(double), 1, fin) &&
         fread(&base_time, 1, sizeof(long), fin) &&
         fread(&n, sizeof(int), 1, fin) &&
         fread(ex_num_iters, sizeof(int), trainset->num_examples, fin));
  assert(b);

  alloc_t = t;
  generalization_errors_by_n = (double*)realloc(generalization_errors_by_n, sizeof(double)*alloc_n);
  optimization_errors_by_n = (double*)realloc(optimization_errors_by_n, sizeof(double)*alloc_n);
  model_errors_by_n = (double*)realloc(model_errors_by_n, sizeof(double)*alloc_n);
  regularization_errors_by_n = (double*)realloc(regularization_errors_by_n, sizeof(double)*alloc_n);
  losses_by_n = (double*)realloc(losses_by_n, sizeof(double)*alloc_n);
  generalization_errors_by_t = (double*)realloc(generalization_errors_by_t, sizeof(double)*alloc_t);
  iter_errors_by_t = (double*)realloc(iter_errors_by_t, sizeof(double)*alloc_t);
  model_errors_by_t = (double*)realloc(model_errors_by_t, sizeof(double)*alloc_t);
  constraint_errors_by_t = (double*)realloc(constraint_errors_by_t, sizeof(double)*alloc_t);
  regularization_errors_by_t = (double*)realloc(regularization_errors_by_t, sizeof(double)*alloc_t);
  losses_by_t = (double*)realloc(losses_by_t, sizeof(double)*alloc_t);
  elapsed_time_by_t = (double*)realloc(elapsed_time_by_t, sizeof(double)*alloc_t);
  iter_examples = (long*)realloc(iter_examples, sizeof(long)*alloc_t);
  examples_by_iteration_number = (int*)realloc(examples_by_iteration_number, sizeof(int)*(M+2));
  examples_by_iteration_next_ind = (int*)realloc(examples_by_iteration_next_ind, sizeof(int)*alloc_n);
  examples_by_iteration_prev_ind = (int*)realloc(examples_by_iteration_prev_ind, sizeof(int)*alloc_n);
  examples_by_iteration_number[M+1] = -1;

  b = (fread(generalization_errors_by_n, sizeof(double), n, fin)); assert(b);
  b = (fread(optimization_errors_by_n, sizeof(double), n, fin)); assert(b);
  b = (fread(model_errors_by_n, sizeof(double), n, fin)); assert(b);
  b = (fread(regularization_errors_by_n, sizeof(double), n, fin)); assert(b);
  b = (fread(losses_by_n, sizeof(double), n, fin)); assert(b);
  b = (fread(generalization_errors_by_t, sizeof(double), t, fin)); assert(b);
  b = (fread(iter_errors_by_t, sizeof(double), t, fin)); assert(b);
  b = (fread(model_errors_by_t, sizeof(double), t, fin)); assert(b);
  b = (fread(constraint_errors_by_t, sizeof(double), t, fin)); assert(b);
  b = (fread(regularization_errors_by_t, sizeof(double), t, fin)); assert(b);
  b = (fread(losses_by_t, sizeof(double), t, fin)); assert(b);
  b = (fread(elapsed_time_by_t, sizeof(double), t, fin)); assert(b);
  b = (fread(iter_examples, sizeof(long), t, fin)); assert(b);
  b = (fread(examples_by_iteration_number, sizeof(int), M+1, fin)); assert(b);
  b = (fread(examples_by_iteration_next_ind, sizeof(int), trainset->num_examples, fin)); assert(b);
  b = (fread(examples_by_iteration_prev_ind, sizeof(int), trainset->num_examples, fin)); assert(b);


  int len;
  assert(fread(&len, sizeof(int), 1, fin));
  if(len) {
    assert(len == t);
    cached_examples = (SVM_cached_sample_set**)malloc(t*sizeof(SVM_cached_sample_set*));
    for(int i = 0; i < t; i++)
      cached_examples[i] = read_SVM_cached_sample_set(fin, this, trainset->examples[iter_examples[i]]->x);
  }

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
  for(i = 0; i < trainset->num_examples; i++)
    CreateTrainingExampleQueues(i);
}

