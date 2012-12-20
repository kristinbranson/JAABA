#include "structured_svm.h"

StructuredDataset *StructuredSVM::LoadDataset(const char *fname) {
  if(debugLevel > 0) fprintf(stderr, "Reading dataset %s...", fname);
  
  Lock();

  FILE *fin = fopen(fname, "r");
  if(!fin) {
    fprintf(stderr, "Couldn't open dataset file %s\n", fname);
    Unlock();
    return NULL;
  }

  StructuredDataset *d = new StructuredDataset();
  char *line = new char[1000000];
  
  Json::Reader reader;
  while(fgets(line, 9999999, fin) && strlen(line) > 1) {
    chomp(line);
    Json::Value r;
    if(!reader.parse(line, r)) {
      fprintf(stderr, "Error parsing dataset example %s\n", line);
      delete d;
      fclose(fin);
      return false;
    }
    StructuredExample *ex = new StructuredExample;
    ex->x = NewStructuredData();
    ex->y = NewStructuredLabel(ex->x);
    if(r.isMember("y_latent")) ex->y_latent = NewStructuredLabel(ex->x);
    if(!r.isMember("x") || !r.isMember("y") || 
       !ex->x->load(r["x"], this) || !ex->y->load(r["y"], this) ||
       (r.isMember("y_latent") && !ex->y_latent->load(r["y_latent"], this))) { 
      fprintf(stderr, "Error parsing values for dataset example %s\n", line);
      delete ex;
      delete d;
      fclose(fin);
      return false;
    }
    d->AddExample(ex);
  }
  fclose(fin);
  Unlock();
  delete [] line;

  if(debugLevel > 0) fprintf(stderr, "done\n");

  return d;
}



bool StructuredSVM::SaveDataset(StructuredDataset *d, const char *fname, int start_from) {
  if(debugLevel > 0 && start_from == 0) fprintf(stderr, "Saving dataset %s...", fname);

  Lock();

  FILE *fout = fopen(fname, start_from>0 ? "a" : "w");
  if(!fout) {
    fprintf(stderr, "Couldn't open dataset file %s for writing\n", fname);
    Unlock();
    return false;
  }

  Json::FastWriter writer;

  char data[100000];
  for(int i = start_from; i < d->num_examples; i++) {
    Json::Value o;
    o["x"] = d->examples[i]->x->save(this);
    o["y"] = d->examples[i]->y->save(this);
    if(d->examples[i]->y_latent) o["y_latent"] =d->examples[i]->y_latent->save(this);
    strcpy(data, writer.write(o).c_str());
    chomp(data);
    fprintf(fout, "%s\n", data);
  }
  fclose(fout);
  Unlock();

  if(debugLevel > 0 && start_from == 0) fprintf(stderr, "done\n");

  return true;
}

bool StructuredSVM::Save(const char *fname, bool saveFull) {
  Lock();
  Json::Value root;
  if(modelfile) free(modelfile);
  modelfile = StringCopy(fname);
  if(sum_w) root["Sum w"] = sum_w->save();
  root["Sum w scale"] = sum_w_scale;
  if(useWeights) {
    Json::Value a(Json::arrayValue);
    for(int i = 0; i < sizePsi; i++)
      a[i] = useWeights[i];
    root["Use Weights"] = a;
  }
  root["Regularization (C)"] = C;
  root["Training accuracy (epsilon)"] = eps;
  root["T"] = (int)t;
  if(trainfile) root["Training Set"] = trainfile;
  root["Custom"] = Save();
  if(saveFull) {
    char full_name[1000];
    sprintf(full_name, "%s.online", fname);
    root["Online Data"] = full_name;
    if(!SaveOnlineData(full_name)) { Unlock(); return false; }
  }

  Json::StyledWriter writer;
  FILE *fout = fopen(fname, "w");
  if(!fout) { fprintf(stderr, "Couldn't open %s for writing\n", fname); Unlock(); return false; }
  fprintf(fout, "%s", writer.write(root).c_str());
  fclose(fout);
  Unlock();

  return true;
}

bool StructuredSVM::Load(const char *fname, bool loadFull) {
  Lock();
  if(modelfile) free(modelfile);
  modelfile = StringCopy(fname);

  char *str = ReadStringFile(fname);
  if(!str) { fprintf(stderr, "Couldn't open %s for reading\n", fname); Unlock(); return false; }

  if(trainfile) free(trainfile);
  trainfile = NULL;

  Json::Reader reader;
  Json::Value root;
  if(!reader.parse(str, root)) { 
    fprintf(stderr, "Couldn't read JSON file %s\n", fname); Unlock(); return false; 
  }

  if(root.isMember("Sum w")) { 
    if(!sum_w) sum_w = new SparseVector; 
    sum_w->load(root["Sum w"]); 
  }
  if(root.isMember("Sum w scale")) 
    sum_w_scale = root.get("Sum w scale", 0).asDouble();
  if(root.isMember("Use Weights")) {
    Json::Value a = root["Use Weights"];
    useWeights = (bool*)malloc(sizeof(bool)*a.size());
    for(int i = 0; i < (int)a.size(); i++)
      useWeights[i] = a[i].asBool();
  }

  if(root.isMember("Regularization (C)")) {
    C = root.get("Regularization (C)", 0).asDouble();
    lambda = 1/C;
  }
  t = root.get("T", 0).asInt();
  if(root.isMember("Training Set")) { 
    char str[1000]; strcpy(str, root.get("Training Set", "").asString().c_str()); trainfile = StringCopy(str); 
  }
  if(!Load(root["Custom"])) { Unlock(); return false; }
     
  Unlock();

  if(loadFull && !root.isMember("Online Data")) {
    
    fprintf(stderr, "Can't load full data from %s\n", fname); //Unlock(); return false; 
  } else if(loadFull) {
    char str[1000]; strcpy(str, root["Online Data"].asString().c_str());
    if(!LoadOnlineData(str)) { return false; }
  }
  return true;
}



StructuredSVM::StructuredSVM() {
  omp_init_lock(&my_lock);
  
  eps = 0;
  C = 0;
  sizePsi = 0;
  debugLevel = 2;
  t = 0;
  sum_w = NULL;
  u_i_buff = NULL;
  trainLatent = false;
  pauseWorkers = false;

  trainfile = modelfile = NULL; 

  lambda = C ? 1/C : 0;
  method = SPO_DUAL_UPDATE;
  base_time = 0;
  runForever = false;
  hasConverged = false;
  currMinIterByExample = 0;
  curr = 0;
  n = 0;
  M = 0;
  finished = false;
  maxCachedSamplesPerExample = 0;
  keepAllEvictedLabels = false;
  isTesting = false;
  savingCachedExamples = false;

  minItersBeforeNewExample = 1;
  isMultiSample = false;
  numMineHardNegativesRound = 500;
  mineHardNegativesRound = 0;
  numHardNegativesPerExample = 10;
  max_samples = 50;
  maxLoss = 0;
  updateFromCacheThread = false;
  numCacheUpdatesPerIteration = 0;
  numMultiSampleIterations = 20;

  dumpModelStartTime = 0;
  dumpModelFactor = 0;
  numModelDumps = 0;

  sum_iter_error = 0;
  sum_iter_error_window = 0;
  window = 1000;

  ex_num_iters = NULL;
  ex_first_iter = NULL;
  validationfile = NULL;
  last_example = 0;
  maxIters = 1000000000;

  alloc_n = alloc_t = 0;
  generalization_errors_by_n = NULL;
  generalization_errors_by_t = iter_errors_by_t = sum_dual_by_t = regularization_errors_by_t = losses_by_t = elapsed_time_by_t = NULL;
  iter_examples = NULL;
  sum_generalization_error = sum_generalization_error_window = 0;

  examples_by_iteration_number = NULL;
  examples_by_iteration_next_ind = NULL;
  examples_by_iteration_prev_ind = NULL;

  trainset = NULL;

  regularize = NULL;
  learnWeights = NULL;
  useWeights = NULL;
  weightConstraints = NULL;
  regularization_error = 0;
  sum_dual = 0;
  sum_alpha_loss = 0;
  sum_w_sqr = 0;
  sum_w_scale = 1;
  canScaleW = false;
  mergeSamples = true;
  relabelingExample = false;

  runMultiThreaded = 0;
  num_thr = omp_get_num_procs();
#ifdef MAX_THREADS
  if(num_thr > MAX_THREADS) num_thr = MAX_THREADS;
#endif
}

StructuredSVM::~StructuredSVM() {

  if(trainset) delete trainset;
  if(trainfile) free(trainfile);
  if(modelfile) free(modelfile);
  if(validationfile) free(validationfile);

  if(regularize)
    free(regularize);
  if(learnWeights)
    free(learnWeights);
  if(useWeights)
    free(useWeights);
  if(weightConstraints)
    free(weightConstraints);

  if(iter_examples) free(iter_examples); 
  if(ex_num_iters) free(ex_num_iters); 
  if(ex_first_iter) free(ex_first_iter); 
  if(examples_by_iteration_number) free(examples_by_iteration_number);
  if(examples_by_iteration_next_ind) free(examples_by_iteration_next_ind);
  if(examples_by_iteration_prev_ind) free(examples_by_iteration_prev_ind);

  if(generalization_errors_by_n) free(generalization_errors_by_n);
  if(generalization_errors_by_t) free(generalization_errors_by_t);
  if(iter_errors_by_t) free(iter_errors_by_t);
  if(sum_dual_by_t) free(sum_dual_by_t);
  if(regularization_errors_by_t) free(regularization_errors_by_t);
  if(losses_by_t) free(losses_by_t);
  if(elapsed_time_by_t) free(elapsed_time_by_t);

  if(sum_w) delete sum_w;
  if(u_i_buff) free(u_i_buff);

  omp_destroy_lock(&my_lock);
}

void StructuredSVM::InferLatentValues(StructuredDataset *d) {
  SparseVector *w = GetCurrentWeights(false);
  omp_lock_t l_lock;
  omp_init_lock(&l_lock);
  double sum_dual_before = sum_dual;

  if(debugLevel)
    fprintf(stderr, "Inferring latent values...\n");

#pragma omp parallel for
  for(int i = 0; i < d->num_examples; i++) {
    if(!d->examples[i]->y)
      d->examples[i]->y = NewStructuredLabel(d->examples[i]->x);
    if(d->examples[i]->y_latent) {
      Inference(d->examples[i]->x, d->examples[i]->y, w, d->examples[i]->y_latent);

      if(d->examples[i]->set) {
	omp_set_lock(&l_lock);
	SparseVector *psi_gt = Psi(d->examples[i]->x, d->examples[i]->y).ptr();
	SparseVector d_gt = (*psi_gt - *d->examples[i]->set->psi_gt);
	double d_sum_w_sqr = SQR(d->examples[i]->set->alpha)*d_gt.dot(d_gt) + 2*d->examples[i]->set->alpha*sum_w->dot(d_gt);
	*sum_w += d_gt*d->examples[i]->set->alpha;
	sum_w_sqr += d_sum_w_sqr;
	regularization_error = sum_w_sqr/SQR(sum_w_scale)*lambda/2;
	sum_dual -= d_sum_w_sqr/(2*sum_w_scale);
	delete d->examples[i]->set->psi_gt;
	d->examples[i]->set->psi_gt = psi_gt;
	omp_unset_lock(&l_lock);
      }
      OnFinishedIteration(d->examples[i]->x, d->examples[i]->y);
    }
  }
  omp_destroy_lock(&l_lock);
  if(debugLevel)
    fprintf(stderr, "done (%f->%f)\n", (float)sum_dual_before, (float)sum_dual);

  sum_w_sqr = sum_w->dot(*sum_w, regularize);
  regularization_error = sum_w_sqr/SQR(sum_w_scale)*lambda/2;
  sum_dual = -sum_w_sqr/(2*sum_w_scale) + sum_alpha_loss;
}

StructuredExample::StructuredExample() { 
  x = NULL; 
  y = NULL; 
  y_latent = NULL;
  set = NULL;
}

StructuredExample::~StructuredExample() {
  if(x) delete x;
  if(y) delete y;
  if(y_latent) delete y_latent;
  if(set) free_SVM_cached_sample_set(set);
}

StructuredDataset::StructuredDataset() { 
  examples = NULL; 
  num_examples = 0; 
}

StructuredDataset::~StructuredDataset() {
  if(examples) {
    for(int i = 0; i < num_examples; i++)
      delete examples[i];
    free(examples);
  }
}

void StructuredDataset::AddExample(StructuredExample *e) {
  examples = (StructuredExample**)realloc(examples, sizeof(StructuredExample*)*(num_examples+1));
  examples[num_examples++] = e;
}


void StructuredDataset::Randomize() {
  int *perm = RandPerm(num_examples);
  StructuredExample **examples_new = (StructuredExample**)malloc(sizeof(StructuredExample*)*(num_examples));
  for(int i = 0; i < num_examples; i++) 
    examples_new[i] = examples[perm[i]];
  free(examples);
  examples = examples_new;
  free(perm);
}

char *StructuredSVM::VisualizeExample(const char *htmlDir, StructuredExample *ex, const char *extraInfo) {
  Json::StyledWriter writer;
  std::string str = writer.write(ex->y->save(this));
  char *retval = (char*)malloc(strlen(str.c_str())+1);
  strcpy(retval, str.c_str());
  return retval;
}
