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
    StructuredExample *ex = new StructuredExample;
    ex->x = NewStructuredData();
    ex->y = NewStructuredLabel(ex->x);
    Json::Value r;
    if(!reader.parse(line, r) || !r.isMember("x") || !r.isMember("y") || 
       !ex->x->load(r["x"], this) || !ex->y->load(r["y"], this)) { 
      fprintf(stderr, "Error parsing dataset example %s\n", line);
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
  root["Regularization (C)"] = C;
  root["Training accuracy (epsilon)"] = eps;
  root["Feature Scale"] = featureScale;
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

  if(root.isMember("Sum w")) { if(!sum_w) sum_w = new SparseVector; sum_w->load(root["Sum w"]); }
  C = root.get("Regularization (C)", 0).asDouble();
  lambda = 1/C;
  featureScale = root.get("Feature Scale", 1).asDouble();
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
  featureScale = 1;
  debugLevel = 2;
  t = 0;
  sum_w = NULL;
 
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
  minItersBeforeNewExample = 0;

  sum_generalization_error = sum_model_error = sum_iter_error = 0;
  sum_generalization_error_window = sum_model_error_window = sum_iter_error_window = 0;
  window = 1000;

  ex_num_iters = NULL;
  cached_examples = NULL;

  alloc_n = alloc_t = 0;
  generalization_errors_by_n = optimization_errors_by_n = model_errors_by_n = regularization_errors_by_n = losses_by_n = NULL;
  generalization_errors_by_t = iter_errors_by_t = model_errors_by_t = constraint_errors_by_t = regularization_errors_by_t = losses_by_t = elapsed_time_by_t = NULL;
  iter_examples = NULL;

  examples_by_iteration_number = NULL;
  examples_by_iteration_next_ind = NULL;
  examples_by_iteration_prev_ind = NULL;

  trainset = NULL;

  regularize = NULL;
  learnWeights = NULL;
  weightConstraints = NULL;

  num_thr = omp_get_num_procs();
#ifdef MAX_THREADS
  if(num_thr > MAX_THREADS) num_thr = MAX_THREADS;
#endif
  if(num_thr < 2) num_thr = 2;
}

StructuredSVM::~StructuredSVM() {
  if(cached_examples) {
    for(int i = 0; i < t; i++)
      if(cached_examples[i])
        free_SVM_cached_sample_set(cached_examples[i]);
    free(cached_examples);
  }

  if(trainset) delete trainset;
  if(trainfile) free(trainfile);
  if(modelfile) free(modelfile);

  if(regularize)
    free(regularize);
  if(learnWeights)
    free(learnWeights);
  if(weightConstraints)
    free(weightConstraints);

  if(iter_examples) free(iter_examples); 
  if(ex_num_iters) free(ex_num_iters); 
  if(examples_by_iteration_number) free(examples_by_iteration_number);
  if(examples_by_iteration_next_ind) free(examples_by_iteration_next_ind);
  if(examples_by_iteration_prev_ind) free(examples_by_iteration_prev_ind);

  if(generalization_errors_by_n) free(generalization_errors_by_n);
  if(generalization_errors_by_t) free(generalization_errors_by_t);
  if(optimization_errors_by_n) free(optimization_errors_by_n);
  if(iter_errors_by_t) free(iter_errors_by_t);
  if(model_errors_by_n) free(model_errors_by_n);
  if(model_errors_by_t) free(model_errors_by_t);
  if(regularization_errors_by_n) free(regularization_errors_by_n);
  if(regularization_errors_by_t) free(regularization_errors_by_t);
  if(losses_by_n) free(losses_by_n);
  if(losses_by_t) free(losses_by_t);
  if(constraint_errors_by_t) free(constraint_errors_by_t);
  if(elapsed_time_by_t) free(elapsed_time_by_t);

  if(sum_w) delete sum_w;

  omp_destroy_lock(&my_lock);
}



StructuredExample::StructuredExample() { 
  x = NULL; 
  y = NULL; 
}

StructuredExample::~StructuredExample() {
  if(x) delete x;
  if(y) delete y;
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
