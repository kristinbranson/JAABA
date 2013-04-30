#include "svm_fly_behavior_sequence.h"
#include "online_interactive_server.h"
#include "matlab_helper.h"

// TODO: consider adding any jaaba specific network routines to this class
class StructuredLearnerJaabaRpc : public StructuredLearnerRpc {
public:
  StructuredLearnerJaabaRpc(StructuredSVM *l) : StructuredLearnerRpc(l) {}
  bool ClassifyExample(const Json::Value& root, Json::Value& response);
  bool RelabelExample(const Json::Value& root, Json::Value& response);
  bool SetBehaviors(const Json::Value& root, Json::Value& response);
  bool SetFrameFeatures(const Json::Value& root, Json::Value& response);
  bool SetBoutFeatures(const Json::Value& root, Json::Value& response);
protected:
  void AddMethods();
};  

void StructuredLearnerJaabaRpc::AddMethods() {
  StructuredLearnerRpc::AddMethods();

  Json::Value set_behaviors_parameters, set_frame_parameters, set_bout_parameters, set_behaviors_returns;
  set_behaviors_parameters["behaviors"] = "An array of behaviors, with each entry behavior a pair {\"name\":name,\"color\",color}";
  set_behaviors_returns["message"] = "A response indicating success";
  server->RegisterMethod(new JsonRpcMethod<StructuredLearnerJaabaRpc>(this, &StructuredLearnerJaabaRpc::SetBehaviors, "set_behaviors", "Initialize the list of behaviors we are training", set_behaviors_parameters, set_behaviors_returns));

  set_frame_parameters["frame_features"] = "An array of frame feature definitions, see SVMBehaviorSequence::LoadFrameFeatureParams() for the appropriate format";
  server->RegisterMethod(new JsonRpcMethod<StructuredLearnerJaabaRpc>(this, &StructuredLearnerJaabaRpc::SetFrameFeatures, "set_frame_features", "Initialize the list of per-frame features", set_frame_parameters, set_behaviors_returns));

  set_bout_parameters["frame_features"] = "An array of bout feature definitions, see SVMBehaviorSequence::LoadBoutFeatureParams() for the appropriate format";
  server->RegisterMethod(new JsonRpcMethod<StructuredLearnerJaabaRpc>(this, &StructuredLearnerJaabaRpc::SetBoutFeatures, "set_bout_features", "Initialize the list of bout-level features", set_bout_parameters, set_behaviors_returns));
}

// Overrides default StructuredLearnerRpc::ClassifyExample() because we want to add support
// to only do inference only on an interval (t0,t1)
bool StructuredLearnerJaabaRpc::ClassifyExample(const Json::Value& root, Json::Value& response) {
  if(learner) {
    FlyBehaviorBoutFeatures *x = (FlyBehaviorBoutFeatures*)learner->NewStructuredData();
    FlyBehaviorBoutSequence *partial_label = NULL;
    if(!x->load(root["x"], learner)) { delete x; JSON_ERROR("Invalid 'x' parameter", -1); }

    if(root.isMember("partial_label")) {
      partial_label = (FlyBehaviorBoutSequence*)learner->NewStructuredLabel(x);
      if(!partial_label->load(root["partial_label"], learner)) {
	delete x; 
	delete partial_label; 
	JSON_ERROR("Invalid 'partial_label' parameter", -1); 
      }
    } else if(root.isMember("t0") && root.isMember("t1") ) {
      partial_label = (FlyBehaviorBoutSequence*)learner->NewStructuredLabel(x);
      ((SVMFlyBehaviorSequence*)learner)->init_bout_label(partial_label, NULL);
      int t0 = root["t0"].asInt();
      if(t0 > x->GetFirstFrame()) partial_label->AddBout(NONE_BEHAVIOR, 0, t0-x->GetFirstFrame());
      int t1 = root["t1"].asInt();
      if(t1 < x->GetLastFrame()) partial_label->AddBout(NONE_BEHAVIOR, t1-x->GetFirstFrame()+1, x->GetLastFrame()-x->GetFirstFrame()+1);
    }

    FlyBehaviorBoutSequence *y = (FlyBehaviorBoutSequence*)learner->NewStructuredLabel(x);
    SparseVector *w = learner->GetCurrentWeights();
    double score = learner->Inference(x, y, w, partial_label);
    
    if(!isnan(score)) response["score"] = score;
    response["y"] = y->save(learner);

    delete x;
    delete y;
    if(partial_label) delete partial_label;
    delete w;

    return true;
  } else
    return false;
}

// Overrides default StructuredLearnerRpc::ClassifyExample() because we want to add support
// for recomputing bout-level mean, viariance, and median statistics
bool StructuredLearnerJaabaRpc::RelabelExample(const Json::Value& root, Json::Value& response) {
  if(learner) {
    int ind = root.get("index",-1).asInt();
    if(ind < 0 || ind >= learner->GetTrainset()->num_examples) {
      JSON_ERROR("Missing or invalid training example 'index' parameter to relabel_example()", -1); 
    }
    StructuredExample *ex = learner->GetTrainset()->examples[ind];

    StructuredLabel *y = learner->NewStructuredLabel(ex->x);
    if(root.isMember("y")) {
      if(!y->load(root["y"], learner)) { 
	delete y;
	JSON_ERROR("Invalid 'y' parameter", -1); 
      }
    } else {
      delete y;
      JSON_ERROR("No 'y' parameter specified", -1); 
    }
    response["index"] = ind;

    bool recomputeStatistics = root.get("recompute_bout_statistics", false).asBool();
    int num_optimization_iters = root.get("num_optimization_iters", 0).asInt();

    if(recomputeStatistics || num_optimization_iters) {
      learner->PauseWorkerThreads(true, true);
      if(recomputeStatistics) learner->FlushCache();
    }
    learner->RelabelExample(ex, y);
    learner->SaveTrainingSet(ind);

    if(recomputeStatistics) {
      ((SVMBehaviorSequence*)learner)->compute_feature_mean_variance_median_statistics(learner->GetTrainset());
    }
    if(num_optimization_iters) {
      learner->Lock();
      learner->OptimizeAllConstraints(num_optimization_iters);
      learner->Unlock();
    }
    if(recomputeStatistics || num_optimization_iters) {
      learner->PauseWorkerThreads(false);
    }
  }

  return true;
}



bool StructuredLearnerJaabaRpc::SetFrameFeatures(const Json::Value& root, Json::Value& response) {
  SVMFlyBehaviorSequence *l = ((SVMFlyBehaviorSequence*)learner);
  if(root.isMember("frame_features")) {
    l->LoadFrameFeatureParams(root["frame_features"]);
    l->Init(l->GetBehaviors(), l->GetFrameFeatureDefs(), l->NumFrameFeatures(), l->GetBoutFeatureDefs(), l->NumBoutFeatures());

    char str[1000];
    sprintf(str, "Initialized %d frame features", l->NumFrameFeatures());
    response["message"] = str;
  } else {
    JSON_ERROR("Missing 'frame_features' parameter in set_frame_features()", -1); 
  }
  return true;
}


bool StructuredLearnerJaabaRpc::SetBoutFeatures(const Json::Value& root, Json::Value& response) {
  SVMFlyBehaviorSequence *l = ((SVMFlyBehaviorSequence*)learner);
  if(root.isMember("bout_features")) {
    l->LoadBoutFeatureParams(root["frame_features"]);
    l->Init(l->GetBehaviors(), l->GetFrameFeatureDefs(), l->NumFrameFeatures(), l->GetBoutFeatureDefs(), l->NumBoutFeatures());

    char str[1000];
    sprintf(str, "Initialized %d bout features", l->NumBoutFeatures());
    response["message"] = str;
  } else {
    JSON_ERROR("Missing 'bout_features' parameter in set_bout_features()", -1); 
  }
  return true;
}


bool StructuredLearnerJaabaRpc::SetBehaviors(const Json::Value& root, Json::Value& response) {
  SVMFlyBehaviorSequence *l = ((SVMFlyBehaviorSequence*)learner);
  if(root.isMember("behaviors")) {
    if(!l->LoadBehaviorDefinitions(root["behaviors"])) {
      JSON_ERROR("No behavior 'none' found in set_behaviors()", -1); 
    }
    l->Init(l->GetBehaviors(), l->GetFrameFeatureDefs(), l->NumFrameFeatures(), l->GetBoutFeatureDefs(), l->NumBoutFeatures());

    char str[1000];
    sprintf(str, "Initialized %d behaviors", l->NumBehaviors());
    response["message"] = str;
  } else {
    JSON_ERROR("Missing 'behaviors' parameter in set_behaviors()", -1); 
  }
  return true;
}

#ifndef AS_LIBRARY
int main(int argc, const char **argv) {
  char pname[1000], debug_name[1000];
  strcpy(pname, "");
  strcpy(debug_name, "");

  int i = 1;
  while(i<argc) {
    if(argv[i][0] == '-') {
      switch ((argv[i])[1]) {
      case 'B': i++; sprintf(pname, "%s", argv[i]); break;
      case 'D': i++; sprintf(debug_name, "%s", argv[i]); break;
      default: break; 
      }
    }
    i++;
  }
  
  //assert(strlen(feat_name) && behaviors);
  SVMFlyBehaviorSequence *svm = new SVMFlyBehaviorSequence();
  if(strlen(debug_name)) svm->SetDebugDir(debug_name);
  if(strlen(pname)) ((StructuredSVM*)svm)->Load((const char*)pname, false);

  StructuredLearnerJaabaRpc v(svm);
  v.main(argc, argv);
}
#endif


/*  
 * Read the trx features from file and store per frame feature and timestamps 
 */
bool FlyBehaviorBoutFeatures::load(const char *pname, SVMBehaviorSequence *svm) {
  char fname[1001], dirname[10001], folder[1000];
  const char *ptr;
  MATFile *pmat;
  
  /*if(!(ptr=strstr(pname, "trx.mat"))) {
    ExtractPathname(pname, folder);
    sprintf(name, "%s/trx.mat", folder);
  } else
    strcpy(name, pname);
  */

  pmat = matOpen(pname, "r");
  if (pmat == NULL) {
    fprintf(stderr, "Error opening mat file %s\n", pname);
    return false;
  }
  
  
  int i;
  mxArray *trx, *tmp, *tmp2;
  double firstframe, lastframe, fps;
  char sex[1001];
  double fly_id;
  int id = this->fly_id-1;
  if(id < 0) id = 0;
  MAT_GET_VARIABLE(trx, "trx");
  MAT_GET_DOUBLE_FIELD(trx, fly_id, "id");
  this->fly_id = (int)fly_id;
  id = this->fly_id-1;
  if(id < 0) { fprintf(stderr, "Invalid fly_id %d\n", id); id = 0; }
  MAT_GET_DOUBLE_FIELD(trx, firstframe, "firstframe");
  MAT_GET_DOUBLE_FIELD(trx, lastframe, "endframe");
  MAT_GET_DOUBLE_FIELD(trx, fps, "fps");
  MAT_GET_STRING_FIELD(trx, this->matname, "moviefile");
  MAT_GET_STRING_FIELD(trx, this->moviename, "moviename");
  MAT_GET_STRING_FIELD(trx, sex, "sex");

  this->sex = *sex;
  this->fps = fps;
  
  this->firstframe = (int)firstframe;
  this->lastframe = (int)lastframe;
  num_frames = (int)lastframe-(int)firstframe+1;
  int num_base_features = svm->NumFrameFeatures();
  FrameFeature *frame_features = svm->GetFrameFeatureDefs();
  AllocateBuffers(svm, false);

  if(mxGetField(trx, id, "timestamps")) {
    MAT_GET_DOUBLE_ARRAY_FIELD(trx, frame_times, "timestamps", num_frames, firstframe);
  } else {
    for(int t = 0; t < num_frames; t++)
      frame_times[t] = (firstframe-1+t)/fps;
  }
  matClose(pmat);
  mxDestroyArray(trx);
  
  ExtractPathname(pname, dirname);
  strcat(dirname, "/perframe"); 
  for(i = 0; i < num_base_features; i++) {
    sprintf(fname, "%s/%s.mat", dirname, frame_features[i].name);
    MAT_READ_DATA_FILE(features[i], fname, num_frames, firstframe);
  }

  return true;
}

int BehaviorBout_start_frame_cmp(const void *a, const void *b) {
  return ( ((BehaviorBout*)a)->start_frame - ((BehaviorBout*)b)->start_frame );
}

/*
 * Read a behavior label from file
 */
bool FlyBehaviorBoutSequence::load(const char *pname) {
  FlyBehaviorBoutFeatures *xx = (FlyBehaviorBoutFeatures*)x;
  this->firstframe = xx->firstframe;
  this->lastframe = xx->lastframe;

  strcpy(this->labelname, pname);
  
  unsigned int sz = (sizeof(int)+sizeof(BehaviorBout*));
  this->bouts = (BehaviorBout*)malloc(sz);
  memset(this->bouts, 0, sz);
  this->num_bouts = 0;
  this->score = this->loss = this->slack = 0;

  MATFile *pmat;
  pmat = matOpen(pname, "r");
  if (pmat == NULL) {
    fprintf(stderr, "Error opening mat file %s\n", pname);
    return false;
  }
  int len = 0, len2 = 0;
  double *t0 = NULL, *t1 = NULL;
  mxArray *names, *names2, *tmp, *tmp2;
  char name[1001];

  this->fly_id = xx->fly_id;
  int id = this->fly_id-1;
  if(id < 0) { fprintf(stderr, "Invalid fly_id %d\n", id); id = 0; }

  tmp2 = matGetVariable(pmat, "t0s"); 
  if(tmp2) tmp = mxGetCell(tmp2, id);  
  len = mxGetM(tmp)*mxGetN(tmp);	       
  t0 = (double*)malloc(sizeof(double)*len);   
  memcpy(t0, ((double*)mxGetData(tmp)), len*sizeof(double)); 
  mxDestroyArray(tmp2);


  MAT_GET_CELL_DOUBLE_ARRAY_VARIABLE(t0, "t0s", len, id);
  MAT_GET_CELL_DOUBLE_ARRAY_VARIABLE(t1, "t1s", len2, id);
  MAT_GET_VARIABLE(names, "names"); 
  names2 = names;//mxGetCell(names, id);
  assert(len == len2);

 
  Behaviors *behaviors = svm->GetBehaviors();
  this->bouts = (BehaviorBout*)malloc(sizeof(BehaviorBout)*(len));
  this->num_bouts = len;
  for(int i = 0; i < len; i++) {
    MAT_GET_CELL_STRING(names2, name, i, "names");
    int j = 0;
    for(j = 0; j < behaviors->num_values; j++) {
      if(!strcasecmp(name, behaviors->values[j].name)) {
	break;
      }
    }
    if(j >= behaviors->num_values) {
      fprintf(stderr, "Could not find behavior %s in %s.%d\n", name, pname, i);
      return false;
    }
    this->bouts[i].start_frame = (int)t0[i]-this->firstframe;
    this->bouts[i].end_frame = (int)t1[i]-this->firstframe;
    this->bouts[i].behavior = j;
  }
  mxDestroyArray(names);
  free(t0);
  free(t1);
  matClose(pmat);
  qsort(this->bouts, len, sizeof(BehaviorBout), BehaviorBout_start_frame_cmp);

  return true;
}

bool FlyBehaviorBoutSequence::save(const char *pname) {
  return false;
}

// Supports reading Kristin's format, which keeps each perframe feature in a separate .mat file
bool FlyBehaviorBoutFeatures::load(const Json::Value &r, StructuredSVM *s) {
  this->fly_id = r.get("fly_id", 0).asInt();  
  return BehaviorBoutFeatures::load(r, s);
}

Json::Value FlyBehaviorBoutFeatures::save(StructuredSVM *s) {
  Json::Value r = BehaviorBoutFeatures::save(s);
  r["fly_id"] = fly_id;
  return r;
}

// Supports reading two formats: 1) Kristin's matlab .trx format, 2) JSON format
bool FlyBehaviorBoutSequence::load(const Json::Value &r, StructuredSVM *s) {
  FlyBehaviorBoutFeatures *xx = (FlyBehaviorBoutFeatures*)x;
  this->firstframe = xx->firstframe;
  this->lastframe = xx->lastframe;
  this->fly_id = r.get("fly_id", 0).asInt();

  return BehaviorBoutSequence::load(r, s);
}

Json::Value FlyBehaviorBoutSequence::save(StructuredSVM *s) {
  Json::Value r = BehaviorBoutSequence::save(s);
  r["fly_id"] = fly_id;

  if(s->IsTesting()) {
    r["fname"] = fname;
    save(fname);
  }

  return r;
}




/*************** The stuff below here is currently identical to svm_fly_behavior_sequence.cpp ************************/


FlyBehaviorBoutFeatures::FlyBehaviorBoutFeatures() : BehaviorBoutFeatures() {
  version = 0;
  strcpy(moviename, "");
  strcpy(matname, "");
  nflies = firstframe = lastframe = 0;
  sex = 0;
  fps = 0;
  fly_id = 0;
}

FlyBehaviorBoutSequence::FlyBehaviorBoutSequence(BehaviorBoutFeatures *x, SVMBehaviorSequence *svm) 
  : BehaviorBoutSequence(x, svm) {
  strcpy(labelname, "");
  strcpy(moviename, "");
  strcpy(matname, "");
  strcpy(trxname, "");
  nflies = firstframe = lastframe = 0;
  is_labeled = 0;
  fly_id = 0;
  this->svm = svm;
}


StructuredLabel *SVMFlyBehaviorSequence::NewStructuredLabel(StructuredData *x) { 
  return new FlyBehaviorBoutSequence((FlyBehaviorBoutFeatures*)x, this); 
}

StructuredData *SVMFlyBehaviorSequence::NewStructuredData() { 
  return new FlyBehaviorBoutFeatures; 
}


SVMFlyBehaviorSequence::SVMFlyBehaviorSequence() : SVMBehaviorSequence()
{
}



// Read a file containing a list of training examples and return a string array of file names
void SVMFlyBehaviorSequence::save_examples(const char *fname, StructuredDataset *dataset) {
  char folder[1000];
  FILE *fout = fopen(fname, "w");
  strcpy(folder, fname);
  for(int i = strlen(folder); i >= 0 && folder[i] != '/' && folder[i] != '\\'; i--) 
    folder[i] = '\0';

  if(!fout) return;

  for(int i = 0; i < dataset->num_examples; i++) {
    char fname2[1000];
    strcpy(fname2, ((FlyBehaviorBoutSequence*)dataset->examples[i]->y)->labelname);
    int j = 0;
    while(fname2[j] == folder[j] || ((fname2[j]=='/'||fname2[j]=='\\')&&(folder[j]=='/'||folder[j]=='\\'))) j++;
    if(fname2[j] == '/' || fname2[j] == '\\')
      j++;
    fprintf(fout, "%s\n", fname2+j);
  }
  fclose(fout);
}

char **SVMFlyBehaviorSequence::load_examples(const char *fname, int *num) {
  char **retval = NULL, line[1000], folder[1000];
  FILE *fin = fopen(fname, "r");
  strcpy(folder, fname);
  for(int i = strlen(folder); i >= 0 && folder[i] != '/' && folder[i] != '\\'; i--) 
    folder[i] = '\0';

  *num = 0;

  if(!fin) return NULL;
  
  while(fgets(line, 999, fin)) {
    chomp(line);
    retval = (char**)realloc(retval, sizeof(char*)*(*num+1));
    retval[*num] = (char*)malloc(strlen(folder)+strlen(line)+2);
    sprintf(retval[*num], "%s%s", folder, line);
    *num = *num+1;
  }
  fclose(fin);

  return retval;
}

