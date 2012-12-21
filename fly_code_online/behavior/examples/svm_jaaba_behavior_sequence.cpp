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
protected:
  void AddMethods();
};  

void StructuredLearnerJaabaRpc::AddMethods() {
  StructuredLearnerRpc::AddMethods();

  Json::Value set_behaviors_parameters, set_behaviors_returns;
  set_behaviors_parameters["behaviors"] = "An array of behaviors, with each entry behavior a pair {\"name\":name,\"color\",color}";
  set_behaviors_returns["message"] = "A response indicating success";
  server->RegisterMethod(new JsonRpcMethod<StructuredLearnerJaabaRpc>(this, &StructuredLearnerJaabaRpc::SetBehaviors, "set_behaviors", "Initialize the list of behaviors we are training", set_behaviors_parameters, set_behaviors_returns));
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
      if(t0 > x->GetFirstFrame()) partial_label->AddBout(0, NONE_BEHAVIOR, 0, t0-x->GetFirstFrame());
      int t1 = root["t1"].asInt();
      if(t1 < x->GetLastFrame()) partial_label->AddBout(0, NONE_BEHAVIOR, t1-x->GetFirstFrame()+1, x->GetLastFrame()-x->GetFirstFrame()+1);
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
      learner->Lock();
      ((SVMBehaviorSequence*)learner)->compute_feature_mean_variance_median_statistics(learner->GetTrainset());
      learner->Unlock();
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



bool StructuredLearnerJaabaRpc::SetBehaviors(const Json::Value& root, Json::Value& response) {
  int noneInd = -1;
  SVMFlyBehaviorSequence *l = ((SVMFlyBehaviorSequence*)learner);

  if(root.isMember("behaviors")) {
    BehaviorGroups *beh = (BehaviorGroups*)malloc(sizeof(BehaviorGroups)); 
    beh->num = 1;
    beh->behaviors = (BehaviorGroup*)malloc(sizeof(BehaviorGroup));
    beh->behaviors[0].num_values = root["behaviors"].size();
    beh->behaviors[0].values = (BehaviorValue*)malloc(beh->behaviors[0].num_values*sizeof(BehaviorValue)); 
    memset(beh->behaviors[0].values, 0, beh->behaviors[0].num_values*sizeof(BehaviorValue));
    for(int i = 0; i < (int)root["behaviors"].size(); i++) {
      strcpy(beh->behaviors[0].values[i].name, root["behaviors"][i].get("name","").asString().c_str());
      strncpy(beh->behaviors[0].values[i].abbreviation, beh->behaviors[0].values[i].name, 3);
      beh->behaviors[0].values[i].abbreviation[3] = '\0';
      beh->behaviors[0].values[i].color = root["behaviors"][i].get("color",rand() & 0x00FFFFFF).asInt();
      if(!strcasecmp(beh->behaviors[0].values[i].name, "none"))
	noneInd = i;
    }
    if(noneInd < 0) {
      JSON_ERROR("No behavior 'none' found in set_behaviors()", -1); 
    } else if(noneInd > 0) {
      // make sure the "none" behavior is at index 0
      BehaviorValue tmp = beh->behaviors[0].values[0];
      beh->behaviors[0].values[0] = beh->behaviors[0].values[noneInd]; 
      beh->behaviors[0].values[noneInd] = tmp;
    }
    l->SetBehaviors(beh);
    l->Init(l->NumBaseFeatures(), beh, 0, l->BaseFeatures());

    char str[1000];
    sprintf(str, "Initialized %d behaviors", beh->behaviors[0].num_values);
    response["message"] = str;
  } else {
    JSON_ERROR("Missing 'behaviors' parameter in set_behaviors()", -1); 
  }
  return true;
}


int main(int argc, const char **argv) {
  char bname[1000], feat_name[1000], debug_name[1000];
  BehaviorGroups *behaviors = NULL;
  strcpy(bname, "");
  strcpy(feat_name, "");
  strcpy(debug_name, "");

  for(int i=1; i<argc; i++) {
    if(argv[i][0] == '-') {
      switch ((argv[i])[1]) {
      case 'B': i++; sprintf(bname, "%s", argv[i]); behaviors = load_behaviors(bname, "."); break;
      case 'F': i++; sprintf(feat_name, "%s", argv[i]); break;
      case 'D': i++; sprintf(debug_name, "%s", argv[i]); break;
      default: i++;
      }
    }
  }
  
  //assert(strlen(feat_name) && behaviors);
  SVMFlyBehaviorSequence *svm = new SVMFlyBehaviorSequence(feat_name, behaviors, -1);
  if(strlen(debug_name)) svm->SetDebugDir(debug_name);
  
  /*
  FlyBehaviorBoutFeatures *f = (FlyBehaviorBoutFeatures*)svm->NewStructuredData();
  f->load("/home/sbranson/code/JAABA/data/movie1_seq1/trx.mat", svm);
  FlyBehaviorBoutSequence *y = (FlyBehaviorBoutSequence*)svm->NewStructuredLabel(f);
  y->load("/home/sbranson/code/JAABA/data/movie1_seq1/labeled_gt_lunge.mat"); 
  delete f;
  delete svm;
  exit(0);
  */

  StructuredLearnerJaabaRpc v(svm);
  v.main(argc, argv);
}



/*  
 * Read the trx features from file and store per frame feature and timestamps 
 */
bool FlyBehaviorBoutFeatures::load(const char *pname, SVMBehaviorSequence *svm) {
  char fname[1001], dirname[10001];
  MATFile *pmat;
  pmat = matOpen(pname, "r");
  if (pmat == NULL) {
    fprintf(stderr, "Error opening mat file %s\n", pname);
    return false;
  }
  
  int id = this->fly_id-1;
  int i;
  mxArray *trx, *tmp, *tmp2;
  double firstframe, lastframe, fly_id, fps;
  char sex[1001];
  MAT_GET_VARIABLE(trx, "trx");
  MAT_GET_DOUBLE_FIELD(trx, firstframe, "firstframe");
  MAT_GET_DOUBLE_FIELD(trx, lastframe, "endframe");
  MAT_GET_DOUBLE_FIELD(trx, fly_id, "id");
  MAT_GET_DOUBLE_FIELD(trx, fps, "fps");
  MAT_GET_STRING_FIELD(trx, this->matname, "moviefile");
  MAT_GET_STRING_FIELD(trx, this->moviename, "moviename");
  MAT_GET_STRING_FIELD(trx, sex, "sex");

  this->sex = *sex;
  this->fps = fps;
  
  this->firstframe = (int)firstframe;
  this->lastframe = (int)lastframe;
  num_frames = (int)lastframe-(int)firstframe+1;
  num_base_features = svm->NumBaseFeatures();
  AllocateBuffers(svm, false);

  MAT_GET_DOUBLE_ARRAY_FIELD(trx, frame_times, "timestamps", num_frames);
  matClose(pmat);
  mxDestroyArray(trx);
  
  ExtractPathname(pname, dirname);
  strcat(dirname, "/perframe"); 
  for(i = 0; i < num_base_features; i++) {
    sprintf(fname, "%s/%s.mat", dirname, ((SVMFlyBehaviorSequence*)svm)->base_feature_names[i]);
    MAT_READ_DATA_FILE(features[i], fname, num_frames);
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

  this->behaviors = behaviors;
  strcpy(this->labelname, pname);
  
  unsigned int sz = (sizeof(int)+sizeof(BehaviorBout*)+2*sizeof(double))*behaviors->num;
  this->bouts = (BehaviorBout**)malloc(sz);
  memset(this->bouts, 0, sz);
  this->num_bouts = (int*)(this->bouts+behaviors->num);
  this->scores = (double*)(this->num_bouts+behaviors->num);
  this->losses = (double*)(this->scores+behaviors->num);
  this->behaviors = behaviors;
  this->score = this->loss = this->slack = 0;

  MATFile *pmat;
  pmat = matOpen(pname, "r");
  if (pmat == NULL) {
    fprintf(stderr, "Error opening mat file %s\n", pname);
    return false;
  }
  int len = 0, len2 = 0, ind = 0;
  double *t0 = NULL, *t1 = NULL;
  mxArray *names, *names2, *tmp, *tmp2;
  char name[1001];
  int id = this->fly_id-1;

  MAT_GET_CELL_DOUBLE_ARRAY_VARIABLE(t0, "t0s", len, id);
  MAT_GET_CELL_DOUBLE_ARRAY_VARIABLE(t1, "t1s", len2, id);
  MAT_GET_VARIABLE(names, "names"); 
  names2 = mxGetCell(names, id);
  assert(len == len2);

  this->bouts[ind] = (BehaviorBout*)malloc(sizeof(BehaviorBout)*(len));
  this->num_bouts[ind] = len;
  for(int i = 0; i < len; i++) {
    MAT_GET_CELL_STRING(names2, name, i, "names");
    int j = 0;
    for(j = 0; j < behaviors->behaviors[ind].num_values; j++) {
      if(!strcasecmp(name, behaviors->behaviors[ind].values[j].name)) {
	break;
      }
    }
    if(j >= behaviors->behaviors[ind].num_values) {
      fprintf(stderr, "Could not find behavior %s in %s.%d\n", name, pname, i);
      return false;
    }
    this->bouts[ind][i].start_frame = (int)t0[i]-this->firstframe;
    this->bouts[ind][i].end_frame = (int)t1[i]-this->firstframe;
    this->bouts[ind][i].behavior = j;
  }
  mxDestroyArray(names);
  free(t0);
  free(t1);
  matClose(pmat);
  qsort(this->bouts[ind], len, sizeof(BehaviorBout), BehaviorBout_start_frame_cmp);

  return true;
}

bool FlyBehaviorBoutSequence::save(const char *pname) {
  return false;
}


bool FlyBehaviorBoutFeatures::load(const Json::Value &r, StructuredSVM *s) {
  this->fly_id = r.get("fly_id", 0).asInt();
  return BehaviorBoutFeatures::load(r, s);
}

Json::Value FlyBehaviorBoutFeatures::save(StructuredSVM *s) {
  Json::Value r = BehaviorBoutFeatures::save(s);
  r["fly_id"] = fly_id;
  return r;
}

bool FlyBehaviorBoutSequence::load(const Json::Value &r, StructuredSVM *s) {
  FlyBehaviorBoutFeatures *xx = (FlyBehaviorBoutFeatures*)x;
  this->firstframe = xx->firstframe;
  this->lastframe = xx->lastframe;
  this->fly_id = r.get("fly_id", 0).asInt();

  if(r.isMember("fname")) {
    char fname[1000];
    strcpy(fname, r.get("fname", "").asString().c_str());
    bool ret = load(fname);
    //remove(fname);
    return ret;
  } else {
    unsigned int sz = (sizeof(int)+sizeof(BehaviorBout*)+2*sizeof(double))*behaviors->num;
    this->bouts = (BehaviorBout**)malloc(sz);
    memset(this->bouts, 0, sz);
    this->num_bouts = (int*)(this->bouts+behaviors->num);
    this->scores = (double*)(this->num_bouts+behaviors->num);
    this->losses = (double*)(this->scores+behaviors->num);
    this->behaviors = behaviors;
    this->score = this->loss = this->slack = 0;

    if((r["num_bouts"].isArray() && (int)r["bouts"].size() == behaviors->num &&
	(int)r["num_bouts"].size() == behaviors->num) || (behaviors->num==1 && !r["num_bouts"].isArray())) {
      ((SVMBehaviorSequence*)s)->init_bout_label(this, NULL);
      
      bool isA = r["num_bouts"].isArray();
      for(int i = 0; i < (isA ? (int)r["num_bouts"].size() : 1); i++) {
	num_bouts[i] = isA ? r["num_bouts"][i].asInt() : r["num_bouts"].asInt();
	this->bouts[i] = (BehaviorBout*)malloc(sizeof(BehaviorBout)*num_bouts[i]);
	Json::Value a = isA ? r["bouts"][i] : r["bouts"];
	for(int j = 0; j < (a.isArray() ? (int)a.size() : 1); j++) {
	  Json::Value c = a.isArray() ? a[j] : a;
	  bouts[i][j].start_frame = c["start_frame"].asInt()-firstframe;
	  bouts[i][j].end_frame = c["end_frame"].asInt()-firstframe;
	  if(c["behavior"].type() == Json::stringValue) {
	    char bname[1001]; strcpy(bname, c["behavior"].asString().c_str());
	    for(bouts[i][j].behavior = 0; bouts[i][j].behavior < behaviors->behaviors[i].num_values; bouts[i][j].behavior++)
	      if(!strcasecmp(bname,  behaviors->behaviors[i].values[bouts[i][j].behavior].name))
		break;
	    if(bouts[i][j].behavior == behaviors->behaviors[i].num_values)
	      return false;  // behavior not found
	  } else {
	    bouts[i][j].behavior = c["behavior"].asInt();
	  }
	}
      }
    }
  }
  return true;
}

Json::Value FlyBehaviorBoutSequence::save(StructuredSVM *s) {
  Json::Value r;
  r["fly_id"] = fly_id;
  if(num_bouts) {
    Json::Value a(Json::arrayValue);
    for(int i = 0; i < behaviors->num; i++)
      a[i] = num_bouts[i];
    r["num_bouts"] = a;
  }

  if(x) {
    FlyBehaviorBoutFeatures *xx = (FlyBehaviorBoutFeatures*)x;
    this->firstframe = xx->firstframe;
    this->lastframe = xx->lastframe;
  }

  if(bouts) {
    Json::Value a(Json::arrayValue);
    for(int i = 0; i < behaviors->num; i++) {
      Json::Value b(Json::arrayValue);
      for(int j = 0; j < num_bouts[i]; j++) {
        Json::Value c;
        c["start_frame"] = bouts[i][j].start_frame+firstframe;
        c["end_frame"] = bouts[i][j].end_frame+firstframe;
        c["behavior"] = behaviors->behaviors[i].values[bouts[i][j].behavior].name;
        c["bout_score"] = bouts[i][j].bout_score;
        b[j] = c;
      }
      a[i] = b;
    }
    r["bouts"] = a;
  }

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
}


StructuredLabel *SVMFlyBehaviorSequence::NewStructuredLabel(StructuredData *x) { return new FlyBehaviorBoutSequence((FlyBehaviorBoutFeatures*)x, this); }
StructuredData *SVMFlyBehaviorSequence::NewStructuredData() { return new FlyBehaviorBoutFeatures; }


SVMFlyBehaviorSequence::SVMFlyBehaviorSequence(const char *feature_params, struct _BehaviorGroups *behaviors, int beh) :
  SVMBehaviorSequence(behaviors, beh)
{
  assert(!behaviors || behaviors->num);
  SVMFeatureParams fparams[MAX_BASE_FEATURES];
  int num_feat = feature_params && strlen(feature_params) ? ReadFeatureParams(feature_params, fparams) : 0;
 Init(num_feat, behaviors, beh, fparams);
}



bool SVMFlyBehaviorSequence::ReadFeatureParam(FILE *modelfl, SVMFeatureParams *p) {
        int num;
        int b[30];
        if((num=fscanf(modelfl, FORMAT__BOUT_FEATURE_PARAMS_READ, &p->feature_sample_smoothness_window, &p->num_temporal_levels, &p->num_bout_max_thresholds,
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

int SVMFlyBehaviorSequence::ReadFeatureParams(const char *fname, SVMFeatureParams *p) {
  int num = 0;
  FILE *fin = fopen(fname, "r");
  if(!fin) { fprintf(stderr, "Couldn't open feature file %s\n", fname); }
  assert(fin);
  while(fscanf(fin, "%s ", feature_defs[num].name) && strlen(feature_defs[num].name) && 
	feature_defs[num].name[strlen(feature_defs[num].name)-1] == ':') {
 	feature_defs[num].name[strlen(feature_defs[num].name)-1] = '\0';

	strcpy(base_feature_names[num], feature_defs[num].name);

	int n = ReadFeatureParam(fin, &p[num++]);   assert(n);
  }
  assert(num);
  fclose(fin);

  return num;
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


const char *SVMFlyBehaviorSequence::get_base_feature_name(int ind) {
  return feature_defs[ind].name;
}


