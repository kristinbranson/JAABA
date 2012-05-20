#include "svm_fly_behavior_sequence.h"
#include "online_interactive_server.h"

char **load_train_list(const char *fname, int *num);

StructuredLabel *SVMFlyBehaviorSequence::NewStructuredLabel(StructuredData *x) { return new FlyBehaviorBoutSequence((FlyBehaviorBoutFeatures*)x, this); }
StructuredData *SVMFlyBehaviorSequence::NewStructuredData() { return new FlyBehaviorBoutFeatures; }


SVMFlyBehaviorSequence::SVMFlyBehaviorSequence(const char *feature_params, struct _BehaviorGroups *behaviors, int beh) :
  SVMBehaviorSequence(behaviors, beh)
{
  assert(behaviors->num);
  SVMFeatureParams fparams[MAX_BASE_FEATURES];
  int num_feat = ReadFeatureParams(feature_params, fparams);
  Init(num_feat, behaviors, beh, fparams);
}


bool SVMFlyBehaviorSequence::read_features(FILE *fin,  FlyBehaviorFeatures *f, double *frames, int num_frames) {
  int len;

  if(!fread(&len, sizeof(int), 1, fin))
    return false;
  assert(len < 1000 && len > 0 && (int)fread(f->name, sizeof(char), len, fin)==len); f->name[len] = '\0';
  assert(fread(&len, sizeof(int), 1, fin) && len < 1000 && len > 0 && (int)fread(f->units_numer, sizeof(char), len, fin)==len); 
  f->units_numer[len] = '\0';
  assert(fread(&len, sizeof(int), 1, fin) && len < 1000 && len > 0 && (int)fread(f->units_denom, sizeof(char), len, fin)==len); 
  f->units_denom[len] = '\0';

  if(frames && num_frames)
    assert(fread(frames, sizeof(double), num_frames, fin));
 
  return true;
}




/*************** Overloaded virtual functions ************************/

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
  assert(fin);
  while(fscanf(fin, "%s ", feature_defs[num].name) && strlen(feature_defs[num].name) && 
	feature_defs[num].name[strlen(feature_defs[num].name)-1] == ':') {
 	feature_defs[num].name[strlen(feature_defs[num].name)-1] = '\0';

	strcpy(base_feature_names[num], feature_defs[num].name);

	assert(ReadFeatureParam(fin, &p[num++]));
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

/*
 * Read the trx features from file and store per frame feature and timestamps 
 */
bool FlyBehaviorBoutFeatures::load(const char *fname, SVMBehaviorSequence *svm, BehaviorBoutSequence *y) {
  FlyBehaviorBoutSequence *fly = (FlyBehaviorBoutSequence*)y;
  int i;

  FlyBehaviorFeatures d;
  double version;
  int len, fly_id, nflies, nfields, firstframe, lastframe;
  char sex; 
  double fps;
  char movie[1001], mat[1001], fname2[1000];

  strcpy(fname2, fname);
  StripFileExtension(fname2);
  strcat(fname2, ".trx");
  FILE *fin = fopen(fname2, "rb");
  assert(fin);
  assert(fread(&version, sizeof(double), 1, fin));
  assert(fread(&len, sizeof(int), 1, fin) && len < 1000 && len > 0 && 
	 (int)fread(movie, sizeof(char), len, fin)==len);
  movie[len] = '\0';
  assert(fread(&len, sizeof(int), 1, fin) && len < 1000 && len > 0 && 
	 (int)fread(mat, sizeof(char), len, fin)==len);
  mat[len] = '\0';
  assert(fread(&nflies, sizeof(int), 1, fin));
  for(i = 0; i < nflies; i++) {
    if(!fly)
      assert(fread(&fly_id, sizeof(int), 1, fin));
    else if(fly->is_labeled)
      assert(fread(&fly_id, sizeof(int), 1, fin) && fly_id == fly->fly_ids[i]);
    else
      assert(fread(&fly->fly_ids[i], sizeof(int), 1, fin));
  }
  assert(fread(&firstframe, sizeof(int), 1, fin) &&
	 fread(&lastframe, sizeof(int), 1, fin) &&
	 fread(&sex, sizeof(char), 1, fin) && fread(&fps, sizeof(double), 1, fin) &&
	 fread(&nfields, sizeof(int), 1, fin) && nfields < MAX_BASE_FEATURES+1);

  if(fly) {
    fly->sex = sex;
    fly->fps = fps;
    if(fly->is_labeled) {
      assert(!strcmp(movie, fly->moviename) && !strcmp(mat, fly->matname) && nflies == fly->nflies && 
	   firstframe == fly->firstframe && lastframe == fly->lastframe);
    } else {
      fly->version = version;
      strcpy(fly->moviename, movie); strcpy(fly->matname, mat); fly->nflies = nflies;
      fly->firstframe = firstframe; fly->lastframe = lastframe;
    }
  }

  num_frames = lastframe-firstframe+1;
  num_base_features = nfields-1;
  assert(num_base_features == svm->NumBaseFeatures());

  AllocateBuffers(svm, false);
  double *buff = (double*)malloc(sizeof(double)*num_frames);
  for(i = 0; i < nfields; i++) {
    assert(((SVMFlyBehaviorSequence*)svm)->read_features(fin, &d, buff, num_frames));
    if(i == 0) assert(!strcmp(d.name, "timestamp"));
    if(!strcmp(d.name, "timestamp")) {
      memcpy(frame_times, buff, sizeof(double)*num_frames);
    } else {
      assert(!strcmp(d.name, ((SVMFlyBehaviorSequence*)svm)->GetFeatureDefs()[i-1].name));
      memcpy(features[i-1], buff, sizeof(double)*num_frames);
      // use the log of certain features: vel, accel, ang_vel, ang_accel, axis_ratio_vel, vel_to_other, accel_to_other, vel_ang_between, accel_ang_between
//       if(strcmp(d.name,"vel") | strcmp(d.name,"accel") | strcmp(d.name,"ang_vel") | strcmp(d.name,"ang_accel") | strcmp(d.name,"axis_ratio_vel") | strcmp(d.name,"vel_to_other") | strcmp(d.name,"accel_to_other") | strcmp(d.name,"vel_ang_between") | strcmp(d.name,"accel_ang_between")) {
// 	for(int f=0; f<num_frames; f++)
// 		features[i-1][f] = log(abs(features[i-1][f])+0.00001);
//       }
    }
  }
  free(buff);

  ComputeCaches(svm);

  fclose(fin);

  return true;
}

FlyBehaviorBoutSequence::FlyBehaviorBoutSequence(BehaviorBoutFeatures *x, SVMBehaviorSequence *svm) 
  : BehaviorBoutSequence(x, svm) {
  version = 0;
  strcpy(labelname, "");
  strcpy(moviename, "");
  strcpy(matname, "");
  strcpy(trxname, "");
  nflies = firstframe = lastframe = 0;
  sex = 0;
  fps = 0;
  is_labeled = 0;
}

/*
 * Read a behavior label from file
 */
bool FlyBehaviorBoutSequence::load(const char *fname) {
  this->behaviors = behaviors;
  strcpy(this->labelname, fname);
  
  unsigned int sz = (sizeof(int)+sizeof(BehaviorBout*)+2*sizeof(double))*behaviors->num;
  this->bouts = (BehaviorBout**)malloc(sz);
  memset(this->bouts, 0, sz);
  this->num_bouts = (int*)(this->bouts+behaviors->num);
  this->scores = (double*)(this->num_bouts+behaviors->num);
  this->losses = (double*)(this->scores+behaviors->num);
  this->behaviors = behaviors;
  this->score = this->loss = this->slack = 0;

  assert(strstr(fname, ".label"));
  
  this->is_labeled = true;

  char behname[1001];
  int len, nbehaviors, i, start_frame, end_frame, num, behavior;
  FILE *fin = fopen(fname, "rb");
  assert(fin);
  assert(fread(&this->version, sizeof(double), 1, fin));
  assert(fread(&len, sizeof(int), 1, fin) && len < 1000 && len > 0 && 
	 (int)fread(this->moviename, sizeof(char), len, fin)==len);
  this->moviename[len] = '\0';
  assert(fread(&len, sizeof(int), 1, fin) && len < 1000 && len > 0 && 
	 (int)fread(this->matname, sizeof(char), len, fin)==len);
  this->matname[len] = '\0';
  assert(fread(&len, sizeof(int), 1, fin) && len < 1000 && len > 0 && 
	 (int)fread(this->trxname, sizeof(char), len, fin)==len);
  this->trxname[len] = '\0';
  assert(fread(&this->nflies, sizeof(int), 1, fin) && this->nflies > 0 && this->nflies < 100);
  for(int i = 0; i < this->nflies; i++)
    assert(fread(&this->fly_ids[i], sizeof(int), 1, fin));
  
  assert(fread(&this->firstframe, sizeof(int), 1, fin) && this->firstframe >= 0 &&
	 fread(&this->lastframe, sizeof(int), 1, fin) && this->lastframe >= this->firstframe);
  assert(fread(&nbehaviors, sizeof(int), 1, fin) && nbehaviors+ADD_DUMMY_BEHAVIORS == behaviors->behaviors[0].num_values);

  
  int ind = 0; // assume for now all behaviors are mutually exclusive
  
  for(i = 0; i < nbehaviors; i++) {
    assert(fread(&len, sizeof(int), 1, fin) && len < 1000 && len > 0 && 
	   (int)fread(behname, sizeof(char), len, fin)==len);
    behname[len] = '\0';
    assert(!strcmp(behname, behaviors->behaviors[ind].values[i+ADD_DUMMY_BEHAVIORS].name));
  }
    
  assert(fread(&num, sizeof(int), 1, fin) && this->num_bouts[ind] >= 0);
  this->bouts[ind] = (BehaviorBout*)malloc(sizeof(BehaviorBout)*(num));
  this->num_bouts[ind] = 0;
  int num_added = 0, lastframe = 0;
  for(i = 0; i < num; i++) {
    assert(fread(&start_frame, sizeof(int), 1, fin));
    assert(fread(&end_frame, sizeof(int), 1, fin) && end_frame >= start_frame);
    assert(fread(&behavior, sizeof(int), 1, fin)); 
    assert(behavior >= 0 && behavior < behaviors->behaviors[ind].num_values); 
    assert(start_frame >= lastframe);
    if(start_frame > lastframe && i > 0) {
      // Pad unlabelled bouts as the "Unknown" class
      this->bouts[ind] = (BehaviorBout*)realloc(this->bouts[ind], sizeof(BehaviorBout)*(num+num_added+1));
      this->bouts[ind][this->num_bouts[ind]].start_frame = lastframe;
      this->bouts[ind][this->num_bouts[ind]].end_frame = start_frame;
      this->bouts[ind][this->num_bouts[ind]++].behavior = 0;
      num_added++;
    }
      
    this->bouts[ind][this->num_bouts[ind]].start_frame = start_frame;
    this->bouts[ind][this->num_bouts[ind]].end_frame = lastframe = end_frame;
    this->bouts[ind][this->num_bouts[ind]].behavior = behavior;
    this->num_bouts[ind]++;
  }
  // EYRUN: since my bouts always start on 0, this is not needed
  // for(i = 0; i < this->num_bouts[ind]; i++) {
  //   this->bouts[ind][i].start_frame -= this->firstframe;
  //   this->bouts[ind][i].end_frame -= this->firstframe;
  // }
  fclose(fin);

  return true;
}


bool FlyBehaviorBoutSequence::save(const char *fname) {
  int len, nbehaviors, i, start_frame, end_frame;

  if(!strlen(moviename)) {
    char path[1000];
    version = .1;
    ExtractPathname(fname, path);
    int len = strlen(path);  if(len) len++;
    strcpy(trxname, fname+len); 
    strcpy(moviename, fname+len); StripFileExtension(moviename); 
    strcpy(labelname, fname); StripFileExtension(labelname); strcat(labelname, ".label");
    strcpy(matname, fname); StripFileExtension(matname); strcat(matname, ".mat");
    nflies = 1;
    fly_ids[0] = 1;
  }

  FILE *fout = fopen(fname, "wb");
  assert(fout);
  assert(fwrite(&this->version, sizeof(double), 1, fout));
  len = strlen(this->moviename);
  assert(fwrite(&len, sizeof(int), 1, fout) && len < 1000 && 
	 (int)fwrite(this->moviename, sizeof(char), len, fout)==len);
  len = strlen(this->matname);
  assert(fwrite(&len, sizeof(int), 1, fout) && len < 1000 && 
	 (int)fwrite(this->matname, sizeof(char), len, fout)==len);
  len = strlen(this->trxname);
  assert(fwrite(&len, sizeof(int), 1, fout) && len < 1000 && 
	 (int)fwrite(this->trxname, sizeof(char), len, fout)==len);
  assert(fwrite(&this->nflies, sizeof(int), 1, fout) && this->nflies < 100);
  for(int i = 0; i < this->nflies; i++)
    assert(fwrite(&this->fly_ids[i], sizeof(int), 1, fout));
  
  assert(fwrite(&this->firstframe, sizeof(int), 1, fout) && this->firstframe >= 0 &&
	 fwrite(&this->lastframe, sizeof(int), 1, fout) && this->lastframe >= this->firstframe);
  nbehaviors = behaviors->behaviors[0].num_values-ADD_DUMMY_BEHAVIORS;
  assert(fwrite(&nbehaviors, sizeof(int), 1, fout));
  
  int ind = 0; // assume for now all behaviors are mutually exclusive
  
  for(i = 0; i < nbehaviors; i++) {
    len = strlen(behaviors->behaviors[ind].values[i+ADD_DUMMY_BEHAVIORS].name);
    assert(fwrite(&len, sizeof(int), 1, fout) && len < 1000 && len > 0 && 
	   (int)fwrite(behaviors->behaviors[ind].values[i+ADD_DUMMY_BEHAVIORS].name, sizeof(char), len, fout)==len);
  }
    
  assert(fwrite(&this->num_bouts[ind], sizeof(int), 1, fout) && this->num_bouts[ind] >= 0);
  for(i = 0; i < this->num_bouts[ind]; i++) {
    //start_frame = this->firstframe + this->bouts[ind][i].start_frame;   // EYRUN
    //end_frame = this->firstframe + this->bouts[ind][i].end_frame;	  
    start_frame = this->bouts[ind][i].start_frame;
    end_frame = this->bouts[ind][i].end_frame;
    assert(fwrite(&start_frame, sizeof(int), 1, fout));
    assert(fwrite(&end_frame, sizeof(int), 1, fout) && end_frame >= start_frame);
    assert(fwrite(&this->bouts[ind][i].behavior, sizeof(int), 1, fout)); 
  }
  fclose(fout);

  return true;
}

int main(int argc, const char **argv) {
  char bname[1000], feat_name[1000];
  BehaviorGroups *behaviors = NULL;
  strcpy(bname, "");
  strcpy(feat_name, "");

  for(int i=1; i<argc; i++) {
    if(argv[i][0] == '-') {
      switch ((argv[i])[1]) {
      case 'B': i++; sprintf(bname, "%s", argv[i]); behaviors = load_behaviors(bname, "."); break;
      case 'F': i++; sprintf(feat_name, "%s", argv[i]); break;
      default: i++;
      }
    }
  }

  assert(strlen(feat_name) && behaviors);
  StructuredLearnerRpc v(new SVMFlyBehaviorSequence(feat_name, behaviors, -1));
  v.main(argc, argv);
}
