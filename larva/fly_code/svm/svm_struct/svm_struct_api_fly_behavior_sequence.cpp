#include "svm_struct_api_fly_behavior_sequence.h"

char **load_train_list(const char *fname, int *num);


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


int SVMFlyBehaviorSequence::ReadFeatureParams(const char *fname, SVMFeatureParams *p) {
  int num = 0;
  FILE *fin = fopen(fname, "r");
  assert(fin);
  while(fscanf(fin, "%s ", feature_defs[num].name) && strlen(feature_defs[num].name) && 
	feature_defs[num].name[strlen(feature_defs[num].name)-1] == ':') {
    feature_defs[num].name[strlen(feature_defs[num].name)-1] = '\0';
    assert(ReadFeatureParam(fin, &p[num++]));
  }
  assert(num);
  fclose(fin);

  return num;
}

// Return the number of frames in a FlyBehaviorBoutSequence
int SVMFlyBehaviorSequence::num_frames(void *b) {
  return ((FlyBehaviorBoutSequence*)b)->lastframe - ((FlyBehaviorBoutSequence*)b)->firstframe + 1;
}

// Read a file containing a list of training examples and return a string array of file names
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

/*
 * Convert a BehaviorBoutSequence struct in y into a custom fly label format
 */
void SVMFlyBehaviorSequence::load_from_bout_sequence(BehaviorBoutSequence *y, void *b) {
  int beh;
  FlyBehaviorBoutSequence *fly = (FlyBehaviorBoutSequence*)b;

  for(beh = 0; beh < behaviors->num; beh++) {
    fly->bouts->num_bouts[beh] = y->num_bouts[beh];
    fly->bouts->bouts[beh] = (BehaviorBout*)realloc(fly->bouts->bouts[beh], sizeof(BehaviorBout)*y->num_bouts[beh]);
    memcpy(fly->bouts->bouts[beh], y->bouts[beh], sizeof(BehaviorBout)*y->num_bouts[beh]);
  }
}

/*
 * Convert fly labels into a BehaviorBoutSequence struct
 */
BehaviorBoutSequence *SVMFlyBehaviorSequence::create_behavior_bout_sequence(void *b, BehaviorGroups *behaviors, bool build_partial_label) {
  FlyBehaviorBoutSequence *fly = (FlyBehaviorBoutSequence*)b;
  unsigned int sz = sizeof(BehaviorBoutSequence)+(sizeof(int)+sizeof(BehaviorBout*)+2*sizeof(double))*behaviors->num;
  BehaviorBoutSequence *bouts = (BehaviorBoutSequence*)my_malloc(sz);
  memset(bouts, 0, sz);
  bouts->bouts = (BehaviorBout**)(bouts+1);
  bouts->num_bouts = (int*)(bouts->bouts+behaviors->num);
  bouts->scores = (double*)(bouts->num_bouts+behaviors->num);
  bouts->losses = (double*)(bouts->scores+behaviors->num);
  bouts->behaviors = behaviors;
  bouts->score = bouts->loss = bouts->slack = 0;
  for(int i = 0; i < behaviors->num; i++) {
    bouts->num_bouts[i] = fly->bouts->num_bouts[i];
    bouts->bouts[i] = (BehaviorBout*)malloc(sizeof(BehaviorBout)*fly->bouts->num_bouts[i]);
    memcpy(bouts->bouts[i], fly->bouts->bouts[i], sizeof(BehaviorBout)*fly->bouts->num_bouts[i]);
  }

  return bouts;
}

const char *SVMFlyBehaviorSequence::get_base_feature_name(int ind) {
  return feature_defs[ind].name;
}

/*
 * Read the trx features from file and store per frame feature and timestamps 
 */
void SVMFlyBehaviorSequence::load_behavior_bout_features(void *b, BehaviorBoutFeatures *feature_cache) {
  FlyBehaviorBoutSequence *fly = (FlyBehaviorBoutSequence*)b;
  int i;

  FlyBehaviorFeatures d;
  double version;
  int len, fly_id, nflies, nfields, firstframe, lastframe;
  char movie[1001], mat[1001], folder[1001], fname[1001];
  
  if(fly->is_labeled) {
    strcpy(folder, fly->labelname);
    for(i = strlen(folder); i >= 0 && folder[i] != '/' && folder[i] != '\\'; i--) 
      folder[i] = '\0';
    sprintf(fname, "%s%s", folder, fly->trxname);
  } else
    sprintf(fname, "%s", fly->trxname);

  FILE *fin = fopen(fname, "rb");
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
    if(fly->is_labeled)
      assert(fread(&fly_id, sizeof(int), 1, fin) && fly_id == fly->fly_ids[i]);
    else
      assert(fread(&fly->fly_ids[i], sizeof(int), 1, fin));
  }
  assert(fread(&firstframe, sizeof(int), 1, fin) &&
	 fread(&lastframe, sizeof(int), 1, fin) &&
	 fread(&fly->sex, sizeof(char), 1, fin) && fread(&fly->fps, sizeof(double), 1, fin) &&
	 fread(&nfields, sizeof(int), 1, fin) && nfields == num_base_features+1 && nfields < MAX_BASE_FEATURES+1);

  if(fly->is_labeled) {
    assert(!strcmp(movie, fly->moviename) && !strcmp(mat, fly->matname) && nflies == fly->nflies && 
	   firstframe == fly->firstframe && lastframe == fly->lastframe);
  } else {
    fly->version = version;
    strcpy(fly->moviename, movie); strcpy(fly->matname, mat); fly->nflies = nflies;
    fly->firstframe = firstframe; fly->lastframe = lastframe;
  }

  if(feature_cache) {
    double *buff = (double*)malloc(sizeof(double)*(lastframe-firstframe+1));
    for(i = 0; i < nfields; i++) {
      assert(read_features(fin, &d, buff, lastframe-firstframe+1));
      if(i == 0) assert(!strcmp(d.name, "timestamp"));
      if(!strcmp(d.name, "timestamp")) {
	memcpy(feature_cache->frame_times, buff, sizeof(double)*(lastframe-firstframe+1));
      } else {
	assert(!strcmp(d.name, feature_defs[i-1].name));
	memcpy(feature_cache->features[i-1], buff, sizeof(double)*(lastframe-firstframe+1));
      }
    }
    free(buff);
  }

  fclose(fin);
}

/*
 * Read a behavior label from file
 */
void *SVMFlyBehaviorSequence::load_training_example(const char *fname, BehaviorGroups *behaviors) {
  FlyBehaviorBoutSequence *fly = (FlyBehaviorBoutSequence*)malloc(sizeof(FlyBehaviorBoutSequence));
  memset(fly, 0, sizeof(FlyBehaviorBoutSequence));
  fly->behaviors = behaviors;
  strcpy(fly->labelname, fname);
  
  unsigned int sz = sizeof(BehaviorBoutSequence)+(sizeof(int)+sizeof(BehaviorBout*)+2*sizeof(double))*behaviors->num;
  fly->bouts = (BehaviorBoutSequence*)my_malloc(sz);
  memset(fly->bouts, 0, sz);
  fly->bouts->bouts = (BehaviorBout**)(fly->bouts+1);
  fly->bouts->num_bouts = (int*)(fly->bouts->bouts+behaviors->num);
  fly->bouts->scores = (double*)(fly->bouts->num_bouts+behaviors->num);
  fly->bouts->losses = (double*)(fly->bouts->scores+behaviors->num);
  fly->bouts->behaviors = behaviors;
  fly->bouts->score = fly->bouts->loss = fly->bouts->slack = 0;

  if(strstr(fname, ".label")) {
    fly->is_labeled = true;

    char behname[1001];
    int len, nbehaviors, i, start_frame, end_frame, num_bouts, behavior;
    FILE *fin = fopen(fname, "rb");
    assert(fin);
    assert(fread(&fly->version, sizeof(double), 1, fin));
    assert(fread(&len, sizeof(int), 1, fin) && len < 1000 && len > 0 && 
	   (int)fread(fly->moviename, sizeof(char), len, fin)==len);
    fly->moviename[len] = '\0';
    assert(fread(&len, sizeof(int), 1, fin) && len < 1000 && len > 0 && 
	   (int)fread(fly->matname, sizeof(char), len, fin)==len);
    fly->matname[len] = '\0';
    assert(fread(&len, sizeof(int), 1, fin) && len < 1000 && len > 0 && 
	   (int)fread(fly->trxname, sizeof(char), len, fin)==len);
    fly->trxname[len] = '\0';
    assert(fread(&fly->nflies, sizeof(int), 1, fin) && fly->nflies > 0 && fly->nflies < 100);
    for(int i = 0; i < fly->nflies; i++)
      assert(fread(&fly->fly_ids[i], sizeof(int), 1, fin));

    assert(fread(&fly->firstframe, sizeof(int), 1, fin) && fly->firstframe >= 0 &&
	   fread(&fly->lastframe, sizeof(int), 1, fin) && fly->lastframe >= fly->firstframe);
    assert(fread(&nbehaviors, sizeof(int), 1, fin) && nbehaviors+2 == behaviors->behaviors[0].num_values);

    int ind = 0; // assume for now all behaviors are mutually exclusive
    
    for(i = 0; i < nbehaviors; i++) {
      assert(fread(&len, sizeof(int), 1, fin) && len < 1000 && len > 0 && 
	     (int)fread(behname, sizeof(char), len, fin)==len);
      behname[len] = '\0';
      assert(!strcmp(behname, behaviors->behaviors[ind].values[i+2].name));
    }
    
    assert(fread(&num_bouts, sizeof(int), 1, fin) && fly->bouts->num_bouts[ind] >= 0);
    fly->bouts->bouts[ind] = (BehaviorBout*)malloc(sizeof(BehaviorBout)*(num_bouts));
    fly->bouts->num_bouts[ind] = 0;
    int num_added = 0, lastframe = 0;
    for(i = 0; i < num_bouts; i++) {
      assert(fread(&start_frame, sizeof(int), 1, fin));
      assert(fread(&end_frame, sizeof(int), 1, fin) && end_frame >= start_frame);
      assert(fread(&behavior, sizeof(int), 1, fin)); 
      assert(behavior >= 0 && behavior < behaviors->behaviors[ind].num_values); 
      assert(start_frame >= lastframe);
      if(start_frame > lastframe && i > 0) {
	// Pad unlabelled bouts as the "Unknown" class
	fly->bouts->bouts[ind] = (BehaviorBout*)realloc(fly->bouts->bouts[ind], sizeof(BehaviorBout)*(num_bouts+num_added+1));
	fly->bouts->bouts[ind][fly->bouts->num_bouts[ind]].start_frame = lastframe;
	fly->bouts->bouts[ind][fly->bouts->num_bouts[ind]].end_frame = start_frame;
	fly->bouts->bouts[ind][fly->bouts->num_bouts[ind]++].behavior = 0;
	num_added++;
      }
      
      fly->bouts->bouts[ind][fly->bouts->num_bouts[ind]].start_frame = start_frame;
      fly->bouts->bouts[ind][fly->bouts->num_bouts[ind]].end_frame = lastframe = end_frame;
      fly->bouts->bouts[ind][fly->bouts->num_bouts[ind]].behavior = behavior;
      fly->bouts->num_bouts[ind]++;
    }
    for(i = 0; i < fly->bouts->num_bouts[ind]; i++) {
      fly->bouts->bouts[ind][i].start_frame -= fly->firstframe;
      fly->bouts->bouts[ind][i].end_frame -= fly->firstframe;
    }
    fclose(fin);
  } else {
    // Extract info from a .trx file instead of a .label file.  Used for testing (not training)
    fly->is_labeled = false;
    strcpy(fly->trxname, fname);
    load_behavior_bout_features(fly, NULL);
  }
  return fly;
}

/*
 * Called after fitting behaviors, to save its results.  Save label bouts
 */
void SVMFlyBehaviorSequence::save_example(void *b, void *d, const char *fname) {
  FlyBehaviorBoutSequence *fly = (FlyBehaviorBoutSequence*)b;
  BehaviorBoutSequence *bouts = (BehaviorBoutSequence*)d;
  load_from_bout_sequence(bouts, fly);
  
  int ind = 0; // assume for now all behaviors are mutually exclusive

  int len, i, nbehaviors=behaviors->behaviors[ind].num_values, s, e;
  FILE *fout = fopen(fname, "wb");
  assert(fout);
  assert(fwrite(&fly->version, sizeof(double), 1, fout));
  len = strlen(fly->moviename);
  assert(fwrite(&len, sizeof(int), 1, fout)  && 
	 (int)fwrite(fly->moviename, sizeof(char), len, fout)==len);
  len = strlen(fly->matname);
  assert(fwrite(&len, sizeof(int), 1, fout) &&
	 (int)fwrite(fly->matname, sizeof(char), len, fout)==len);
  len = strlen(fly->trxname);
  assert(fwrite(&len, sizeof(int), 1, fout) && len < 1000 && len > 0 && 
	 (int)fwrite(fly->trxname, sizeof(char), len, fout)==len);
  
  assert(fwrite(&fly->nflies, sizeof(int), 1, fout));
  for(int i = 0; i < fly->nflies; i++)
    assert(fwrite(&fly->fly_ids[i], sizeof(int), 1, fout));

  assert(fwrite(&fly->firstframe, sizeof(int), 1, fout) &&
	 fwrite(&fly->lastframe, sizeof(int), 1, fout));
  assert(fwrite(&nbehaviors, sizeof(int), 1, fout));


  for(i = 0; i < nbehaviors; i++) {
    len = (int)strlen(behaviors->behaviors[ind].values[i].name);
    assert(fwrite(&len, sizeof(int), 1, fout) &&
	   (int)fwrite(behaviors->behaviors[ind].values[i].name, sizeof(char), len, fout)==len);
  }
 
  assert(fwrite(&fly->bouts->num_bouts[ind], sizeof(int), 1, fout));
  for(i = 0; i < fly->bouts->num_bouts[ind]; i++) {
    s = fly->bouts->bouts[ind][i].start_frame+fly->firstframe;
    e = fly->bouts->bouts[ind][i].end_frame+fly->firstframe;
    assert(fwrite(&s, sizeof(int), 1, fout));
    assert(fwrite(&e, sizeof(int), 1, fout));
    assert(fwrite(&fly->bouts->bouts[ind][i].behavior, sizeof(int), 1, fout)); 
  }
  fclose(fout);
}

/*
 * Free the result of load_training_example (if necessary)
 */
void SVMFlyBehaviorSequence::free_data(void *b) {
  FlyBehaviorBoutSequence *fly = (FlyBehaviorBoutSequence*)b;
  free_behavior_bout_sequence(fly->bouts, behaviors->num);
  free(fly);
}
