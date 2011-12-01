#include "blob.h"
#include "train.h"

#include <cv.h>
#include <highgui.h>
#include <ml.h>  

#ifdef WIN32
#include <windows.h>
#endif

//#include <wx/dir.h>

CvStatModel *allocate_classifier(BehaviorGroups *behaviors, int classifier_method) {
  return NULL;
}

void deallocate_classifier(CvStatModel *m, int classifier_method) {
  if(m) delete m;
}

char *get_classifier_name(char *classifier_name, const char *classifier_dir, const char *b_name, int method) {
  sprintf(classifier_name, "%s/%s.%s", classifier_dir, method == CLASSIFIER_SVM_STRUCT_BEHAVIOR ? "All" : b_name,
          g_classifer_extensions[method]);
  return classifier_name;
}

const char *g_classifer_extensions[400] = {
  "none", "dec.tree", "rand.forest", "boost", "tree.boost", "near.neighbor", "svm", "nbayes", "svmlight", "svmbehavior"
};

void strip_extension(char *fname) {
  int i;

  for(i = (int)strlen(fname)-1; i > 0; i--) {
    if(fname[i] == '.') {
      strcpy(fname+i, "");
      break;
    } else if(fname[i] == '/') {
      break;
    }
    fname[i] = '\0';
  }
}


void CreateDirectoryIfNecessary(const char *dirName) {
#ifndef WIN32
  char str[400]; 
  sprintf(str, "mkdir %s", dirName);
  system(str);
#else
  CreateDirectory(dirName, NULL);
#endif
}

// Read a blob outline file.  These are outputs of the MultiwormTracker, which
// contain a list of x,y coordinates defining the contour of a blob in world c oordinates
BlobSequence *import_blob_sequence(const char *fname, int num_spine_points) {
  BlobSequence *s = (BlobSequence*)malloc(sizeof(BlobSequence));
  FILE *fin;
  char line[100000];
  double spine_pts[1000];
  int n, i, j, ind, nn;
  char *ptr;
  Blob *b;
  char sname[400];
  double w, mm, tm;

  memset(s, 0, sizeof(BlobSequence));
  fprintf(stderr, "Importing %s...\n", fname);

  // By default, the imported file will be saved to the same location as the import file
  strcpy(s->from, fname);
  strcpy(s->fname, fname);
  strip_extension(s->fname);
  strcat(s->fname, ".anno");
  

  fin = fopen(fname, "r");
  if(!fin) {
    fprintf(stderr, "Blob file %s not found\n", fname);
    return NULL;
  }
  
  while(fgets(line, 99999, fin)) {
    s->frames = (Blob*)realloc(s->frames, (sizeof(Blob)*(s->num_frames+1)));
    b = &s->frames[s->num_frames];
    memset(b, 0, sizeof(Blob));

    // Read the frame time
    ptr = strtok(line, " ");
    if(!ptr || sscanf(ptr, "%lf", &b->frame_time) < 1) {
      fprintf(stderr, "Error reading line %d in %s: \n %s", b->num_pts,
	      fname, line);
      free_blob_sequence(s);
      fclose(fin);
      return NULL;
    }
    ptr = strtok(NULL, " "); 

    // Read a sequence of x,y coordinates
    n = 0;
    while(ptr) {
      b->contour = (double*)realloc(b->contour, 2*sizeof(double)*(b->num_pts+1));
      if(sscanf(ptr, "%lf", &b->contour[2*b->num_pts+n]) < 1) {
	fprintf(stderr, "Error reading line %d in %s: \n %s", b->num_pts,
		fname, line);
	free_blob_sequence(s);
	fclose(fin);
	return NULL;
      }
      if(n == 1) {
	b->num_pts++;
	n = 0;
      } else
	n = 1;
      ptr = strtok(NULL, " ");
    }

    s->num_frames++;
  }
  fclose(fin);

  // Read in the .spine file
  if(num_spine_points) {
    ind = 0;
    strcpy(sname, fname);
    strip_extension(sname);
    strcat(sname, ".spine");
    fin = fopen(sname, "r");
    if(fin) {
      while(fgets(line, 99999, fin)) {
	ptr = line;
	sscanf(ptr, "%lf", &tm);
	if(!(ptr=strstr(ptr, " ")))
	  break;
	
	while(ind < s->num_frames && s->frames[ind].frame_time < tm) 
	  ind++;
    
	if(s->frames[ind].frame_time == tm) {
	  s->frames[ind].num_model_pts = num_spine_points;
	  s->frames[ind].model_pts =(SpinePointLocation*)malloc(sizeof(SpinePointLocation)*s->frames[ind].num_model_pts);
	  memset(s->frames[ind].model_pts, 0, sizeof(sizeof(SpinePointLocation)*s->frames[ind].num_model_pts));

	  i = 0;
	  while(sscanf(ptr+1, "%lf %lf", &spine_pts[2*i], &spine_pts[2*i+1]) == 2) {
	    if(!(ptr=strstr(ptr+1, " ")))
	      break;
	    if(!(ptr=strstr(ptr+1, " ")))
	      break;
	    i++;
	  }
      
	  for(j = 0; j < num_spine_points; j++) {
	    mm = .5+j*(i-1)/(double)(num_spine_points);
	    nn = (int)mm;
	    w = mm-nn;
	    s->frames[ind].model_pts[num_spine_points-j-1].x = spine_pts[2*nn]*(1-w) + (w ? w*spine_pts[2*(nn+1)] : 0);
	    s->frames[ind].model_pts[num_spine_points-j-1].y = spine_pts[2*nn+1]*(1-w) + (w ? w*spine_pts[2*(nn+1)+1] : 0);
	  }
	}
      }
      fclose(fin);
    }
  }

  return s;
}


// Free the data associated with a call to read_blob_file
void free_blob_sequence(BlobSequence *s) {
  int i;

  for(i = 0; i < s->num_frames; i++) 
      free_blob(&s->frames[i]);

  free(s->frames);
  free(s);
}



MultiBlobSequence *import_multi_blob_sequence(const char *dir_name, int num_spine_points) {
  int l;
  FILE *pin;
  MultiBlobSequence *m = (MultiBlobSequence*)malloc(sizeof(MultiBlobSequence));
  char fname[400], sys[400], line[1000];

  m->num_blobs = 0;
  m->blobs = NULL;
  strcpy(m->fname, dir_name);

#ifdef WIN32
  WIN32_FIND_DATA FileInformation; 
  char strPattern[1000];
  sprintf(strPattern, "%s/*.outline", dir_name);
  HANDLE hFile = ::FindFirstFile(strPattern, &FileInformation);
  if(hFile != INVALID_HANDLE_VALUE) {
	  do {
		  m->blobs = (BlobSequence**)realloc(m->blobs, (m->num_blobs+1)*sizeof(BlobSequence*));
	      sprintf(fname, "%s/%s", dir_name, FileInformation.cFileName);
		  m->blobs[m->num_blobs++] = import_blob_sequence(fname, num_spine_points);
	  } while(::FindNextFile(hFile, &FileInformation) == TRUE);
  }

  /*
  wxArrayString files;
  size_t nfiles = wxDir::GetAllFiles(wxString(dir_name, wxConvUTF8), &files, wxT("*.outline"), wxDIR_FILES|wxDIR_DIRS);
  for (size_t i = 0; i < nfiles; i++) {
    m->blobs = (BlobSequence**)realloc(m->blobs, (m->num_blobs+1)*sizeof(BlobSequence*));
    sprintf(fname, "%s", files.Item(i).mb_str());
    m->blobs[m->num_blobs++] = import_blob_sequence(fname, num_spine_points);
  }*/
#else

  sprintf(sys, "ls %s/*.outline", dir_name);
  pin = popen(sys, "r");
  while(fgets(line, 999, pin)) {
    l = strlen(line);
    if(l && line[l-1] == '\n')
      line[l-1] = '\0';

    m->blobs = (BlobSequence**)realloc(m->blobs, (m->num_blobs+1)*sizeof(BlobSequence*));
    sprintf(fname, "%s", line);
    m->blobs[m->num_blobs++] = import_blob_sequence(fname, num_spine_points);
  }
  pclose(pin);
#endif

  return m;
}

void save_multi_blob_sequence(MultiBlobSequence *m, BehaviorGroups *behaviors, int num_model_pts, int num_orientations)  {
  int i;
  
  for(i = 0; i < m->num_blobs; i++) 
    save_blob_sequence(m->blobs[i]->fname, m->blobs[i], behaviors, num_model_pts, num_orientations);
}

void free_multi_blob_sequence(MultiBlobSequence *m) {
  int i;
  
  for(i = 0; i < m->num_blobs; i++) 
    free_blob_sequence(m->blobs[i]);
  
  free(m);
}


char *save_contour(double *c, int num, char *str) {
  int i;

  for(i = 0; i < num; i++) {
    sprintf(str, "%f\t%f\t", (float)c[2*i], (float)c[2*i+1]);
    str += strlen(str);
  }
  return str;
}

char *save_spine_point(SpinePointLocation *pt, char *str) {
  sprintf(str, "%f\t%f\t%d\t", (float)pt->x, (float)pt->y, pt->orientation);
  return str+strlen(str);
}

char *save_blob(Blob *b, BehaviorGroups *behaviors, int num_model_pts, char *str) {
  int i;
  char beh_str[1000], tmp[400];

  strcpy(beh_str, "");
  for(i = 0; i < behaviors->num; i++) {
    sprintf(tmp, "%s%s:%s", i ? ";" : "", behaviors->behaviors[i].name, 
	    behaviors->behaviors[i].values[b->behaviors[i]].abbreviation);
    strcat(beh_str, tmp);
  }

  sprintf(str, "%f\t%d\t%s\t%d\t#\t%d\t", (float)b->frame_time, b->is_manual, beh_str, b->is_good, b->num_model_pts);
  str = str + strlen(str);
  for(i = 0; i < num_model_pts; i++) {
    if(i < b->num_model_pts)
      str = save_spine_point(&b->model_pts[i], str);
    else {
      strcpy(str, "\t\t\t"); 
      str += strlen(str);
    }
  }

  sprintf(str, "##\t%d\t", b->num_pts);
  str += strlen(str);
  return save_contour(b->contour, b->num_pts, str);
}

void save_blob_sequence(const char *fname, BlobSequence *b, BehaviorGroups *behaviors, 
			int num_model_pts, int num_orientations) {
  char line[100000];
  FILE *fout = fopen(fname, "w");
  int i;
  assert(fout);

  if(b->num_frames) {
    fprintf(fout, "%%\tNum Orientations\tImported From\tFile\n%%\t%d\t%s\t%s\n%%\n", 
	    num_orientations, b->from, b->fname);
    fprintf(fout, "Frame Num\tFrame Time\tIs Manual\tBehavior\tIs Valid\t#\tNum Model Pts\t");
    for(i = 0; i < num_model_pts; i++) {
      fprintf(fout, "x%d\ty%d\ttheta%d\t",i+1,i+1,i+1);
    }
    fprintf(fout, "##\tNum Contour Points\tx1\ty1\t...");
  }

  for(i = 0; i < b->num_frames; i++) {
    save_blob(&b->frames[i], behaviors, num_model_pts, line);
    fprintf(fout, "\n%d\t%s", i+1, line);
  }
  fclose(fout);
}

char *load_contour(char *str, double **contour, int *num) {
  int i;
  char *ptr;

  if(sscanf(str, "%d", num) < 1) return NULL;
  *contour = (double*)malloc(2*sizeof(double)*(*num));
  if(!(str = strstr(str, "\t"))) return NULL;
  for(i = 0; i < (*num); i++) {
    if(sscanf(str+1, "%lf", &(*contour)[2*i]) < 1) return NULL;
    if(!(str = strstr(str+1, "\t"))) return NULL;
    if(sscanf(str+1, "%lf", &(*contour)[2*i+1]) < 1) return NULL;
    if(i < (*num-1)) {
      if(!(str = strstr(str+1, "\t"))) return NULL;
    } else {
      if(!(ptr = strstr(str+1, "\t")) && !(ptr = strstr(str+1, "\n"))) 
	return str + strlen(str);
      else str = ptr;
    }
  }
  return str+1;
}

char *load_spine_point(char *str, SpinePointLocation *pt) {
  if(sscanf(str, "%lf", &pt->x) < 1) return NULL;
  if(!(str = strstr(str, "\t"))) return NULL;
  if(sscanf(str+1, "%lf", &pt->y) < 1) return NULL;
  if(!(str = strstr(str+1, "\t"))) return NULL;
  if(sscanf(str+1, "%d", &pt->orientation) < 1) return NULL;
  if(!(str = strstr(str+1, "\t"))) return NULL;
  return str+1;
}
  

char *load_blob(char *str, BehaviorGroups *behaviors,  Blob *b) {
  int i, frame;
  char beh_str[400], *ptr;

  memset(b, 0, sizeof(Blob));

  if(sscanf(str, "%d", &frame) < 1) return NULL;
  if(!(str = strstr(str, "\t"))) return NULL;
  if(sscanf(str+1, "%lf", &b->frame_time) < 1) return NULL;
  if(!(str = strstr(str+1, "\t"))) return NULL;
  if(sscanf(str+1, "%d", &b->is_manual) < 1) return NULL;
  if(!(ptr = strstr(str+1, "\t"))) return NULL;
  if(!(str = strstr(ptr+1, "\t"))) return NULL;
  strncpy(beh_str, ptr+1, str-ptr-1);
  beh_str[str-ptr-1] = '\0';

  if(!find_behavior(behaviors, b->behaviors, beh_str)) 
    return NULL;

  if(sscanf(str+1, "%d", &b->is_good) == 1) {
    if(!(str = strstr(str+1, "\t"))) return NULL;
  } else
    b->is_good = 1;

  if(str[1] != '#') return NULL;
  if(!(str = strstr(str+2, "\t"))) return NULL;
  if(sscanf(str+1, "%d", &b->num_model_pts) < 1) return NULL;
  if(!(str = strstr(str+1, "\t"))) return NULL;
  str++;
  
  if(b->num_model_pts <= 0)
    b->is_good = 0;
  b->model_pts = (SpinePointLocation*)malloc(sizeof(SpinePointLocation)*b->num_model_pts);
  for(i = 0; i < b->num_model_pts; i++) {
    if(!(str = load_spine_point(str, &b->model_pts[i])))
      return NULL;
  }
  while(str && str[0] != '#' && str[1] != '#') {
    if(!(str = strstr(str, "\t"))) return NULL;
    str++;
  }

  if(str[0] != '#' && str[1] != '#') return NULL;
  if(!(str = strstr(str+2, "\t"))) return NULL;
  if(sscanf(str+1, "%d", &b->num_pts) < 1) return NULL;
  if(!(str = strstr(str, "\t"))) return NULL;
  return load_contour(str+1, &b->contour, &b->num_pts);
}

BlobSequence *load_blob_sequence(const char *fname, BehaviorGroups *behaviors) {
  int num_orientations, i, j;
  char line[100000];
  FILE *fin = fopen(fname, "r");
  BlobSequence *b = (BlobSequence*)malloc(sizeof(BlobSequence));
  memset(b, 0, sizeof(BlobSequence));
  strcpy(b->fname, fname);
  assert(fin);

  while(fgets(line, 99999, fin)) {
    if(line[0] == '%') {
      sscanf(line, "%d\t%s\t%s", &num_orientations, b->from, b->fname);
      continue;
    } else if(line[0] == 'F')
      continue;
    b->frames = (Blob*)realloc(b->frames, sizeof(Blob)*(b->num_frames+1));
    assert(load_blob(line, behaviors, &b->frames[b->num_frames]));

    b->num_frames++;
  }
  fclose(fin);

  // HACK to make it ignore excel imported labels by default
  for(i = 0; i < b->num_frames; i++) {
    if(b->frames[i].is_manual == 2) {
      j = i;
      while(j < i && b->frames[j].is_manual == 2)
	j++;
      if(j-i <= 2) {
	//b->frames[i].is_manual = b->frames[j-1].is_manual = 0; // CSC 20110210:
	b->frames[i].is_manual = b->frames[j].is_manual = 0; 
      }
    }
  }

  return b;
}

void free_blob(Blob *b) {
  if(b) {
    if(b->contour) free(b->contour);
    if(b->model_pts) free(b->model_pts);
    if(b->features) free(b->features);
    if(b->fixed_contour) free(b->fixed_contour);

    assert(!b->data);  // the caller should free this and set to NULL before freeing a Blob
  }
}


// Get the bounding box around all points in b->contour
void blob_get_bounding_box(Blob **blobs, int num_blobs, double *min_x, double *min_y, double *max_x, double *max_y) {
  int i, j;
  Blob *b;

  *min_x = INFINITY;
  *max_x = -INFINITY;
  *min_y = INFINITY;
  *max_y = -INFINITY;
  
  // Find the bounding box around the set of contour points
  for(j = 0; j < num_blobs; j++) {
    b = blobs[j];
    for(i = 0; i < b->num_pts; i++) {
      if(b->contour[2*i] < *min_x) *min_x = b->contour[2*i];
      if(b->contour[2*i] > *max_x) *max_x = b->contour[2*i];
      if(b->contour[2*i+1] < *min_y) *min_y = b->contour[2*i+1];
      if(b->contour[2*i+1] > *max_y) *max_y = b->contour[2*i+1];
    }
  }
}


void free_behaviors(BehaviorGroups *behs) {
  int i, j;

  if(behs) {
    if(behs->behaviors) {
      for(i = 0; i < behs->num; i++) {
	for(j = 0; j < behs->behaviors[i].num_values; j++)
	  if(behs->behaviors[i].values[j].classifier && behs->behaviors[i].values[j].classifier != behs->behaviors[i].classifier)
	    deallocate_classifier(behs->behaviors[i].values[j].classifier, behs->behaviors[i].values[j].classifier_method);
	free(behs->behaviors[i].values);
	if(behs->behaviors[i].classifier)
	  deallocate_classifier(behs->behaviors[i].classifier, behs->behaviors[i].classifier_method);
      }
      free(behs->behaviors);
    }
    if(behs->classifier)
      deallocate_classifier(behs->classifier, behs->classifier_method);
    free(behs);
  }
}

  
BehaviorGroups *load_behaviors(const char *fname, const char *classifier_dir) {
  BehaviorGroups *behs = (BehaviorGroups*)malloc(sizeof(BehaviorGroups)); 
  BehaviorGroup *b;
  char line[10000], classifier_name[400];
  char *ptr, *ptr2;
  FILE *fin = fopen(fname, "r");
  assert(fin);
  memset(behs, 0, sizeof(BehaviorGroups));

  behs->num = 0;
  behs->behaviors = NULL;

  fgets(line, 999, fin);
  if(sscanf(line, "%d %d", &behs->classifier_method, &behs->is_multiclass) != 2) {
    behs->classifier_method = 0;
    behs->is_multiclass = 0;
  } 
  while(fgets(line, 999, fin)) {
    chomp(line);
    if(line[0] == '*') {
      // Add a new behavior group
      behs->behaviors = (BehaviorGroup*)realloc(behs->behaviors, sizeof(BehaviorGroup)*(behs->num+1));
      b = behs->behaviors + behs->num++;
      memset(b, 0, sizeof(BehaviorGroup));
      assert(sscanf(line+1, "%s %d %d", b->name, &b->classifier_method, &b->is_multiclass) == 3);
      if(b->classifier_method) {
	get_classifier_name(classifier_name, classifier_dir, b->name, b->classifier_method);
	if(FileExists(classifier_name)) {
	  b->classifier = allocate_classifier(behs, b->classifier_method);
	  b->classifier->load(classifier_name);
	} else 
	  fprintf(stderr, "ERROR: classifier %s doesn't exist\n", classifier_name);
      }

#if ADD_DUMMY_BEHAVIORS
      b->num_values = 2;
      b->values = (BehaviorValue*)malloc(2*sizeof(BehaviorValue)); memset(b->values, 0, 2*sizeof(BehaviorValue));
      strcpy(b->values[0].name, "*Unknown*");
      strcpy(b->values[0].abbreviation, "");
      b->values[0].color = 0x010101;
      strcpy(b->values[1].name, "*Tracker Failure*");
      strcpy(b->values[1].abbreviation, "BAD");
      b->values[1].color = 0x6e6e6e;
#endif
    } else if(!strlen(line))
      continue;
    else {
      b->values = (BehaviorValue*)realloc(b->values, (b->num_values+1)*sizeof(BehaviorValue)); 
      memset(&b->values[b->num_values], 0, sizeof(BehaviorValue));
      assert((ptr = strtok(line, "\t")) != NULL);
      assert((ptr2 = strtok(NULL, "\t")) != NULL);
      strcpy(b->values[b->num_values].name, ptr);
      strcpy(b->values[b->num_values].abbreviation, ptr2);   
      chomp(b->values[b->num_values].abbreviation);

      assert((ptr2 = strtok(NULL, "\t")) != NULL && sscanf(ptr2, "%d", &b->values[b->num_values].classifier_method));
      if(b->values[b->num_values].classifier_method) {
	sprintf(classifier_name, "%s/%s_%s.%s", classifier_dir, b->name, b->values[b->num_values].name, 
		g_classifer_extensions[b->values[b->num_values].classifier_method]);
	if(FileExists(classifier_name)) {
	  b->values[b->num_values].classifier = allocate_classifier(behs, b->values[b->num_values].classifier_method);
	  b->values[b->num_values].classifier->load(classifier_name);
	} else 
	  fprintf(stderr, "ERROR: classifier %s doesn't exist\n", classifier_name);
      }
      ptr2 = strtok(NULL, "\t");
      if(ptr2)
	sscanf(ptr2, "%x", &b->values[b->num_values].color);
  
      b->num_values++;
    }
  }

  if(behs->classifier_method > 0) {
    get_classifier_name(classifier_name, classifier_dir, "All", 
			behs->classifier_method);
    if(FileExists(classifier_name)) {
      behs->classifier = allocate_classifier(behs, behs->classifier_method);
      behs->classifier->load(classifier_name);
    } else 
      fprintf(stderr, "ERROR: classifier %s doesn't exist\n", classifier_name);
  }

  return behs;
}

int save_behaviors(const char *fname, BehaviorGroups *groups) {
  BehaviorGroup *b;
  int i, j;
  FILE *fout = fopen(fname, "w");
  if(!fout)
    return 0;
  
  fprintf(fout, "%d %d\n", groups->classifier_method, groups->is_multiclass);
  for(i = 0; i < groups->num; i++) {
    b = &groups->behaviors[i];
    fprintf(fout, "*%s %d %d\n", b->name, b->classifier_method, b->is_multiclass);
    for(j = 2; j < b->num_values; j++) {
      fprintf(fout, "%s\t%s\t%d\t%x\n", b->values[j].name, b->values[j].abbreviation,
	      b->values[j].classifier_method, b->values[j].color);
    }
    fprintf(fout, "\n");
  }
  fclose(fout);

  return 1;
}


int find_behavior(BehaviorGroups *behaviors, int *ids, const char *str) {
  int i, j;
  char tmp[1000], *ptr, *ptr2, abbr[400], name[400];
  char *tmpP = tmp;
  strcpy(tmp, str);

  memset(ids, 0, sizeof(int)*behaviors->num);

  while((ptr=strtok(tmpP, ";")) != NULL) {
    tmpP = NULL;
    if((ptr2 = strstr(ptr, ":")) == NULL) {
      return 0;

      /*int id;
      char *abbrevs[36] = { "", "BAD", "SH", "WH", "LHW", "RHW", "LHS", "RHS",
				 "LHES", "RHES", "LUHW", "RUHW", "LUHS", "RUHS", "LUHES", 
				 "RUHES", "RTT", "RHWTT", "LTT", "LHWTT", "SRL", "SLR",
				 "MF", "MB", "MFH", "ML", "MFRHES", "MFRUHS", "MFRHS", "MFRUHW",
				 "MFLUHW", "MBLHS", "MBTWUL", "MFRTT", "R", "D" };
      if(sscanf(ptr, "%d", &id) == 1) 
	ptr = abbrevs[id];
      if(!strcmp(ptr, "")) { ids[0]=ids[1]=ids[2]=0; return 1; }
      else if(!strcmp(ptr, "BAD")) { ids[0]=ids[1]=ids[2]=1; return 1; }
      else if(!strcmp(ptr, "SH")) { ids[0]=12;   ids[1]=4;   ids[2]=2; return 1; }
      else if(!strcmp(ptr, "WH")) { ids[0]=12;   ids[1]=3;   ids[2]=2; return 1; }
      else if(!strcmp(ptr, "LHW")) { ids[0]=2;   ids[1]=3;   ids[2]=2; return 1; }
      else if(!strcmp(ptr, "RHW")) { ids[0]=13;   ids[1]=3;   ids[2]=2; return 1; }
      else if(!strcmp(ptr, "LHS")) { ids[0]=3;   ids[1]=3;   ids[2]=2; return 1; }
      else if(!strcmp(ptr, "RHS")) { ids[0]=14;   ids[1]=3;   ids[2]=2; return 1; }
      else if(!strcmp(ptr, "LHES")) { ids[0]=4;   ids[1]=3;   ids[2]=2; return 1; }
      else if(!strcmp(ptr, "RHES")) { ids[0]=15;   ids[1]=3;   ids[2]=2; return 1; }
      else if(!strcmp(ptr, "LUHW")) { ids[0]=2;   ids[1]=2;   ids[2]=2; return 1; }
      else if(!strcmp(ptr, "RUHW")) { ids[0]=13;   ids[1]=2;   ids[2]=2; return 1; }
      else if(!strcmp(ptr, "LUHS")) { ids[0]=3;   ids[1]=2;   ids[2]=2; return 1; }
      else if(!strcmp(ptr, "RUHS")) { ids[0]=14;   ids[1]=2;   ids[2]=2; return 1; }
      else if(!strcmp(ptr, "LUHES")) { ids[0]=4;   ids[1]=2;   ids[2]=2; return 1; }
      else if(!strcmp(ptr, "RUHES")) { ids[0]=15;   ids[1]=2;   ids[2]=2; return 1; }
      else if(!strcmp(ptr, "RTT")) { ids[0]=19;   ids[1]=2;   ids[2]=2; return 1; }
      else if(!strcmp(ptr, "RHWTT")) { ids[0]=19;   ids[1]=3;   ids[2]=2; return 1; }
      else if(!strcmp(ptr, "LTT")) { ids[0]=8;   ids[1]=2;   ids[2]=2; return 1; }
      else if(!strcmp(ptr, "LHWTT")) { ids[0]=8;   ids[1]=3;   ids[2]=2; return 1; }
      else if(!strcmp(ptr, "SRL")) { ids[0]=22;   ids[1]=2;   ids[2]=2; return 1; }
      else if(!strcmp(ptr, "SLR")) { ids[0]=11;   ids[1]=2;   ids[2]=2; return 1; }
      else if(!strcmp(ptr, "MF")) { ids[0]=12;   ids[1]=2;   ids[2]=3; return 1; }
      else if(!strcmp(ptr, "MB")) { ids[0]=12;   ids[1]=2;   ids[2]=4; return 1; }
      else if(!strcmp(ptr, "MFH")) { ids[0]=12;   ids[1]=3;   ids[2]=3; return 1; }
      else if(!strcmp(ptr, "MBH")) { ids[0]=12;   ids[1]=3;   ids[2]=4; return 1; }
      else if(!strcmp(ptr, "ML")) { ids[0]=2;   ids[1]=2;   ids[2]=3; return 1; }
      else if(!strcmp(ptr, "MFRHES")) { ids[0]=15;   ids[1]=3;   ids[2]=3; return 1; }
      else if(!strcmp(ptr, "MFRUHS")) { ids[0]=14;   ids[1]=2;   ids[2]=3; return 1; }
      else if(!strcmp(ptr, "MFRHS")) { ids[0]=14;   ids[1]=3;   ids[2]=3; return 1; }
      else if(!strcmp(ptr, "MFRUHW")) { ids[0]=13;   ids[1]=2;   ids[2]=3; return 1; }
      else if(!strcmp(ptr, "MFLUHW")) { ids[0]=2;   ids[1]=2;   ids[2]=3; return 1; }
      else if(!strcmp(ptr, "MBLHS")) { ids[0]=3;   ids[1]=3;   ids[2]=4; return 1; }
      else if(!strcmp(ptr, "MBTWUL")) { ids[0]=2;   ids[1]=2;   ids[2]=4; return 1; }
      else if(!strcmp(ptr, "MFRTT")) { ids[0]=19;   ids[1]=2;   ids[2]=3; return 1; }
      else if(!strcmp(ptr, "R")) { ids[0]=0;   ids[1]=0;   ids[2]=5; return 1; }
      else if(!strcmp(ptr, "D")) { ids[0]=0;   ids[1]=0;   ids[2]=6; return 1; }
      else { assert(0); }*/
    }

    strncpy(name, ptr, ptr2-ptr); name[ptr2-ptr] = '\0';
    strcpy(abbr, ptr2+1);
    for(i = 0; i < behaviors->num; i++) {
      if(!strcmp(behaviors->behaviors[i].name, name)) {
	for(j = 0; j < behaviors->behaviors[i].num_values; j++) {
	  if(!strcmp(behaviors->behaviors[i].values[j].abbreviation, abbr))
	    ids[i] = j;
	}
      }
    }
  }
  return 1;
}


