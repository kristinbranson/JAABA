#include "structured_svm_multiclass.h"



/**
 * @file structured_svm_multiclass.cpp
 * @example structured_svm_multiclass.cpp
 *
 * This is an example of how to use the structured learning API to implement a custom structured learner.  This
 * example implements a multiclass SVM learner and classification with custom loss function
 *
 */

StructuredLabel *MulticlassStructuredSVM::NewStructuredLabel(StructuredData *x) { return new MulticlassStructuredLabel(x); }

StructuredData *MulticlassStructuredSVM::NewStructuredData() { return new MulticlassStructuredData; }


MulticlassStructuredSVM::MulticlassStructuredSVM() {
  num_classes = 0;
  num_features = 0;
  classConfusionCosts = NULL;
  eps = .01;
  C = 5000.0;
  window = 30000;
  featureScale = 1;
}

MulticlassStructuredSVM::~MulticlassStructuredSVM() {
  if(classConfusionCosts) 
    free(classConfusionCosts);
}

double MulticlassStructuredSVM::Inference(StructuredData *x, StructuredLabel *ybar, SparseVector *w, 
					  StructuredLabel *y_partial, StructuredLabel *y_gt) {
  int bestclass=-1, first=1;
  double score,bestscore=-1;

  MulticlassStructuredLabel *m_ybar = (MulticlassStructuredLabel*)ybar;
  MulticlassStructuredLabel *m_y_partial = y_partial ? (MulticlassStructuredLabel*)y_partial : NULL;

  // Loop through every possible class y and compute its score <w,Psi(x,y)>
  for(int class_id = 1; class_id <= num_classes; class_id++) {
    // By default, compute ybar = max_y <w,Psi(x,y)>, but it y_partial is non-null,
    // only consider labels that agree with y_partial 
    m_ybar->class_id=class_id;
    if(y_partial && m_y_partial->class_id != class_id)
      score = -INFINITY;
    else 
      score = w->dot(Psi(x, ybar));
   
    // If y_gt is non-null, compute ybar = max_y <w,Psi(x,y)>+Loss(y_gt,y) 
    if(y_gt)  
      score += Loss(y_gt, ybar);

    if(score > bestscore || first) {
      bestscore=score;
      bestclass=class_id;
      first=0;
    }
  }
  m_ybar->class_id = bestclass;

  return bestscore;
}

SparseVector MulticlassStructuredSVM::Psi(StructuredData *x, StructuredLabel *y) {
  // The dimensionality of Psi(x,y) is num_featuresXnum_classes, by concatenating
  // num_features features for each class. The entries for Psi are equal to x->psi for 
  // the true class y and 0 for all classes other than y.  
  MulticlassStructuredData *m_x = (MulticlassStructuredData*)x;
  MulticlassStructuredLabel *m_y = (MulticlassStructuredLabel*)y;
  return m_x->psi->shift((m_y->class_id-1)*num_features);
}

double MulticlassStructuredSVM::Loss(StructuredLabel *y_gt, StructuredLabel *y_pred) {
  // Computes the loss of prediction y_pred against the correct label y_gt. 
  MulticlassStructuredLabel *m_y_gt = (MulticlassStructuredLabel*)y_gt;
  MulticlassStructuredLabel *m_y_pred = (MulticlassStructuredLabel*)y_pred;
  if(classConfusionCosts)
    return classConfusionCosts[m_y_gt->class_id][m_y_pred->class_id];
  else 
    return m_y_gt->class_id == m_y_pred->class_id ? 0 : 1;
}



Json::Value MulticlassStructuredSVM::Save() {
  Json::Value root;
  
  root["version"] = VERSION;
  root["Num Classes"] = num_classes;
  root["Num Features"] = num_features;

  //if(classConfusionCosts) {
    Json::Value c;
    int n = 0;
    for(int i = 1; i <= num_classes; i++) {
      for(int j = 1; j <= num_classes; j++) {
	if((!classConfusionCosts && i != j) || (classConfusionCosts && classConfusionCosts[i][j])) {
	  Json::Value o;
	  o["c_gt"] = i;
	  o["c_pred"] = j;
	  o["loss"] = classConfusionCosts ? classConfusionCosts[i][j] : 1;
	  c[n++] = o;
	}
      }
    }
    root["Class Confusion Costs"] = c;
  //}

  return root;
}


bool MulticlassStructuredSVM::Load(const Json::Value &root) {
  fprintf(stdout, "loading parameters\n");
  if(strcmp(root.get("version", "").asString().c_str(), VERSION)) {
    fprintf(stderr, "Version of parameter file does not match version of the software"); 
    return false;
  }
  num_classes = root.get("Num Classes",0).asInt();
  num_features = root.get("Num Features",0).asInt();
  
  sizePsi = num_features*num_classes;

  if(root.isMember("Class Confusion Costs") && root["Class Confusion Costs"].isArray()) {
    classConfusionCosts = (double**)malloc((num_classes+1)*(sizeof(double*)+(num_classes+1)*sizeof(double)));
    double *ptr = (double*)(classConfusionCosts+(num_classes+1));
    for(int i = 0; i <= num_classes; i++, ptr += (num_classes+1)) {
      classConfusionCosts[i] = ptr;
      for(int j = 0; j <= num_classes; j++)
	classConfusionCosts[i][j] = 0;
    }
    Json::Value a = root["Class Confusion Costs"];
    for(int i = 0; i < (int)a.size(); i++) {
      int c_gt = a[i].get("c_gt",-1).asInt();
      int c_pred = a[i].get("c_pred",-1).asInt();
      double l = a[i].get("loss",0).asDouble();
      if(c_gt > num_classes || c_pred > num_classes || c_gt <= 0 || c_pred <= 0) {
	fprintf(stderr, "Error reading Class Confusion Costs\n");
	return false;
      }
      classConfusionCosts[c_gt][c_pred] = l;
    }
  }

  return true;
}



// Ordinarily, one need not override this function and should use the default StructuredSVM::LoadDataset() function
// instead.  We override because we want to import SVM^light format files
StructuredDataset *MulticlassStructuredSVM::LoadDataset(const char *fname) {
  if(debugLevel > 0) fprintf(stderr, "Reading dataset %s...", fname);
  
  Lock();

  bool detectNumFeatures = num_features == 0;
  bool detectNumClasses = num_classes == 0;

  FILE *fin = fopen(fname, "r");
  if(!fin) {
    fprintf(stderr, "Couldn't open dataset file %s\n", fname);
    Unlock();
    return NULL;
  }

  StructuredDataset *d = new StructuredDataset();
  char line[100000];
  
  while(fgets(line, 99999, fin) && strlen(line) > 1) {
    chomp(line);
    StructuredExample *ex = new StructuredExample;
    ex->x = NewStructuredData();
    ex->y = NewStructuredLabel(ex->x);
    SparseVector *psi = ((MulticlassStructuredData*)ex->x)->psi = new SparseVector();

    // Assume each example is a line containing a class id followed by %d:%f pairs 
    if(!sscanf(line, "%d", &((MulticlassStructuredLabel*)ex->y)->class_id) || !psi->from_string(strstr(line, " ")+1)) { 
      fprintf(stderr, "Error parsing dataset example %s\n", line);
      delete ex;
      delete d;
      fclose(fin);
      return false;
    }

    if(detectNumFeatures) num_features = my_max(num_features, psi->Length());
    if(detectNumClasses) num_classes = my_max(num_classes, ((MulticlassStructuredLabel*)ex->y)->class_id);

    d->AddExample(ex);
  }
  fclose(fin);
  Unlock();

  if(detectNumFeatures || detectNumClasses)
    sizePsi = num_features*num_classes;

  if(debugLevel > 0) fprintf(stderr, "done\n");

  return d;
}

// Ordinarily, one need not override this function and should use the default StructuredSVM::SaveDataset() function
// instead.  We override because we want to import SVM^light format files
bool MulticlassStructuredSVM::SaveDataset(StructuredDataset *d, const char *fname, int start_from) {
  if(debugLevel > 0 && start_from == 0) fprintf(stderr, "Saving dataset %s...", fname);

  Lock();

  FILE *fout = fopen(fname, start_from>0 ? "a" : "w");
  if(!fout) {
    fprintf(stderr, "Couldn't open dataset file %s for writing\n", fname);
    Unlock();
    return false;
  }

  char data[100000];
  for(int i = start_from; i < d->num_examples; i++) {
    ((MulticlassStructuredData*)d->examples[i]->x)->psi->to_string(data);
    fprintf(fout, "%d %s\n", ((MulticlassStructuredLabel*)d->examples[i]->y)->class_id, data);
  }
  fclose(fout);
  Unlock();

  if(debugLevel > 0 && start_from == 0) fprintf(stderr, "done\n");

  return true;
}


#ifndef NO_SERVER

#include "online_interactive_server.h"

/** 
 * @brief Run the structured learning server
 *
 * To run as a standalone structured learning algorithm with a fixed dataset file, run something like
 *   Train: ./online_interactive_server -p data/params.txt -d data/train.dat -o data/learned_model.txt
 *   Test:  ./online_interactive_server -i data/learned_model.txt -t data/test.dat data/predictions.dat
 *
 * where 
 *   data/params.txt is in the format of StructuredSVM::Load()
 *   data/train.dat is a training set in the format of  StructuredSVM::LoadDataset() (the default implementation 
 *       reads each training example as a line "<y> <x>", where <y> is in the format of StructuredLabel::read()
 *       and <x> is in the format of StructuredData::read()
 *   data/learned_model.txt is the file where the learned model is written
 *   data/test.dat is a testset in the format of  StructuredSVM::LoadDataset() 
 *   data/predictions.txt is the file where predictions for all labels are written, where each line
 *      corresponds to a test example in the format 
 *          "<y_predicted> <y_ground_truth> <loss> <score_prediction> <score_ground_truth>"
 *
 *
 *
 * To run as a server that trains in online fashion, allowing a client to interactively classify examples
 * and add new training examples, run something like
 *   ./online_interactive_server -P 8086 -p data/params.txt -d data/initial_train.txt
 *
 * where data/initial_train.dat is optional and 8086 is the port in which the serve listens on.
 *
 * See StructuredLearnerRpc for info on the network protocol
 *
 **/ 
int main(int argc, const char **argv) {
  StructuredLearnerRpc v(new MulticlassStructuredSVM);
  v.main(argc, argv);
}

#endif
