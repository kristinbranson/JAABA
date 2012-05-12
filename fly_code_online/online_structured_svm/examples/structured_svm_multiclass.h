#ifndef __STRUCTURED_SVM_MULTICLASS_H
#define __STRUCTURED_SVM_MULTICLASS_H

#include "structured_svm.h"

#define VERSION       "V0.0"  /**< Software version that is saved to parameter definition files */

/**
 * @file structured_svm_multiclass.h
 * @brief Simple example of how to use this structured SVM package: implements routines for loss-sensitive multiclass SVM training and classification
 */


/**
 * @class MulticlassStructuredSVM
 * @brief Implements routines for loss-sensitive multiclass SVM training and classification
 */
class MulticlassStructuredSVM : public StructuredSVM {
 public:
  /******************* Functions overridden from StructuredSVM *************************/
  bool Load(const Json::Value &root);
  Json::Value Save();
  SparseVector Psi(StructuredData *x, StructuredLabel *y);
  double Inference(StructuredData *x, StructuredLabel *ybar, SparseVector *w, StructuredLabel *y_partial=NULL, StructuredLabel *y_gt=NULL, double w_scale=1);
  double Loss(StructuredLabel *y_gt, StructuredLabel *y_pred);
  double ImportanceSample(StructuredData *x, SparseVector *w, StructuredLabel *y_gt, struct _SVM_cached_sample_set *set, double w_scale=1);
  void MultiSampleUpdate(SVM_cached_sample_set *set, StructuredExample *ex, int R);

  /**
   * @brief returns a new MulticlassStructuredLabel label
   * @param x aMulticlassStructuredData  object associated with this label
   */
  StructuredLabel *NewStructuredLabel(StructuredData *x);

  
  /**
   * @brief returns a new MulticlassStructuredData label
   */
  StructuredData *NewStructuredData();

  StructuredDataset *LoadDataset(const char *fname);
  bool SaveDataset(StructuredDataset *d, const char *fname, int start_from);

  /******************* Customized functions and member variables *************************/
 public:
  MulticlassStructuredSVM();
  ~MulticlassStructuredSVM();

  /**
   * @brief Get the number of possible classes
   */
  int NumClasses() { return num_classes; }

  /**
   * @brief Get the number of features, where the dimensionality of psi is NumFeatures()*NumClasses()
   */
  int NumFeatures() { return num_features; }

 private:
  int num_classes;   /**< The number of possible classes */
  int num_features;  /**< The number of features in the original feature space.  The dimenisionality of the structured featurespace Psi(x,y) is num_classes*num_features */ 
  double **classConfusionCosts;  /**< A num_classesXnum_classes matrix storing the confusion cost for each pair of classes */

  double **alphas;
  int alphas_alloc_size;

  void ExportModel(const char *fname);
  void CreateSamples(struct _SVM_cached_sample_set *set, StructuredData *x, StructuredLabel *y_gt);
};

/**
 * @class MulticlassStructuredData
 * @brief Stores data x for a training example for a multiclass classifier.  This is just a
 * a regular feature vector
 */
class MulticlassStructuredData : public StructuredData {
  SparseVector *psi;  /**< A feature vector x of size num_features, where the structured feature space is of size num_features*num_classes */

  friend class MulticlassStructuredSVM;

public:
  MulticlassStructuredData() { 
    psi = NULL; 
  }

  ~MulticlassStructuredData() { 
    if(psi) delete psi; 
  }

  bool load(const Json::Value &x, StructuredSVM *s) {
    // Read a vector of features from a JSON encoding x
    if(psi) { delete psi; psi = NULL; } 
    psi = new SparseVector;
    if(!psi->load(x) || psi->Length() > ((MulticlassStructuredSVM*)s)->NumFeatures()) {
      fprintf(stderr, "Feature vector exceeds feature space dimensionality\n");
      return false;
    }
    return true;
  }

  Json::Value save(StructuredSVM *s) {
    // Encode a vector of features into a JSON object
    return psi ? psi->save() : Json::Value();
  }
};


/**
 * @class MulticlassStructuredLabel
 * @brief Stores a label y for a training example for a multiclass classifier.  This is just a
 * class id label
 */
class MulticlassStructuredLabel : public StructuredLabel {
  int class_id;  // The id of the class (this should vary from 1...num_classes)

  friend class MulticlassStructuredSVM;

public:
  /**
   * @brief Create a new MulticlassStructuredLabel
   * @param x the MulticlassStructuredData this label is associated with
   */
  MulticlassStructuredLabel(StructuredData *x) : StructuredLabel(x) {}

  bool load(const Json::Value &y, StructuredSVM *s) {
    if(!y.isMember("class")) { 
      fprintf(stderr, "Error reading MulticlassStructuredLabel, missing class\n"); 
      return false; 
    } else {
      class_id = y["class"].asInt();
      if(class_id < 1 || class_id > ((MulticlassStructuredSVM*)s)->NumClasses()) {
	fprintf(stderr, "Invalid class id %d\n", class_id);
	return false;
      }
      return true;
    }
  }

  Json::Value save(StructuredSVM *s) {
    // Encode the class id into a JSON object
    Json::Value y;
    y["class"] = class_id;
    return y;
  }
};

#endif


