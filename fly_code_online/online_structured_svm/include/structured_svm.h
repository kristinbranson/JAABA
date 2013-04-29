#ifndef __STRUCTURED_SVM_H
#define __STRUCTURED_SVM_H

#include "sparse_vector.h"

#include <omp.h>
#include <time.h>

/**
 * @file structured_svm.h
 * @brief Implements an abstract class that should be extended to define one's own customstructured learner
 */

/**
 * @enum StructuredPredictionOptimizationMethod
 * @brief Optimization method used for structured learning
 */
typedef enum {
  SPO_CUTTING_PLANE,   /**< SVM^struct */
  SPO_CUTTING_PLANE_SMALL_BATCH_SIZE,   /**< SVM^struct with batch size of 1% of the dataset */

  // Online Algorithms, which process one training example i per iteration:
  SPO_SGD,             /**< stochastic gradient descent: w^t = w^{t-1} - step_size*(w^{t-1} + grad_i(w^{t-1})),   step_size=1/(lambda*t) */
  SPO_SGD_PEGASOS,     /**< Do SGD step, then ensure that ||w^t||^2 <= 1/lambda   (downscale w^t if this doesn't hold) */
  SPO_DUAL_UPDATE,     /**< w^t = w^{t-1}(t-1)/t - alpha*grad_i(w^{t-1}),   where alpha is chosen to maximize dual objective */
  SPO_DUAL_UPDATE_WITH_CACHE,     /**< Same as SPO_DUAL_UPDATE, but also run dual update steps in a background process on labels from earlier iterations */
  SPO_DUAL_MULTI_SAMPLE_UPDATE,   /**< Sample multiple labels per iteration (instead of just getting the "most violated constraint"), then optimize parameters jointly */
  SPO_DUAL_MULTI_SAMPLE_UPDATE_WITH_CACHE,  /**< Same as SPO_DUAL_MULTI_SAMPLE_UPDATE, but also run multi-sample dual update steps in a background process  */

  SPO_MAP_TO_BINARY,    /**< Map into a binary classification problem, by randomly selecting samples for negative examples */
  SPO_MAP_TO_BINARY_MINE_HARD_NEGATIVES,    /**< Map into a series of binary classification problems, at each round randomly selecting hard negative examples using the current classifier */
  SPO_FIXED_SAMPLE_SET  /**< Train using a fixed pre-determined sample set, without ever calling Inference() or ImportanceSample() */
} StructuredPredictionOptimizationMethod;


class StructuredExample;
class StructuredDataset;
class StructuredSVM;
struct _SVM_cached_sample;
struct _SVM_cached_sample_set;

///@cond
#define MAX_DRIFT_BITS 10.0

inline double get_runtime() { 
  return ((double)clock()/(double)CLOCKS_PER_SEC); 
}

inline void OptimizationMethodToString(StructuredPredictionOptimizationMethod method, char *str) {
  switch(method) {
  case SPO_CUTTING_PLANE: strcpy(str, "cutting_plane"); break;
  case SPO_CUTTING_PLANE_SMALL_BATCH_SIZE: strcpy(str, "cutting_plane_small_batch_size"); break;
  case SPO_SGD: strcpy(str, "SGD"); break;
  case SPO_SGD_PEGASOS: strcpy(str, "SGD_pegasos"); break;
  case SPO_DUAL_UPDATE: strcpy(str, "online_dual_ascent"); break;
  case SPO_DUAL_UPDATE_WITH_CACHE: strcpy(str, "online_dual_ascent_with_cache"); break;
  case SPO_DUAL_MULTI_SAMPLE_UPDATE: strcpy(str, "multi_sample"); break;
  case SPO_DUAL_MULTI_SAMPLE_UPDATE_WITH_CACHE: strcpy(str, "multi_sample_with_cache"); break;
  case SPO_MAP_TO_BINARY: strcpy(str, "binary_classification"); break;
  case SPO_MAP_TO_BINARY_MINE_HARD_NEGATIVES: strcpy(str, "mine_hard_negatives"); break;
  case SPO_FIXED_SAMPLE_SET: strcpy(str, "fixed_sample_set"); break;
  default:  strcpy(str, "unknown"); break;
  }
}
inline StructuredPredictionOptimizationMethod OptimizationMethodFromString(const char *str) {
  if(!strcmp(str, "cutting_plane")) return SPO_CUTTING_PLANE; 
  else if(!strcmp(str, "cutting_plane_small_batch_size")) return SPO_CUTTING_PLANE_SMALL_BATCH_SIZE; 
  else if(!strcmp(str, "SGD")) return SPO_SGD; 
  else if(!strcmp(str, "SGD_pegasos")) return SPO_SGD_PEGASOS; 
  else if(!strcmp(str, "online_dual_ascent")) return SPO_DUAL_UPDATE; 
  else if(!strcmp(str, "online_dual_ascent_with_cache")) return SPO_DUAL_UPDATE_WITH_CACHE; 
  else if(!strcmp(str, "multi_sample")) return SPO_DUAL_MULTI_SAMPLE_UPDATE; 
  else if(!strcmp(str, "multi_sample_with_cache")) return SPO_DUAL_MULTI_SAMPLE_UPDATE_WITH_CACHE; 
  else if(!strcmp(str, "binary_classification")) return SPO_MAP_TO_BINARY; 
  else if(!strcmp(str, "mine_hard_negatives")) return SPO_MAP_TO_BINARY_MINE_HARD_NEGATIVES; 
  else if(!strcmp(str, "fixed_sample_set")) return SPO_FIXED_SAMPLE_SET; 
  else return (StructuredPredictionOptimizationMethod)-1;
}

///@endcond

/** @mainpage Online Structured SVM Optimizer
 *
 * This software package implements an optimizer for training max-margin structured prediction algorithms (structured SVMs)
 * using online subgradient methods.  The goal is to train a function 
 * \f[ g : X \to Y \f]
 * using training examples \f$ (X_1,Y_1), ..., (X_n,Y_n) \f$, where \f$ Y \f$ can be some multi-dimensional structured output.
 * It has a network interface such that annotators can add new training examples in online
 * fashion as it trains.  It also supports active labeling (using the current learned model to accelerate labeling of a new
 * example) via the same network interface.  The worst case time complexity is 
 * \f$ O(\frac{R^2 T}{\lambda \epsilon} + T n) \f$, which is an improvement over the worst case time complexity 
 * for SVM^struct \f$ O(n \frac{R^2 T}{\lambda \epsilon} ) \f$.  Here \f$ \epsilon \f$ is the approximation level from
 * the minimal achievable training error, \f$ \lambda \f$ is a regularization constant, \f$ T \f$ is the amount of time
 * needed to solve an inference problem, and \f$ R \f$ is a bound on the magnitude of the feature space.  This optimization 
 * package currently only supports linear kernels.  
 * 
 * At test time, the predicted label is 
 * \f[ g(X;\mathbf{w}) = \arg\max_{Y} \langle \mathbf{w}, \Psi(X,Y) \rangle \f]
 * where \f$ \mathbf{w} \f$ is a vector of model parameters, and \f$ \Psi(X,Y) \f$ is a vector of features extracted
 * from input \f$ X \f$ with respect to a candidate label \f$ Y \f$.  For example, for sliding window object detection
 * \f$ X \f$ is an image and  \f$ \Psi(X,Y) \f$ extracts features from a bounding box defined by \f$ Y \f$.  Training
 * solves for the value of \f$ \mathbf{w} \f$ that minimizes the training error
 * \f[ f(\mathbf{w}) = \frac{\lambda}{2} \|\mathbf{w}\|^2 + \frac{1}{n} \sum_{i=1}^n \max_{\bar{Y}_i} \left( \langle \mathbf{w}, \Psi(X_i,\bar{Y}_i) \rangle + \Delta(\bar{Y}_i,Y_i) \right) - \langle \mathbf{w}, \Psi(X_i,Y_i) \rangle, \f]
 * which is a convex upper bound on \f$ \Delta(\bar{Y}_i,Y_i) \f$, a customizable function defining the loss associated with predicting \f$ \bar{Y}_i \f$ 
 * when the true label is \f$ Y_i \f$.  
 * 
 * The main thing a person implementing a customized structured SVM learning method needs to implement is an algorithm that efficiently solves
 * \f[ \bar{Y}_i = \arg\max_{Y} \left( \langle \mathbf{w}, \Psi(X_i,Y) \rangle + \Delta(Y,Y_i) \right) \f]
 * Specifically, one should use the following procedure (see Examples tab for examples):
 *  -# Create your own class StructuredLabelCustom which inherits from StructuredLabel.  This is
 *     used to read and write labels y from file and store any custom data for a structured label
 *  -# Create your own class StructuredDataCustom which inherits from StructuredData.  This is
 *     used to read and examples x from file and store any custom data 
 *  -# Create your own class StructuredSVMCustom which inherits from StructuredSVM and defines the
 *     methods Psi(), Inference(), Loss(), Load(), Save(), NewStructuredLabel(), NewStructuredData()
 *  -# Optionally, to allow interactive classification of new test examples and dynamically adding
 *     training examples while training is in progress via commands over the network
 *     (see StructuredLearnerRpc for info on the network protocol), add code 
 *     - int main(int argc, char **argv) { StructuredLearnerRpc v(new StructuredSVMCustom); v.main(argc, argv); }    
 */

/**
 * @class StructuredData
 * @brief Stores data x (e.g., feature data) for a training example (x,y)
 *
 * A person implementing their own custom structured SVM should extend this class, adding
 * custom data storage, and implementing methods read() and write(),
 * for loading and saving a StructuredData object, respectively. 
 */
class StructuredData {
public:
  virtual ~StructuredData() {};

  /**
   * @brief Reads a StructuredData object from a string encoding of the data. 
   * 
   * This method will be called either 
   *   -# when reading a dataset from file, or
   *   -# when adding a new training example via requests over the network
   *
   * @param x A JSON encoding of the StructuredData object
   * @param s A StructuredSVM object defining the structural model parameters
   * @return A pointer to a location in str coming after the last parsed character
   *  in str.  This is used when multiple StructuredData or StructuredLabel objects
   *  are encoded sequentially into the same string.
   */
  virtual bool load(const Json::Value &x, StructuredSVM *s) = 0;

  /**
   * @brief Writes a StructuredData object into a string encoding of the data
   *
   * This method will be called either 
   *   -# when saving a dataset to a file, or
   * @param s A StructuredSVM object defining the structural model parameters
   * @return A JSON encoding of the StructuredData object
   */
  virtual Json::Value save(StructuredSVM *s) = 0;
};

/**
 * @class StructuredLabel
 *
 * @brief Stores a structured label y (e.g., class labels, part locations, segmentation, etc.) for a training example (x,y).  
 *
 * A person implementing their own custom structured SVM should extend 
 * this class, adding their own custom data variables and implementing methods read() and write(),
 * for loading and saving a StructuredLabel object, respectively.  
 */
class StructuredLabel {
 protected:
  StructuredData *x;  /**< A pointer to the data associated with this label */
public:

  /**
   * @brief Create a new structured label
   * @param x A pointer to the data associated with this label
   */
  StructuredLabel(StructuredData *x) { this->x = x; };

  /**
   * @brief Get a pointer to the data associated with this label
   */
  StructuredData *GetData() { return x; }

  virtual ~StructuredLabel() {};

  /**
   * @brief Reads a StructuredLabel object from a JSON encoding of the data.  
   *
   * This method will be called either 
   *   -# when reading a dataset from file, or
   *   -# when adding a new training example via requests over the network.  
   * @param x A JSON encoding of this StructuredLabel
   * @param s A StructuredSVM object defining the structural model parameters
   * @return True if this label was parsed correctly
   */
  virtual bool load(const Json::Value &x, StructuredSVM *s) = 0;

  /**
   * @brief Writes a StructuredLabel object into a JSON encoding of the data
   * This method will be called either 1) when saving a dataset to a file, or
   *   2) when sending this example to a client via requests over the network
   * @param s A StructuredSVM object defining the structural model parameters
   * @return a JSON encoding of the data
   */
  virtual Json::Value save(StructuredSVM *s) = 0;

  virtual bool IsValid() { return true; }
};


/**
 * @class StructuredSVM
 *
 * @brief A class for structured SVM learning and classification.  Supports online learning
 * and interactive classification
 *
 * A person implementing a customized structured SVM should extend this class and then 
 * define overriden methods for Psi(), Inference(), Loss(), Load(), Save(), 
 * NewStructuredLabel(), NewStructuredData().  One can optionally define custom methods for a 
 * constructor/destructor and methods for OnFinishedIteration(), LoadDataset(), SaveDataset().
 *
 * In summary, the steps for creating your own custom structured learning method are (see Examples tab for examples):
 *  -# Create your own class StructuredLabelCustom which inherits from StructuredLabel.  This is
 *     used to read and write labels y from file and store any custom data for a structured label
 *  -# Create your own class StructuredDataCustom which inherits from StructuredData.  This is
 *     used to read and examples x from file and store any custom data 
 *  -# Create your own class StructuredSVMCustom which inherits from StructuredSVM and defines the
 *     routines mentioned above
 *  -# Optionally, to allow interactive classification of new test examples and dynamically adding
 *     training examples while training is in progress via commands over the network
 *     (see StructuredLearnerRpc for info on the network protocol), add code 
 *     - int main(int argc, char **argv) { StructuredLearnerRpc v(new StructuredSVMCustom); v.main(argc, argv); }    
 */
class StructuredSVM {
public:
  StructuredSVM();
  virtual ~StructuredSVM();

  /**
   * @brief Extract features from a training example x with respect to label y
   * 
   * @param x The data for this training example
   * @param y The label at which to extract features
   * @return A vector of features
   */
  virtual SparseVector Psi(StructuredData *x, StructuredLabel *y) = 0;

  /**
   * @brief Solves an inference problem.  
   *
   * There are 3 different ways of invoking this function: 
   * -# A regular inference or classification problem selects the highest scoring label:
   *   Inference(x, ybar, w) solves
   *   \f[ 
   *      \bar{Y} = \arg\max_Y \mathbf{w} \cdot \Psi(X,Y) 
   *   \f]
   * -# During training, the label that is the most violated contraint is
   *   Inference(x, ybar, w, NULL, y_gt) solves
   *   \f[ 
   *      \bar{Y} = \arg\max_Y \mathbf{w} \cdot \Psi(X,Y) + \Delta(Y,Y_{\textrm{gt}})
   *   \f]
   * -# During interactive labeling
   *   Inference(x, ybar, w, y_partial) solves
   *   \f[ 
   * \bar{Y} = \arg\max_Y \mathbf{w} \cdot \Psi(X,Y)\\
   *   \f]
   *   \f[ 
   *    \ \ \ \mathrm{s.t.\ }Y\mathrm{\ is\ consistent\ with\ } Y_{\textrm{partial}}
   *   \f]
   *   where y_partial is a user-specified assignment to some of the variables in ybar
   * 
   * @param x The data for this training example
   * @param ybar The returned predicted label with the highest score (this is modified by the function)
   * @param w A vector of model weights
   * @param y_partial An optional partial assignment to ybar, which constrains which labels are possible
   * @param y_gt The ground truth label of x, which is used only during training when finding the most violated label
   * @return The score of the predicted label (which includes the loss for option 2)
   */
  virtual void saveBoutFeatures(StructuredDataset *dataset, const char *filename, bool sphered=true, bool addRandBouts=true, int window_size = 32){}
  virtual double Inference(StructuredData *x, StructuredLabel *ybar, SparseVector *w, 
			   StructuredLabel *y_partial=NULL, StructuredLabel *y_gt=NULL, double w_scale = 1) = 0;


  /**
   * @param x The data for this training example
   * @param w A vector of model weights
   * @param y_gt The ground truth label of x, which is used only during training when finding the most violated label
   * @param psi_y_gt The features extracted at the ground truth label of x
   * @param set A set of samples that should be extracted by this function.  For correctness, it should include
   *            the same ybar sample found by Inference(x, ybar, w, NULL, y_gt)
   * @return The score of the predicted label (which includes the loss for option 2)
   */
  virtual double ImportanceSample(StructuredData *x, SparseVector *w, StructuredLabel *y_gt, struct _SVM_cached_sample_set *set, double w_scale = 1);

  /**
   * @brief Computes the loss associated with predicting y_pred when the true label is y_gt
   * @param y_gt The true label
   * @param y_pred The predicted label
   * @return The loss associated with predicting y_pred when the true label is y_gt
   */
  virtual double Loss(StructuredLabel *y_gt, StructuredLabel *y_pred) = 0;

  /**
   * @brief Read all info for a structured SVM.  This may include variables for the model definition
   *  (e.g., things like the number of possible classes or features) as well as the learned model
   *  weights, learning parameters C and eps, etc.
   * @param root A JSON object from which to read all values from
   * @return true if the model was loaded successfully
   */
  virtual bool Load(const Json::Value &root) = 0;

  /**
   * @brief Save all info for a structured SVM.  This may include variables for the model definition
   *  (e.g., things like the number of possible classes or features).  
   * @return A JSON encoding of this structured SVM
   */
  virtual Json::Value Save() = 0;


  /**
   * @brief Create a new empty label.  Typically, this just calls new StructuredLabelCustom, where StructuredLabelCustom
   *   is the API user's custom class implementing a label y
   */
  virtual StructuredLabel *NewStructuredLabel(StructuredData *x) = 0;

  /**
   * @brief Create a new empty data example.  Typically, this just calls new StructuredDataCustom, where StructuredDataCustom
   *   is the API user's custom class implementing an example x
   */
  virtual StructuredData *NewStructuredData() = 0;

  /**
   * @brief Function invoked by online learning algorithm after it finishes iterating over an example (x,y).  This is
   *   an optional function which is useful if the entire training set can't fit in memory, and we want to clear
   *   allocated memory stored in x and y
   * @param x The data example we just processed
   * @param y The label we just processed
   */
  virtual void OnFinishedIteration(StructuredData *x, StructuredLabel *y, SparseVector *w = NULL, StructuredLabel *ybar=NULL) {}

  /**
   * @brief Load a dataset of examples from file.  By default, this assumes each line in the file will correspond
   *   to one example in the format "x y", where x and y are in the same format as StructuredData::read()
   *   and StructuredLabel::read(), but the API user can optionally override this function
   * @param fname The filename from which to read the dataset
   */
  virtual StructuredDataset *LoadDataset(const char *fname);

  /**
   * @brief Save a dataset of examples to file.  By default, this assumes each line in the file will correspond
   *   to one example in the format "x y", where x and y are in the same format as StructuredData::read()
   *   and StructuredLabel::read(), but the API user can optionally override this function
   * @param d The dataset to save
   * @param fname The filename in which to save the dataset
   * @param start_from If non-zero, saves only training examples beginning at this index, appending the output file
   */
  virtual bool SaveDataset(StructuredDataset *d, const char *fname, int start_from = 0);

  /**
   * @brief Read all info for a structured SVM.  This may include variables for the model definition
   *  (e.g., things like the number of possible classes or features) as well as the learned model
   *  weights, learning parameters C and eps, etc.  Typically, one should not override this function.  
   * @param fname The filename from which to load the structured SVM model
   * @param loadFull If true, loads data that is needed to continue online learning again
   * @return true if the model was loaded successfully
   */
  virtual bool Load(const char *fname, bool loadFull=false);

  /**
   * @brief Save all info for a structured SVM.  This may include variables for the model definition
   *  (e.g., things like the number of possible classes or features) as well as the learned model
   *  weights, learning parameters C and eps, etc.  Typically, one should not override this function
   * @param fname The filename from which to save the structured SVM model
   * @param saveFull If true, saves data that is needed to continue online learning again
   * @return true if the model was saved successfully
   */
  virtual bool Save(const char *fname, bool saveFull=false);

  
  /**
   * @brief the number of columns for HTML visualizations of a dataset, where each column contains one training example
   */
  virtual int NumHTMLColumns() { return 4; }
  
  /**
   * @brief Create an HTML visualization of this training example
   * @param htmlDir A directory where html for the dataset visualization is stored.  One can save things like images into this directory, which can be referenced by the returned html code
   * @param extraInfo custom info associated with this label (e.g., the classification score)
   * @return a string allocated using malloc containing html code to visualize this label
   */
  virtual char *VisualizeExample(const char *htmlDir, StructuredExample *ex, const char *extraInfo=NULL);

  virtual void OnFinishedPassThroughTrainset() {} 

  /************************** Routines used for training **************************/
public:

  /**
   * @brief Learn the structured model weights.  While this function is running, a client is free to
   *   add new examples or use the current model parameters to classify examples
   * @param modelfile If non-null, file where to store the learned model 
   * @param runForever If true, keep training running forever (because the client might add more
   * training examples).  Otherwise, stop training once the optimization error is below eps
   */
  void Train(const char *modelfile=NULL, bool runForever=false, const char *initial_sample_set=NULL);  

  /**
   * @brief Evaluate the structured model on a testset
   * @param testfile The filename of the dataset, in the format of LoadDataset()
   * @param predictionsFile A file where the predicted labels will be written, with lines in the format:
\htmlonly <div style="padding: 0.5em 1em; border-top: 1px solid #ddd; border-bottom: 1px solid #ddd; background-color: #eaeafa;">
y_predicted y_ground_truth loss score_prediction score_ground_truth
</div> \endhtmlonly
   *  where y_predicted and y_ground_truth are in the format of StructuredLabel::read()
   * @return The average classification loss
   */
  VFLOAT Test(const char *testfile, const char *predictionsFile=NULL, const char *htmlDir=NULL, double *svm_err=NULL);

  /**
   * @brief Add a new training example.  The data is copied (as opposed to storing pointers)
   * @param x The data for the example
   * @param y The label for the example
   * @return The index of the newly added training example
   */
  int AddExample(StructuredData *x, StructuredLabel *y);  

  void RelabelExample(StructuredExample *ex, StructuredLabel *y);

  /**
   * @brief Load a training set from file in the format of LoadDataset()
   * @param fname The filename of the training set
   */
  void LoadTrainset(const char *fname);

  /**
   * @brief Get the current model weights w.  
   * @return The vector of weights.  It is dynamically allocated and should be freed using delete
   * @param lock if true, calls Lock() to synchronize access to w 
   */
  SparseVector *GetCurrentWeights(bool lock=true);

  /**
   * @brief Acquire the lock for this structured svm, to make access to common datastructures thread safe
   */
  void Lock() { omp_set_lock(&my_lock); }

  /**
   * @brief Release the lock for this structured svm, to make access to common datastructures thread safe
   */
  void Unlock() { omp_unset_lock(&my_lock); }

  /**
   * @brief Signal the learner (Iniitially started by calling Train()) to stop running
   */
  void Shutdown() { finished = true; }

  /**
   * @brief Change the lambda regularization parameter (where lambda=1/C).  This can be called as Train()
   * is in progress if the optimization method is one of the _WITH_CACHE methods.
   * @param l The parameter lambda
   * @param num_iter Number of times to iterate over the cache set to adjust the learned model parameters 
   * with the new lambda taken into account
   */
  void SetLambda(double l, int num_iter=0);

  /**
   * @brief Change the regularization parameter C (where lambda=1/C).  This can be called as Train()
   * is in progress if the optimization method is one of the _WITH_CACHE methods.
   * @param c The parameter C
   * @param num_iter Number of times to iterate over the cache set to adjust the learned model parameters 
   * with the new C taken into account
   */
  void SetC(double c, int num_iter=0) { SetLambda(1.0/c, num_iter); }

  /**
   * @brief Set the optimization method used for training
   * @param m The optimization method.  Options are:
   *  - SPO_SGD: stochastic gradient descent: \f$ \mathbf{w}^t = \mathbf{w}^{t-1} - \gamma (\lambda \mathbf{w}^{t-1} + \nabla \ell_i(\mathbf{w}^{t-1})) \f$,   \f$ \gamma=\frac{1}{\lambda*t} \f$
   *  - SPO_SGD_PEGASOS: do SGD step, then ensure that \f$ \|\mathbf{w}^t\|^2 \le \frac{1}{\lambda} \f$   (downscale \f$ \mathbf{w}^t \f$ if this doesn't hold) 
   *  - SPO_DUAL_UPDATE: \f$ \mathbf{w}^t = \frac{t-1}{t} \mathbf{w}^{t-1} - \alpha \nabla \ell_i(\mathbf{w}^{t-1}) \f$,   where \f$ \alpha=\min \left(1,\max \left(0,\frac{\lambda (t-1)\ \langle \mathbf{w}^{t-1}, \mathbf{v}_t^{\bar{Y}_t} \rangle  + \lambda t \ \Delta(\bar{Y}_t,Y_t)}{\| \mathbf{v}_t^{\bar{Y}_t} \|^2  }\right) \right) \f$ is chosen to maximize dual objective
   *  - SPO_DUAL_UPDATE_WITH_CACHE: Same as SPO_DUAL_UPDATE, but also run dual update steps in a background process on labels from earlier iterations
   *  - SPO_DUAL_MULTI_SAMPLE_UPDATE: Sample multiple labels per iteration (instead of just getting the "most violated constraint"), then optimize parameters jointly
   *  - SPO_DUAL_MULTI_SAMPLE_UPDATE_WITH_CACHE: Same as SPO_DUAL_MULTI_SAMPLE_UPDATE, but also run multi-sample dual update steps in a background process
   */ 
  void SetMethod(StructuredPredictionOptimizationMethod m) { method = m; }

  /**
   * @brief Set the target training accuracy epsilon.  Training ends when the training error is within epsilon of the minimum achievable training error
   * @param epsilon The target training accuracy epsilon
   */
  void SetEpsilon(double epsilon) { eps = epsilon; }
  /**
   * @brief Get the target training accuracy epsilon.
   */
  double GetEpsilon() { return eps; }

  /**
   * @brief Get the lambda regularization parameter (where lambda=1/C).  
   */
  double GetLambda() { return lambda; }


  /**
   * @brief Save the training set to the same location read by LoadTrainset()
   * @param start_from If greater than 0, begins saving from exampple number start_from, appending
   * the training set file.  This is useful if adding new examples gradually
   */
  bool SaveTrainingSet(int start_from);

  /**
   * @brief Get statistics useful for plotting the progression of different types of 
   * training error as a function of training computation time
   */
  void GetStatisticsByIteration(int ave, long *tt, double **gen_err_buff, double **emp_err_buff, double **model_err_buff, 
				double **reg_err_buff, double **loss_buff, double **time_buff);

  /**
   * @brief Get statistics useful for plotting the progression of different types of 
   * training error as a function of number of training examples
   */
  void GetStatisticsByExample(int ave, long *nn, double **gen_err_buff, double **emp_err_buff, double **model_err_buff, 
			      double **reg_err_buff, double **loss_buff);

  /**
   * @brief Get the dimensionality of the structured feature space Psi(x,y)
   */
  int GetSizePsi() { return sizePsi; }

  bool *GetUseWeights() { return useWeights; }

  /**
   * @brief Create a new example.  It allocates and copies a new version of x and y, rather than storing a pointer to those objects
   * @param x the data
   * @param y the label
   */
  StructuredExample *CopyExample(StructuredData *x, StructuredLabel *y, StructuredLabel *y_latent=NULL);

  /**
   * @brief get a pointer to the training set object
   */
  StructuredDataset *GetTrainset() { return trainset; }

  /**
   * @brief set the training set object
   */
  void SetTrainset(StructuredDataset *t);

  /**
   * @brief Set whether or not to train in a multi-threaded algorithm
   * @param r if non-zero, train using a multi-threaded algorithm where r is the number of threads (or setting to -1 auto-detects the number of CPU cores)
   */
  void RunMultiThreaded(int r) { runMultiThreaded = r; }

  /**
   * @brief Set verbosity controlling how much debug information is outputted
   * @param level The verbosity level
   */
  void SetVerbosity(int level) { debugLevel = level; }

  /**
   * @brief Get the elapsed training time
   */
  double GetElapsedTime() { double t = (get_runtime() - start_time + base_time); if(num_thr > 1) { t /= num_thr; } return t; }

  /**
   * @brief Get the number of training iterations
   */
  long GetNumIterations() { return t; }

  /**
   * @brief Get the number of training examples
   */
  long GetNumExamples() { return n; }

  const char *GetModelFile() { return modelfile; }
  void SetValidationFile(const char *v) { validationfile = StringCopy(v); }
  void SetSumWScale(double sum_w_scale_new);

  void DumpModelIfNecessary(const char *modelfile, bool force=false);
  void DumpWeights(double start, double factor) { dumpModelStartTime=start;  dumpModelFactor=factor; }

  void SetWindow(int w) { window = w; }
  void SetMinItersBeforeNewExample(int m) { minItersBeforeNewExample = m; }
  void SetModelId(int m) { modelId = m; }
  void SetIterId(int i) { iterId = i; }
  void VisualizeDataset(StructuredDataset *dataset, const char *htmlDir, int maxExamples=-1);
  void SetWeights(double *d) { if (sum_w) { delete sum_w; } sum_w = new SparseVector(d, sizePsi, true);  *sum_w *= sum_w_scale; }
  void InferLatentValues(StructuredDataset *d);

  void SaveCachedExamples(const char *output_name, bool saveFull=true);
  void LoadCachedExamples(const char *fname, bool loadFull=true);
  bool IsTesting() { return isTesting; }

 protected:

  SparseVector *sum_w; /**< the unnormalized learned model weights w^t = sum_w*(t*lambda) */ 
  int sizePsi;         /**< maximum number of weights in w */
  double eps;          /**< precision for which to solve optimization problem */
  double C;            /**< regularization parameter */
  StructuredPredictionOptimizationMethod method;  /**< optimization method */
  int window;  /**< default window size used for computing statistics for decomposing the test error */

  bool *regularize;  /**< if non-null, a sizePsi array where a value of true means we regularize this weight entry */
  bool *learnWeights;  /**< if non-null, a sizePsi array where a value of true means we learn this weight entry */
  int *weightConstraints;  /**< if non-null, a sizePsi array where a value of 1 means a weight must be >= 0, -1 means <= 0, 0 means anything */
  bool *useWeights;
  bool canScaleW;   /**< If true, Inference() should support scaling w.  This results in a speed improvement for some models, by avoiding normalizing w */

  VFLOAT lambda;       /**< regularization constant */
  double sum_dual;  /**< The dual objective D_t(alpha_1,...alpha_t) over each training iteration */
  double sum_alpha_loss;
  double sum_w_sqr; /**< (lambda*t*w)^2 */
  double regularization_error;      /**< Current regularization error  */ 
  double sum_w_scale;
  int maxIters;

  int numMineHardNegativesRound;
  int mineHardNegativesRound;
  int numHardNegativesPerExample;
  bool useFixedSampleSet;
  bool mergeSamples;
  bool updateFromCacheThread;
  int numCacheUpdatesPerIteration;
  double maxLoss;
  int numMultiSampleIterations;
  int maxCachedSamplesPerExample;
  bool savingCachedExamples;
  bool relabelingExample;

  int modelId;
  int iterId;

protected:
  int debugLevel; /**< debug level defining how many print statements are outputted by the learner */ 
  long t;  /**< number of iterations run so far */ 
  long n;  /**< number of examples that have been processed in at least one iteration */ 
  char *modelfile;   /**< name of the file to save the model to */ 
  char *trainfile;   /**< name of the file for the training set */ 
  char *validationfile;   
  int runMultiThreaded;  /**< if non-zero, train using a multi-threaded algorithm where runMultiThreaded is the number of threads (or setting to -1 auto-detects the number of CPU cores) */ 
  int max_samples;


  virtual void MultiSampleUpdate(struct _SVM_cached_sample_set *set, StructuredExample *ex, int R=1);

 private:
  omp_lock_t my_lock;
  

  
  /************************ Variables used for online learning ******************************/
  double *u_i_buff;     // memory buffer used for a non-sparse version of u_i, u_i = \sum_y alpha_{i,y} (psi(x_i,ybar)-psi(x_i,y_i))
  bool runForever;
  bool hasConverged;
  long curr;  /**< the next example to process */ 
  long just_added;  /**< if >= 0, a new example was just added and we may want to focus on it */ 
  int M;   /**< number of online updates per example (determined automatically) */ 
  bool finished;  /**< If set to true, then the Train() threads will all finish up after the next iteration */ 
  int num_thr;    /**< number of threads */ 

  bool isMultiSample;
  double dumpModelStartTime;
  double dumpModelFactor;
  int numModelDumps;
  bool trainLatent;

  int alloc_n, alloc_t;
  bool cache_old_examples;
  int *ex_num_iters;  // an array of size sample.n, where each entry stores the number of iterations each training example has been processed */ 
  int *ex_first_iter; /**<  an array of size sample.n, where each entry stores the iteration (value of t) when each training example was first iterated over; */
  int last_example;  /**< The index of the last training example iterated over during training */

  double sum_generalization_error;  /**< Total error measured on new examples as we stream them in and before processing them. */
  double sum_iter_error;            /**< Total error measured on unseen labels.  This is different than sum_generalization_error if we have 
				      processed each example more than once (it is effectively the sum loss associated with each call 
				      to Inference(x_i, ybar_i, w_t, NULL, y_i) */
  double sum_iter_error_window, sum_generalization_error_window;
  double *generalization_errors_by_n;
  double *generalization_errors_by_t, *iter_errors_by_t, *sum_dual_by_t, *regularization_errors_by_t, *losses_by_t, *elapsed_time_by_t;
  long *iter_examples;
  double base_time, start_time;
  
  int *examples_by_iteration_number; /**< An array of size M, where each entry stores an index to the start of a linked list (implemented in 
					examples_by_iteration_next_ind).  Each linked list stores the indices of all examples that have been iterated over
					i times */
  int *examples_by_iteration_next_ind;  /**< An array of size sample.n.  For each example, contains the index of the next example in a linked list of examples
				       that have been iterated over the same number of times as this example */ 
  int *examples_by_iteration_prev_ind;  /**< An array of size sample.n.  For each example, contains the index of the previous example in a linked list of examples
				       that have been iterated over the same number of times as this example */ 
  long minItersBeforeNewExample;
  int currMinIterByExample;

 protected:
  StructuredDataset *trainset;
  bool keepAllEvictedLabels;
  bool isTesting;
  
  void SVM_cached_sample_set_compute_features(struct _SVM_cached_sample_set *set, StructuredExample *ex);

 private:

  int ChooseNextExample();
  void UpdateWeights(struct _SVM_cached_sample_set *ex, int iterInd);
  void UpdateStatistics(struct _SVM_cached_sample_set *set, int iter);
  void UpdateExampleIterationQueue(struct _SVM_cached_sample_set *set, int iter);
  void CheckConvergence();

  long UpdateWeightsAddStatisticsBefore(struct _SVM_cached_sample_set *ex, int iterInd);
  void UpdateWeightsAddStatisticsAfter(struct _SVM_cached_sample_set *ex, int iterInd, long tt);
  void RecomputeWeights(bool full=true);
  void OptimizeAllConstraints(int num_iter);
  bool SaveOnlineData(const char *fname);
  bool LoadOnlineData(const char *fname);
  void CreateTrainingExampleQueues(int ind);
  void SingleSampleUpdate(struct _SVM_cached_sample_set *set, bool useSmartStepSize);
  VFLOAT Test(StructuredDataset *testset, const char *predictionsFile, const char *htmlDir, double *svm_err, bool getLock);

  void ExtractSampleSet(int num_samples, bool augment);
  void ConvertCachedExamplesToBinaryTrainingSet();
  void CondenseSamples(struct _SVM_cached_sample_set *set);
  void UncondenseSamples(struct _SVM_cached_sample_set *set);
  void UpdateFromCache(bool lock=true, int *num=NULL, int i=-1);
  void TrainMain(const char *modelfile=NULL, bool runForever=false, const char *initial_sample_set=NULL); 
  void TrainBinary(const char *modelfile=NULL, bool runForever=false, const char *initial_sample_set=NULL); 
};



/**
 * @class StructuredExample
 * @brief Simple class implementing a structured example (x,y).  This is just a StructuredData and StructuredLabel object
 */
class StructuredExample {
 public:
  StructuredData *x;   /**< The data for this example */
  StructuredLabel *y;  /**< The label for this example */

  StructuredLabel *y_latent;  /**< The label for this example, without latent values filled.  If this is non-null, y will be modified over time, such that at any instant it will have been assigned a particular value for all latent variables */

  struct _SVM_cached_sample_set *set;   /**< For some optimization methods, cache the set of extracted samples for this example */

  StructuredExample();
  ~StructuredExample(); 
};


/**
 * @class StructuredDataset
 * @brief Simple class implementing a dataset of examples.  This is just an array of StructuredExample
 */
class StructuredDataset {
 public:
  StructuredExample **examples;  /**< An array of num_examples examples */
  int num_examples;   /**< The number of examples in the dataset */

  StructuredDataset();
  ~StructuredDataset(); 

  /**
   * @brief Append a new example to the array of dataset examples
   * @param e The example to be added.  This function does not copy e; it just stores a pointer
   */
  void AddExample(StructuredExample *e); 

  /**
   * @brief Randomly permute all examples in this dataset
   */
  void Randomize();
};


/**
 * @struct _SVM_cached_sample
 *
 * @brief Helper struct used by StructuredSVMOnlineLearner.  Encodes data associated with a call to Inference(x_i, ybar_i, w_t, NULL, y_i) 
 */
typedef struct _SVM_cached_sample {
  StructuredLabel *ybar;       /**< the label \f$ \bar{Y}_i \f$ */
  SparseVector *psi;    /**< feature vector of the ith example at ybar: \f$ \Psi(X,\bar{Y}_i) \f$ */
  VFLOAT loss;       /**< the loss: \f$ \Delta(\bar{Y}_i, Y_i) \f$ */
  VFLOAT alpha;      /**< dual parameter: \f$ \alpha_i^{\bar{Y}_i} \f$ for this sample */
  VFLOAT sqr;        /**< Squared version of dpsi: \f$ \| \Psi(X,\bar{Y}_i) - \Psi(X,Y_i) \|^2 \f$ */
  VFLOAT slack;      /**< \f$ \langle w, \Psi(X,\bar{Y}_i) - \Psi(X,Y_i) \rangle + \Delta(\bar{Y}_i, Y_i) \f$ */
  VFLOAT dot_psi_gt_psi;  /**< \f$ \langle \Psi(X,Y_i), \Psi(X,\bar{Y}_i) \rangle \f$ */
  VFLOAT dot_w;      /**< \f$ \langle w, \Psi(X,\bar{Y}_i) - \Psi(X,Y_i) \rangle  \f$ */
} SVM_cached_sample;

/**
 * @struct _SVM_cached_sample_set
 *
 * @brief Holds a set of sample labels for the same example
 *
 * Helper struct used by StructuredSVMOnlineLearner.  Encodes data associated with a call to Inference(x_i, ybar_i, w_t, NULL, y_i). 
 */
typedef struct _SVM_cached_sample_set {
  SVM_cached_sample *samples;  /**< A set of samples: \f$ \bar{Y}_i^1,...,\bar{Y}_i^L \f$ */
  int num_samples;             /**< Number of samples L */
  int i;                       /**< training example index (index into sample.examples) */
  double alpha;                /**< sum_i(samples[i].alpha) */
  double loss;                 /**< sum_i(samples[i].alpha*samples[i].loss) */
  double slack_before;         /**< (<w_{t-1},samples[0].dpsi>) + loss */
  double slack_after;          /**< sum_i(samples[i].alpha*(<w_t,samples[i].dpsi>) + loss*/
  double score_gt;             /**< (<w_{t-1},y_i>) */
  SparseVector *psi_gt;        /**< Psi(x_i,y_i) */
  double psi_gt_sqr;           /**< |Psi(x_i,y_i)|^2 */
  StructuredLabel *ybar;
  double dot_sum_w_psi_gt;
  double *alphas;

  SparseVector *u_i;           /**< u_i = \sum_ybar alpha_{i,ybar} Psi(x_i,ybar) */
  VFLOAT D_i;                  /**< D_i = \sum_ybar alpha_{i,ybar} loss(y_i,ybar) */
  VFLOAT dot_u_psi_gt;         /**< dot_u_psi_gt = <u_i,y_i> */
  double u_i_sqr;              /**< u_i_sqr = <u_i,u_i> */
  double drift_bits;           /**< Used to figure out how frequently to recompute u_i_sqr to avoid numerical precision errors */
  bool lock;
  int numIters;
  double sumSlack;

  StructuredLabel **evicted_labels;
  int num_evicted_labels;
} SVM_cached_sample_set;


///@cond
void free_SVM_cached_sample(SVM_cached_sample *s);
void read_SVM_cached_sample(SVM_cached_sample *s, FILE *fin, StructuredSVM *svm, bool readFull);
void write_SVM_cached_sample(SVM_cached_sample *s, FILE *fout, StructuredSVM *svm, bool writeFull);
SVM_cached_sample_set *new_SVM_cached_sample_set(int i, SparseVector *psi_gt=NULL);
void free_SVM_cached_sample_set(SVM_cached_sample_set *s);
SVM_cached_sample_set *read_SVM_cached_sample_set(FILE *fin, StructuredSVM *svm, StructuredData *x, bool readFull);
void write_SVM_cached_sample_set(SVM_cached_sample_set *s, FILE *fout, StructuredSVM *svm);
void write_SVM_cached_sample_set(SVM_cached_sample_set *s, FILE *fout, StructuredSVM *svm, bool writeFull);
SVM_cached_sample *SVM_cached_sample_set_add_sample(SVM_cached_sample_set *s, StructuredLabel *ybar);
void clear_SVM_cached_sample(SVM_cached_sample *s);
int SVM_cached_sample_cmp(const void * a, const void * b);
int SVM_cached_sample_set_alpha_cmp(const void *a, const void *b);
int SVM_cached_sample_set_ave_slack_cmp(const void *a, const void *b);
///@endcond



#endif
