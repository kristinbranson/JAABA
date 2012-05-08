#ifndef __STRUCTURED_SVM_H
#define __STRUCTURED_SVM_H

#include "sparse_vector.h"

#include <omp.h>
#include <time.h>


typedef enum {
  SPO_CUTTING_PLANE,   /**< SVM^struct */

  // Online Algorithms, which process one training example i per iteration:
  SPO_SGD,             /**< stochastic gradient descent: w^t = w^{t-1} - step_size*(w^{t-1} + grad_i(w^{t-1})),   step_size=1/(lambda*t) */
  SPO_SGD_PEGASOS,     /**< Do SGD step, then ensure that ||w^t||^2 <= 1/lambda   (downscale w^t if this doesn't hold) */
  SPO_DUAL_UPDATE,     /**< w^t = w^{t-1}(t-1)/t - alpha*grad_i(w^{t-1}),   where alpha is chosen to maximize dual objective */
  SPO_DUAL_UPDATE_WITH_CACHE,     /**< Same as SPO_DUAL_UPDATE, but also run dual update steps in a background process on labels from earlier iterations */
  SPO_DUAL_MULTI_SAMPLE_UPDATE,   /**< Sample multiple labels per iteration (instead of just getting the "most violated constraint"), then optimize parameters jointly */
  SPO_DUAL_MULTI_SAMPLE_UPDATE_WITH_CACHE,  /**< Same as SPO_DUAL_MULTI_SAMPLE_UPDATE, but also run multi-sample dual update steps in a background process  */
} StructuredPredictionOptimizationMethod;


class StructuredExample;
class StructuredDataset;
class StructuredSVM;
struct _SVM_cached_sample;
struct _SVM_cached_sample_set;


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
   * @param str String where the encoding of the StructuredData object will be written
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
  StructuredData *x;
public:
  StructuredLabel(StructuredData *x) { this->x = x; };
  StructuredData *GetData() { return x; }
  virtual ~StructuredLabel() {};

  /**
   * @brief Reads a StructuredLabel object from a string encoding of the data.  
   *
   * This method will be called either 
   *   -# when reading a dataset from file, or
   *   -# when adding a new training example via requests over the network.  
   * @param str String storing the encoding of the StructuredLabel object
   * @param s A StructuredSVM object defining the structural model parameters
   * @return A pointer to a location in str coming after the last parsed character
   *  in str.  This is used when multiple StructuredData or StructuredLabel objects
   *  are encoded sequentially into the same string.
   */
  virtual bool load(const Json::Value &x, StructuredSVM *s) = 0;

  /**
   * @brief Writes a StructuredLabel object into a string encoding of the data
   * This method will be called either 1) when saving a dataset to a file, or
   *   2) when sending this example to a client via requests over the network
   * @param str String where the encoding of the StructuredLabel object will be written
   * @param s A StructuredSVM object defining the structural model parameters
   * @return A pointer to a location in str coming after the last written character
   *  in str.  This is used when multiple StructuredData or StructuredLabel objects
   *  are encoded sequentially into the same string.
   */
  virtual Json::Value save(StructuredSVM *s) = 0;
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
   *      \bar{y} = \arg\max_y w \cdot \Psi(x,y) 
   *   \f]
   * -# During training, the label that is the most violated contraint is
   *   Inference(x, ybar, w, NULL, y_gt) solves
   *   \f[ 
   *      \bar{y} = \arg\max_y w \cdot \Psi(x,y) + Loss(y_{gt},y)
   *   \f]
   * -# During interactive labeling
   *   Inference(x, ybar, w, y_partial) solves
   *   \f[ 
   * \bar{y} = \arg\max_y w \cdot \Psi(x,y)\\
   *   \f]
   *   \f[ 
   *    \ \ \ \mathrm{s.t.\ }y\mathrm{\ is\ consistent\ with\ } y_{partial}
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
  virtual void saveBoutFeatures(StructuredDataset *dataset, const char *filename, bool sphered=true, bool addRandBouts=true){}
  virtual double Inference(StructuredData *x, StructuredLabel *ybar, SparseVector *w, 
			   StructuredLabel *y_partial=NULL, StructuredLabel *y_gt=NULL) = 0;

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
  virtual void OnFinishedIteration(StructuredData *x, StructuredLabel *y) {}

  /**
   * @brief Load a dataset of examples from file.  By default, this assumes each line in the file will correspond
   *   to one example in the format "<x> <y>", where <x> and <y> are in the same format as StructuredData::read()
   *   and StructuredLabel::read(), but the API user can optionally override this function
   * @param fname The filename from which to read the dataset
   */
  virtual StructuredDataset *LoadDataset(const char *fname);

  /**
   * @brief Save a dataset of examples to file.  By default, this assumes each line in the file will correspond
   *   to one example in the format "<x> <y>", where <x> and <y> are in the same format as StructuredData::read()
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

  
  /************************** Routines used for training **************************/
public:

  /**
   * @brief Learn the structured model weights.  While this function is running, a client is free to
   *   add new examples or use the current model parameters to classify examples
   * @param modelfile If non-null, file where to store the learned model 
   * @param runForever If true, keep training running forever (because the client might add more
   * training examples).  Otherwise, stop training once the optimization error is below eps
   */
  void Train(const char *modelfile=NULL, bool runForever=false);  

  /**
   * @brief Evaluate the structured model on a testset
   * @param testfile The filename of the dataset, in the format of LoadDataset()
   * @param predictionsFile A file where the predicted labels will be written, with lines in the format:
   *   - <y_predicted> <y_ground_truth> <loss> <score_prediction> <score_ground_truth>
   *  where <y_predicted> and <y_ground_truth> are in the format of StructuredLabel::read()
   * @return The average classification loss
   */
  VFLOAT Test(const char *testfile, const char *predictionsFile=NULL);

  /**
   * @brief Add a new training example.  The data is copied (as opposed to storing pointers)
   * @param x The data for the example
   * @param y The label for the example
   * @return The index of the newly added training example
   */
  int AddExample(StructuredData *x, StructuredLabel *y);  

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
   * @brief Change the featureScale parameter, which is a constant used to scale Psi(x,y) and should
   * typically be left at 1.  This can be called as Train()
   * is in progress if the optimization method is one of the _WITH_CACHE methods.
   * @param fs The parameter featureScale
   * @param num_iter Number of times to iterate over the cache set to adjust the learned model parameters 
   * with the new featureScale taken into account
   */
  void SetFeatureScale(double fs, int num_iter=0);

  /**
   * @brief Set the optimization method used for training
   * @param m The optimization method.  Options are:
   *  - SPO_SGD: stochastic gradient descent: w^t = w^{t-1} - step_size*(w^{t-1} + grad_i(w^{t-1})),   step_size=1/(lambda*t)
   *  - SPO_SGD_PEGASOS: do SGD step, then ensure that ||w^t||^2 <= 1/lambda   (downscale w^t if this doesn't hold)
   *  - SPO_DUAL_UPDATE: w^t = w^{t-1}(t-1)/t - alpha*grad_i(w^{t-1}),   where alpha is chosen to maximize dual objective
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
   * @brief Get the featureScale parameter, which is a constant used to scale Psi(x,y) and should
   * typically be left at 1.  
   */
  double GetFeatureScale() { return featureScale; }

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
  void GetStatisticsByIteration(int ave, long *tt, long *tm, double **gen_err_buff, double **opt_err_buff, double **model_err_buff, 
				double **reg_err_buff, double **train_err_buff, double **test_err_buff, double **time_buff);

  /**
   * @brief Get statistics useful for plotting the progression of different types of 
   * training error as a function of number of training examples
   */
  void GetStatisticsByExample(int ave, long *nn, double **gen_err_buff, double **opt_err_buff, double **model_err_buff, 
			      double **reg_err_buff, double **train_err_buff, double **test_err_buff);

  /**
   * @brief Get the dimensionality of the structured feature space Psi(x,y)
   */
  int GetSizePsi() { return sizePsi; }

  StructuredExample *CopyExample(StructuredData *x, StructuredLabel *y);

  StructuredDataset *GetTrainset() { return trainset; }

 protected:

  SparseVector *sum_w; /**< the unnormalized learned model weights w^t = sum_w*(t*lambda) */ 
  int sizePsi;         /**< maximum number of weights in w */
  double eps;          /**< precision for which to solve optimization problem */
  double C;            /**< regularization parameter */
  double featureScale; /**< factor in which to scale all features */
  StructuredPredictionOptimizationMethod method;  /* optimization method */
  int window;

  bool *regularize;  /**< if non-null, a sizePsi array where a value of true means we regularize this weight entry */
  bool *learnWeights;  /**< if non-null, a sizePsi array where a value of true means we learn this weight entry */
  int *weightConstraints;  /**< if non-null, a sizePsi array where a value of 1 means a weight must be >= 0, -1 means <= 0, 0 means anything */

protected:
  int debugLevel;
  long t;  /**< number of iterations run so far */ 
  long n;  /**< number of examples that have been processed in at least one iteration */ 


 private:
  omp_lock_t my_lock;
  VFLOAT lambda;       /**< regularization constant */
  

  
  /************************ Variables used for online learning ******************************/
  StructuredDataset *trainset;
  bool runForever;
  bool hasConverged;
  long curr;  /**< the next example to process */ 
  long just_added;  /**< if >= 0, a new example was just added and we may want to focus on it */ 
  int M;   /**< number of online updates per example (determined automatically) */ 
  bool finished;  /**< If set to true, then the Train() threads will all finish up after the next iteration */ 
  int num_thr;    /**< number of threads */ 


  int alloc_n, alloc_t;
  bool cache_old_examples;
  struct _SVM_cached_sample_set **cached_examples; // an array of 't' cached examples, one per iteration of the optimization algorithm */ 
  int *ex_num_iters;  // an array of size sample.n, where each entry stores the number of iterations each training example has been processed */ 

  double sum_generalization_error;  /**< Total error measured on unseen examples as they are streamed in */ 
  double sum_iter_error;            /**< Total error measured on unseen labels.  This is different than sum_generalization_error if we have *
				      processed each example more than once (it is effectively the sum loss associated with each call 
				      to find_most_violated_constraint()  */ 
  double sum_model_error;            /**< Approximate lower-bound for model error.  This is measured by computing the error  
				      on ybar after each call to find_most_violated_constraint(), but after the model has been updated
				      with respect to ybar  */ 
  double sum_generalization_error_window, sum_iter_error_window, sum_model_error_window;
  double regularization_error;      /**< Current regularization error  */ 
  double *generalization_errors_by_n, *optimization_errors_by_n, *model_errors_by_n, *regularization_errors_by_n, *losses_by_n;
  double *generalization_errors_by_t, *iter_errors_by_t, *model_errors_by_t, *constraint_errors_by_t, *regularization_errors_by_t, *losses_by_t, *elapsed_time_by_t;
  long *iter_examples;
  long base_time, start_time;
  char *modelfile;
  char *trainfile;
  
  int *examples_by_iteration_number; /**< An array of size M, where each entry stores an index to the start of a linked list (implemented in 
					examples_by_iteration_next_ind).  Each linked list stores the indices of all examples that have been iterated over
					i times */
  int *examples_by_iteration_next_ind;  /**< An array of size sample.n.  For each example, contains the index of the next example in a linked list of examples
				       that have been iterated over the same number of times as this example */ 
  int *examples_by_iteration_prev_ind;  /**< An array of size sample.n.  For each example, contains the index of the previous example in a linked list of examples
				       that have been iterated over the same number of times as this example */ 
  long minItersBeforeNewExample;
  int currMinIterByExample;

 private:
  long GetElapsedTime() { return (long)(time(NULL) - start_time + base_time); }

  int ChooseNextExample();
  double UpdateWeights(struct _SVM_cached_sample_set *ex, int iterInd);
  long UpdateWeightsAddStatisticsBefore(struct _SVM_cached_sample_set *ex, int iterInd, double e);
  void UpdateWeightsAddStatisticsAfter(struct _SVM_cached_sample_set *ex, int iterInd, double e, long tt);
  void RecomputeWeights();
  void OptimizeAllConstraints(int num_iter);
  bool SaveOnlineData(const char *fname);
  bool LoadOnlineData(const char *fname);
  void CreateTrainingExampleQueues(int ind);
};



/**
 * @class StructuredExample
 * @brief Simple class implementing a structured example (x,y).  This is just a StructuredData and StructuredLabel object
 */
class StructuredExample {
 public:
  StructuredData *x;   /**< The data for this example */
  StructuredLabel *y;  /**< The label for this example */

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
  void Randomize();
};


/**
 * @struct _SVM_cached_sample
 *
 * @brief Helper struct used by StructuredSVMOnlineLearner.  Encodes data associated with a call to find_most_violated_constraint()
 */
typedef struct _SVM_cached_sample {
  StructuredLabel *ybar;       /**< label */
  SparseVector *dpsi;    /**< gradient w.r.t. ith example at ybar: psi(ybar,h,x)-psi(y_i,h,x) */
  VFLOAT loss;       /**< loss(y_i,ybar) */
  VFLOAT alpha;      /**< dual parameter alpha_i^{ybar} for this sample */
  VFLOAT sqr;        /**< ||dpsi||^2 */
  VFLOAT slack;      /**< slack <w_t, psi(ybar,h,x)-psi(y_i,h,x)> + loss */
  VFLOAT slack_orig; /**< slack <w_i, psi(ybar,h,x)-psi(y_i,h,x)> + loss */
} SVM_cached_sample;

/**
 * @struct _SVM_cached_sample
 *
 * @brief Helper struct used by StructuredSVMOnlineLearner.  Encodes data associated with a call to find_most_violated_constraint(). When
 * sparm->method==SPO_DUAL_MULTI_SAMPLE_..., this may include a set of L labels that we want to jointly optimize over.  Otherwise,
 * it will just contain a single label
 * when
 */
typedef struct _SVM_cached_sample_set {
  SVM_cached_sample *samples;  /**< ybar_i^1...ybar_i^1^L */
  int num_samples;             /**< Number of samples L */
  int i;                       /**< training example index (index into sample.examples) */
} SVM_cached_sample_set;


void free_SVM_cached_sample(SVM_cached_sample *s);
void read_SVM_cached_sample(SVM_cached_sample *s, FILE *fin, StructuredSVM *svm);
void write_SVM_cached_sample(SVM_cached_sample *s, FILE *fout, StructuredSVM *svm);
SVM_cached_sample_set *new_SVM_cached_sample_set(int i);
void free_SVM_cached_sample_set(SVM_cached_sample_set *s);
SVM_cached_sample_set *read_SVM_cached_sample_set(FILE *fin, StructuredSVM *svm);
void write_SVM_cached_sample_set(SVM_cached_sample_set *s, FILE *fout, StructuredSVM *svm);
void write_SVM_cached_sample_set(SVM_cached_sample_set *s, FILE *fout, StructuredSVM *svm);
SVM_cached_sample *SVM_cached_sample_set_add_sample(SVM_cached_sample_set *s, StructuredLabel *ybar, SparseVector *dpsi, VFLOAT l, bool *regularize);
void SVM_cached_sample_set_update_parameters(SVM_cached_sample *c, SparseVector *sum_w, VFLOAT s, VFLOAT d, VFLOAT lambda, long t, VFLOAT featureScale,
                                             SVM_cached_sample_set *set, SparseVector *w_i, VFLOAT *sum_alpha, VFLOAT *L_i);
void SVM_cached_sample_optimize_dual(SVM_cached_sample *s, SparseVector *sum_w, VFLOAT lambda, long t, VFLOAT featureScale,
                                     SVM_cached_sample_set *set=NULL, SparseVector *w_i=NULL, VFLOAT *sum_alpha=NULL, VFLOAT *L_i=NULL);
void SVM_cached_sample_set_optimize_dual(SVM_cached_sample_set *s, SparseVector *sum_w, VFLOAT lambda, long t, VFLOAT featureScale, int R, bool *regularize);



#endif
