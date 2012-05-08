#ifndef ONLINE_STRUCTURED_LEARNER_H
#define ONLINE_STRUCTURED_LEARNER_H

#include "svm_struct_api.h"

#include <omp.h>
#include <time.h>

#define SAVE_PROGRESS

struct _SVM_cached_sample;
struct _SVM_cached_sample_set;


/**
 * @class StructuredSVMOnlineLearner
 *
 * @brief A class used to train a structured prediction model in online fashion.  Supports
 * a few different online algorithms, and training examples can either come from a predefined
 * dataset or can be streamed in incrementally
 */
class StructuredSVMOnlineLearner {
  omp_lock_t my_lock;
  double lambda; /**< regularization constant */
  double eps;    /**< desired approximation factor for detecting convergence */
  int R;        /**<  number of iterations per multi-sample dual update */

  SVMStructMethod *m;
  LEARN_PARM *learn_parm;
  STRUCT_LEARN_PARM *sparm;
  KERNEL_PARM *kernel_parm;
  STRUCTMODEL sm;
  SAMPLE sample;
  char *modelfile;
  char *trainfile;

  SVECTOR *sum_w;  /**< the unnormalized model weights w^t = sum_w*(t*lambda) */ 
  long curr;  /**< the next example to process */ 
  long t;  /**< number of iterations run so far */ 
  long n;  /**< number of examples that have been processed in at least one iteration */ 
  long just_added;  /**< if >= 0, a new example was just added and we may want to focus on it */ 
  int M;   /**< number of online updates per example (determined automatically) */ 
  bool finished;  /**< If set to true, then the Train() threads will all finish up after the next iteration */ 
  int num_thr;    /**< number of threads */ 

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
  double regularization_error;      /**< Current regularization error  */ 
  double *generalization_errors_by_n, *optimization_errors_by_n, *model_errors_by_n, *regularization_errors_by_n, *losses_by_n;
  double *generalization_errors_by_t, *iter_errors_by_t, *model_errors_by_t, *constraint_errors_by_t, *regularization_errors_by_t, *losses_by_t, *elapsed_time_by_t;
  long *iter_examples;
  long base_time, start_time;
  
  int *examples_by_iteration_number; /**< An array of size M, where each entry stores an index to the start of a linked list (implemented in 
					examples_by_iteration_next_ind).  Each linked list stores the indices of all examples that have been iterated over
					i times */
  int *examples_by_iteration_next_ind;  /**< An array of size sample.n.  For each example, contains the index of the next example in a linked list of examples
				       that have been iterated over the same number of times as this example */ 
  int *examples_by_iteration_prev_ind;  /**< An array of size sample.n.  For each example, contains the index of the previous example in a linked list of examples
				       that have been iterated over the same number of times as this example */ 
  long minItersBeforeNewExample;
  int currMinIterByExample;
  
public:

  /**
   * @brief Constructor
   * @param learn_parm Currently only used to specify regularization and target approximation level: learn_parm.svm_c and learn_parm.eps
   * @param kernel_parm Not currently used
   * @param sparm Struct used to store custom data for a specific structured prediction algorithm
   * @param trainfile File containing a list of training examples
   * @param modelfile File where to output learned models
   */
  StructuredSVMOnlineLearner(SVMStructMethod *m, LEARN_PARM *learn_parm, KERNEL_PARM *kernel_parm, STRUCT_LEARN_PARM *sparm, SAMPLE *sample, const char *trainfile, const char *modelfile);

  /**
   * @brief Constructor
   * @param modelfile Read model from this file.  This allows us to resume from a previously time a StructuredSVMOnlineLearner was run and saved
   */
  StructuredSVMOnlineLearner(SVMStructMethod *m, const char *modelfile);

  /**
   * @brief Destructor
   * @param modelfile Read model from this file.  This allows us to resume from a previously time a StructuredSVMOnlineLearner was run and saved
   */
  ~StructuredSVMOnlineLearner();

  /**
   * @brief Main loop implementing online learning algorithm.  Continuously iterates over each training example
   */
  void Train();  

  /**
   * @brief Dynamically add a new example to be processed by the online learning algorithm
   * @param fname File name of the new training example, to be read by read_struct_example
   */
  void AddExample(const char *fname);

  /**
   * @brief Save entire state of this StructuredSVMOnlineLearner.  Learning should be resumable via a call to "new StructuredSVMOnlineLearner(fname)"
   * @param fname File name of where to store learned models
   */
  const char *SaveModel(const char *fname=NULL);

  /**
   * @brief Read entire state of this StructuredSVMOnlineLearner.  Typically this is invoked from the constructor and shouldn't be called directly
   * @param fname File name of where to read learned models.  This should be the file name passed to SaveModel(), or if the model was saved
   * automatically, it could be the modelfile parameter passed to the constructor
   */
  void ReadModel(const char *fname=NULL);

  long GetElapsedTime() { return time(NULL) - start_time + base_time; }

  
  void SetLambda(double l) { omp_set_lock(&my_lock); lambda = l;  sparm->C = (double)(1.0/lambda); omp_unset_lock(&my_lock); }
  void SetC(double c) { SetLambda(1.0/c); }

  
  SVMStructMethod *GetStructMethod() { return m; }
  STRUCT_LEARN_PARM *GetStructLearnParms() { return sparm; }
  STRUCTMODEL *GetStructModel() { return &sm; }
  SAMPLE *GetSample() { return &sample; }
  void Shutdown() { finished = true; }
  void GetStatisticsByIteration(int ave, long *tt, long *tm, double **gen_err_buff, double **opt_err_buff, double **model_err_buff, 
				double **reg_err_buff, double **train_err_buff, double **test_err_buff, double **time_buff);
  void GetStatisticsByExample(int ave, long *nn, double **gen_err_buff, double **opt_err_buff, double **model_err_buff, 
			      double **reg_err_buff, double **train_err_buff, double **test_err_buff);
  double *GetCurrentWeights();

private:
  void Init();

  /**
   * @brief Main function in the inner loop of online learning algorithm, used to update weights
   * @param ex The training example we are using to update the model weights.  
   * @param isNewIter If true, increment the iteration number t.  Otherwise, it is assumed we are just updating an old example
   * @param tt The iteration number for this update
   */
  void UpdateWeights(struct _SVM_cached_sample_set *ex, int isNewIter=-1);
  
  /**
   * @brief Sanity check: recompute weights from dual parameters.  Could help avoid drifting due to numerical precision errors
   */
  void RecomputeWeights(); 

  int ChooseNextExample();
};




/**
 * @struct _SVM_cached_sample
 * 
 * @brief Helper struct used by StructuredSVMOnlineLearner.  Encodes data associated with a call to find_most_violated_constraint()
 */
typedef struct _SVM_cached_sample {
  LABEL ybar;       /**< label */
  SVECTOR *dpsi;    /**< gradient w.r.t. ith example at ybar: psi(ybar,h,x)-psi(y_i,h,x) */
  double loss;       /**< loss(y_i,ybar) */
  double alpha;      /**< dual parameter alpha_i^{ybar} for this sample */
  double sqr;        /**< ||dpsi||^2 */
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

void free_SVM_cached_sample(SVM_cached_sample *s, SVMStructMethod *m);
void read_SVM_cached_sample(SVM_cached_sample *s, FILE *fin, STRUCT_LEARN_PARM *sparm);
void write_SVM_cached_sample(SVM_cached_sample *s, FILE *fout, STRUCT_LEARN_PARM *sparm);

SVM_cached_sample_set *new_SVM_cached_sample_set(int i);
void free_SVM_cached_sample_set(SVM_cached_sample_set *s, SVMStructMethod *m);
SVM_cached_sample_set *read_SVM_cached_sample_set(FILE *fin, STRUCT_LEARN_PARM *sparm);
void write_SVM_cached_sample_set(SVM_cached_sample_set *s, FILE *fout, STRUCT_LEARN_PARM *sparm);
void SVM_cached_sample_set_add_sample(SVM_cached_sample_set *s, LABEL ybar, SVECTOR *dpsi, double l);

// Update weights sum_w when the alpha parameters for sample set 'set' get scaled by 's', and the alpha parameter for sample 'c' gets 
// incremented by 'd'
SVECTOR *SVM_cached_sample_set_update_parameters(SVM_cached_sample *c, SVECTOR *sum_w, double s, double d, double lambda, long t, 
						 SVM_cached_sample_set *set=NULL, SVECTOR **w_i=NULL, double *sum_alpha=NULL, double *L_i=NULL);

// Optimize alpha parameter for a single sample 's'.  Implements 'Online Dual Update Step' of writeup
SVECTOR *SVM_cached_sample_optimize_dual(SVM_cached_sample *s, SVECTOR *sum_w, double lambda, long t, 
					 SVM_cached_sample_set *set=NULL, SVECTOR **w_i=NULL, double *sum_alpha=NULL, double *L_i=NULL);

// Iteratively optimize alpha parameters for a set of samples 's'.  Implements 'Multi-Sample Dual Update Step' of writeup
SVECTOR *SVM_cached_sample_set_optimize_dual(SVM_cached_sample_set *s, SVECTOR *sum_w, double lambda, long t, int R=1);

int svm_learn_online_main(SVMStructMethod *m, const char *trainfile, const char *modelfile, LEARN_PARM *learn_parm, KERNEL_PARM *kernel_parm, 
		      STRUCT_LEARN_PARM *sparm, const char *background_trainfile);

#endif
