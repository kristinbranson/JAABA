/***********************************************************************/
/*                                                                     */
/*   svm_struct_api_types.h                                            */
/*                                                                     */
/*   Definition of API for attaching implementing SVM learning of      */
/*   structures (e.g. parsing, multi-label classification, HMM)        */ 
/*                                                                     */
/*   Author: Thorsten Joachims                                         */
/*   Date: 13.10.03                                                    */
/*                                                                     */
/*   Copyright (c) 2003  Thorsten Joachims - All rights reserved       */
/*                                                                     */
/*   This software is available for non-commercial use only. It must   */
/*   not be modified and distributed without prior permission of the   */
/*   author. The author is not responsible for implications from the   */
/*   use of this software.                                             */
/*                                                                     */
/***********************************************************************/

#ifndef svm_struct_api_types
#define svm_struct_api_types

# include "svm_light/svm_common.h"
# include "svm_light/svm_learn.h"

# define INST_NAME          "Generic and empty API"
# define INST_VERSION       "V0.00"
# define INST_VERSION_DATE  "??.??.??"


/* default precision for solving the optimization problem */
# define DEFAULT_EPS         0.1 
/* default loss rescaling method: 1=slack_rescaling, 2=margin_rescaling */
# define DEFAULT_RESCALING   2
/* default loss function: */
# define DEFAULT_LOSS_FCT    0
/* default optimization algorithm to use: */
# define DEFAULT_ALG_TYPE    3
/* store Psi(x,y) (for ALG_TYPE 1) instead of recomputing it every time: */
# define USE_FYCACHE         1
/* decide whether to evaluate sum before storing vectors in constraint
   cache: 
   0 = NO, 
   1 = YES (best, if sparse vectors and long vector lists), 
   2 = YES (best, if short vector lists),
   3 = YES (best, if dense vectors and long vector lists) */
# define COMPACT_CACHED_VECTORS 1
/* minimum absolute value below which values in sparse vectors are
   rounded to zero. Values are stored in the FVAL type defined in svm_common.h 
   RECOMMENDATION: assuming you use FVAL=float, use 
     10E-15 if COMPACT_CACHED_VECTORS is 1 
     10E-10 if COMPACT_CACHED_VECTORS is 2 or 3 
*/
# define COMPACT_ROUNDING_THRESH 10E-15


typedef enum {
  SPO_CUTTING_PLANE,   // SVM^struct 

  // Online Algorithms, which process one training example i per iteration:
  SPO_SGD,             // stochastic gradient descent: w^t = w^{t-1} - step_size*(w^{t-1} + grad_i(w^{t-1})),   step_size=1/(lambda*t)
  SPO_SGD_PEGASOS,     // Do SGD step, then ensure that ||w^t||^2 <= 1/lambda   (downscale w^t if this doesn't hold)
  SPO_DUAL_UPDATE,     // w^t = w^{t-1}(t-1)/t - alpha*grad_i(w^{t-1}),   where alpha is chosen to maximize dual objective
  SPO_DUAL_UPDATE_WITH_CACHE,     // Same as SPO_DUAL_UPDATE, but also run dual update steps in a background process on labels from earlier iterations
  SPO_DUAL_MULTI_SAMPLE_UPDATE,   // Sample multiple labels per iteration (instead of just getting the "most violated constraint"), then optimize parameters jointly
  SPO_DUAL_MULTI_SAMPLE_UPDATE_WITH_CACHE,  // Same as SPO_DUAL_MULTI_SAMPLE_UPDATE, but also run multi-sample dual update steps in a background process 
} StructuredPredictionOptimizationMethod;


typedef struct pattern {
  /* this defines the x-part of a training example, e.g. the structure
     for storing a natural language sentence in NLP parsing */
  void *data;
} SPATTERN;

typedef struct label {
  /* this defines the y-part (the label) of a training example,
     e.g. the parse tree of the corresponding sentence. */
  void *data;
} LABEL;

typedef struct structmodel {
  double *w;          /* pointer to the learned weights */
  MODEL  *svm_model;  /* the learned SVM model */
  long   sizePsi;     /* maximum number of weights in w */
  double walpha;
  /* other information that is needed for the stuctural model can be
     added here, e.g. the grammar rules for NLP parsing */
  int add_your_variables_here;
  bool compactUnaryCosts; /* CSC: true iff case Unary costs are stored as same transition costs (and no further weights are given) */
  bool extraUnaryCosts; /* CSC: true iff case Unary costs are given in addition to same transition costs; equals !compactUnaryCosts but explicitely given for sanity checks (CHECK THE ASSIGNMENT WHENEVER ADDITIONAL WEIGHTS (other than feature & transition weights) ARE ADDED) */
} STRUCTMODEL;

typedef struct struct_learn_parm {
  StructuredPredictionOptimizationMethod method;

  double epsilon;              /* precision for which to solve
				  quadratic program */
  double newconstretrain;      /* number of new constraints to
				  accumulate before recomputing the QP
				  solution (used in w=1 algorithm) */
  int    ccache_size;          /* maximum number of constraints to
				  cache for each example (used in w=4
				  algorithm) */
  double batch_size;           /* size of the mini batches in percent
				  of training set size (used in w=4
				  algorithm) */
  double C;                    /* trade-off between margin and loss */
  char   custom_argv[50][300]; /* storage for the --* command line options */
  int    custom_argc;          /* number of --* command line options */
  int    slack_norm;           /* norm to use in objective function
                                  for slack variables; 1 -> L1-norm, 
				  2 -> L2-norm */
  int    loss_type;            /* selected loss type from -r
				  command line option. Select between
				  slack rescaling (1) and margin
				  rescaling (2) */
  int    loss_function;        /* select between different loss
				  functions via -l command line
				  option */

  /* further parameters that are passed to init_struct_model() */
  char debugdir[400];
  bool debug_predictions, debug_weights, debug_features, debug_model;
  int iter;

} STRUCT_LEARN_PARM;

typedef struct struct_test_stats {
  /* you can add variables for keeping statistics when evaluating the
     test predictions in svm_struct_classify. This can be used in the
     function eval_prediction and print_struct_testing_stats. */
} STRUCT_TEST_STATS;

#endif
