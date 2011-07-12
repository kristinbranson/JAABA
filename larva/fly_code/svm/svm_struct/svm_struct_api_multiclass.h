#ifndef  __SVM_STRUCT_API_MULTICLASS__
#define __SVM_STRUCT_API_MULTICLASS__

#include "svm_struct_api_types.h"
#include "svm_struct_api.h"

// When SVMMulticlass is used, LABEL->data can be cast to a (CLASS_LABEL*)
// and SPATTERN->data can be cast to a DOC*

typedef struct _class_label {
  /* this defines the y-part (the label) of a training example,
     e.g. the parse tree of the corresponding sentence. */
  int class_label;       /* class label */
  int num_classes; /* total number of classes */
  double *scores;  /* value of linear function of each class */
} CLASS_LABEL;


class SVMMulticlass : public SVMStructMethod {
  int num_classes;
  int num_features;

 public:
  void        svm_struct_learn_api_init(int argc, const char* argv[]);
  void        svm_struct_learn_api_exit();
  void        svm_struct_classify_api_init(int argc, const char* argv[]);
  void        svm_struct_classify_api_exit();
  SAMPLE      read_struct_examples(char *file, STRUCT_LEARN_PARM *sparm);
  void        init_struct_model(SAMPLE sample, STRUCTMODEL *sm, 
			      STRUCT_LEARN_PARM *sparm, LEARN_PARM *lparm, 
			      KERNEL_PARM *kparm);
  CONSTSET    init_struct_constraints(SAMPLE sample, STRUCTMODEL *sm, 
				    STRUCT_LEARN_PARM *sparm);
  LABEL       find_most_violated_constraint_slackrescaling(SPATTERN x, LABEL y, 
						     STRUCTMODEL *sm, 
						     STRUCT_LEARN_PARM *sparm);
  LABEL       find_most_violated_constraint_marginrescaling(SPATTERN x, LABEL y, 
						     STRUCTMODEL *sm, 
						     STRUCT_LEARN_PARM *sparm);
  LABEL       classify_struct_example(SPATTERN x, STRUCTMODEL *sm, 
				    STRUCT_LEARN_PARM *sparm);
  int         empty_label(LABEL y);
  SVECTOR     *psi(SPATTERN x, LABEL y, STRUCTMODEL *sm, 
			  STRUCT_LEARN_PARM *sparm);
  double      loss(LABEL y, LABEL ybar, STRUCT_LEARN_PARM *sparm);
  int         finalize_iteration(double ceps, int cached_constraint,
			       SAMPLE sample, STRUCTMODEL *sm,
			       CONSTSET cset, double *alpha, 
			       STRUCT_LEARN_PARM *sparm);
  void        print_struct_learning_stats(SAMPLE sample, STRUCTMODEL *sm,
					CONSTSET cset, double *alpha, 
					STRUCT_LEARN_PARM *sparm);
  void        print_struct_testing_stats(SAMPLE sample, STRUCTMODEL *sm,
				       STRUCT_LEARN_PARM *sparm,
				       STRUCT_TEST_STATS *teststats);
  void        eval_prediction(long exnum, EXAMPLE ex, LABEL prediction, 
			    STRUCTMODEL *sm, STRUCT_LEARN_PARM *sparm,
			    STRUCT_TEST_STATS *teststats);
  void        write_struct_model(const char *file,STRUCTMODEL *sm, 
			       STRUCT_LEARN_PARM *sparm);
  STRUCTMODEL read_struct_model(const char *file, STRUCT_LEARN_PARM *sparm);
  void        write_label(FILE *fp, LABEL y);
  void        free_pattern(SPATTERN x);
  void        free_label(LABEL y);
  void        free_struct_model(STRUCTMODEL sm);
  void        free_struct_sample(SAMPLE s);
  void        print_struct_help();
  void        parse_struct_parameters(STRUCT_LEARN_PARM *sparm);
  void        print_struct_help_classify();
  void        parse_struct_parameters_classify(STRUCT_LEARN_PARM *sparm);


  SAMPLE struct_examples_from_array(float *feat, int *classes, int num_examples, int num_features, int num_classes, unsigned char *mask);
  void svm_struct_classify_from_array(float *feat, int *preds, float *scores, int num_examples, STRUCTMODEL *model, STRUCT_LEARN_PARM *sparm) ;
  void svm_struct_learn_from_array(STRUCT_LEARN_PARM *struct_parm, STRUCTMODEL *structmodel, 
				   float *feat, int *classes, unsigned char *samples, int num_examples, int num_features, int num_classes, int argc, const char *argv[]);
  void struct_learn_read_input_parameters(int argc,const char *argv[],
				      long *verbosity,long *struct_verbosity, 
				      STRUCT_LEARN_PARM *struct_parm,
				      LEARN_PARM *learn_parm, KERNEL_PARM *kernel_parm,
				      int *alg_type);
  void struct_learn_print_help();
  void struct_learn_wait_any_key();
};


#endif
