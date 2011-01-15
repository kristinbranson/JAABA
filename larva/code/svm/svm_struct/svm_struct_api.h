/***********************************************************************/
/*                                                                     */
/*   svm_struct_api.h                                                  */
/*                                                                     */
/*   Definition of API for attaching implementing SVM learning of      */
/*   structures (e.g. parsing, multi-label classification, HMM)        */ 
/*                                                                     */
/*   Author: Thorsten Joachims                                         */
/*   Date: 03.07.04                                                    */
/*                                                                     */
/*   Copyright (c) 2004  Thorsten Joachims - All rights reserved       */
/*                                                                     */
/*   This software is available for non-commercial use only. It must   */
/*   not be modified and distributed without prior permission of the   */
/*   author. The author is not responsible for implications from the   */
/*   use of this software.                                             */
/*                                                                     */
/***********************************************************************/

#ifndef svm_struct_api
#define svm_struct_api

#include "svm_struct_api_types.h"
#include "svm_struct/svm_struct_common.h"

class SVMStructMethod {
 protected:
  char load_const_set_fname[400], save_const_set_fname[400], save_constraints_fname[400], load_constraints_fname[400];

 public:

  const char *save_constraints_file() { return save_constraints_fname; }
  void set_constraint_files(char *l, char *s) { strcpy(load_constraints_fname, l); strcpy(save_constraints_fname, s); }

  virtual void        svm_struct_learn_api_init(int argc, const char* argv[]) {};
  virtual void        svm_struct_learn_api_exit() {};
  virtual void        svm_struct_classify_api_init(int argc, const char* argv[]) {};
  virtual void        svm_struct_classify_api_exit() {};
  virtual SAMPLE      read_struct_examples(char *file, STRUCT_LEARN_PARM *sparm) { SAMPLE s; memset(&s,0,sizeof(s)); return s;};
  virtual void        init_struct_model(SAMPLE sample, STRUCTMODEL *sm, 
			      STRUCT_LEARN_PARM *sparm, LEARN_PARM *lparm, 
			      KERNEL_PARM *kparm) {};
  virtual CONSTSET    init_struct_constraints(SAMPLE sample, STRUCTMODEL *sm, 
					      STRUCT_LEARN_PARM *sparm) {CONSTSET c; memset(&c,0,sizeof(c)); return c;};
  virtual LABEL       find_most_violated_constraint_slackrescaling(SPATTERN x, LABEL y, 
						     STRUCTMODEL *sm, 
						     STRUCT_LEARN_PARM *sparm) { LABEL yy; memset(&yy,0,sizeof(yy)); return yy; };
  virtual LABEL       find_most_violated_constraint_marginrescaling(SPATTERN x, LABEL y, 
						     STRUCTMODEL *sm, 
						     STRUCT_LEARN_PARM *sparm) { LABEL yy; memset(&yy,0,sizeof(yy)); return yy; };
  virtual LABEL       classify_struct_example(SPATTERN x, STRUCTMODEL *sm, 
					      STRUCT_LEARN_PARM *sparm) { LABEL y; memset(&y,0,sizeof(y)); return y; };
  virtual int         empty_label(LABEL y) { return 0; };
  virtual SVECTOR     *psi(SPATTERN x, LABEL y, STRUCTMODEL *sm, 
			   STRUCT_LEARN_PARM *sparm) { return NULL; };
  virtual double      loss(LABEL y, LABEL ybar, STRUCT_LEARN_PARM *sparm) { return 0; };
  virtual int         finalize_iteration(double ceps, int cached_constraint,
			       SAMPLE sample, STRUCTMODEL *sm,
			       CONSTSET cset, double *alpha, 
					 STRUCT_LEARN_PARM *sparm) {return 0;};
  virtual void        print_struct_learning_stats(SAMPLE sample, STRUCTMODEL *sm,
					CONSTSET cset, double *alpha, 
					STRUCT_LEARN_PARM *sparm) {};
  virtual void        print_struct_testing_stats(SAMPLE sample, STRUCTMODEL *sm,
				       STRUCT_LEARN_PARM *sparm,
				       STRUCT_TEST_STATS *teststats) {};
  virtual void        eval_prediction(long exnum, EXAMPLE ex, LABEL prediction, 
			    STRUCTMODEL *sm, STRUCT_LEARN_PARM *sparm,
			    STRUCT_TEST_STATS *teststats) {};
  virtual void        write_struct_model(const char *file,STRUCTMODEL *sm, 
			       STRUCT_LEARN_PARM *sparm) {};
  virtual STRUCTMODEL read_struct_model(const char *file, STRUCT_LEARN_PARM *sparm) {STRUCTMODEL m; memset(&m,0,sizeof(m)); return m; };
  virtual void        write_label(FILE *fp, LABEL y) {};
  virtual void        free_pattern(SPATTERN x) {};
  virtual void        free_label(LABEL y) {};
  virtual void        free_struct_model(STRUCTMODEL sm) {};
  virtual void        free_struct_sample(SAMPLE s) {};
  virtual void        print_struct_help() {};
  virtual void        parse_struct_parameters(STRUCT_LEARN_PARM *sparm) {};
  virtual void        print_struct_help_classify() {};
  virtual void        parse_struct_parameters_classify(STRUCT_LEARN_PARM *sparm) {};
  virtual void on_finished_iteration(CONSTSET c, STRUCTMODEL *sm, 
						STRUCT_LEARN_PARM *sparm, int iter_num) {};
  virtual void write_label(LABEL y, FILE *fout, int iter, int num) {};
  virtual LABEL read_label(FILE *fin, char *fname, int *iter, int *num) {LABEL y; memset(&y,0,sizeof(y)); return y; };
  virtual void save_example(void *b, const char *fname) {}
  virtual char **load_examples(const char *fname, int *num) { return NULL; }
};



SVMStructMethod   *read_input_parameters_train(int, const char **, char *, char *,long *, long *,
					       STRUCT_LEARN_PARM *, LEARN_PARM *, KERNEL_PARM *,
					       int *, SVMStructMethod *m = NULL);


#endif
