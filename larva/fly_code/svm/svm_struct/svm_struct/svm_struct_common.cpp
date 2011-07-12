/***********************************************************************/
/*                                                                     */
/*   svm_struct_common.h                                               */
/*                                                                     */
/*   Functions and types used by multiple components of SVM-struct.    */
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

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "svm_struct_common.h"
#include "../svm_struct_api.h"
#include "svm_struct_learn.h"
#include "../../../blob.h"
#include "../../../train.h"
#include "../../../fit_model.h"
#include "../../../spine_features.h"
#include "../svm_struct_api_blob_behavior_sequence.h"
#include "../svm_struct_api_fly_behavior_sequence.h"
#include "../svm_struct_api_behavior_sequence.h"
#include "../svm_struct_api_multiclass.h"
#include "../online_structured_learning.h"


#define DATA_DIR "../data"

#ifdef USE_OPENMP
#include <omp.h>
#endif

void wait_any_key();

long struct_verbosity;                   /* verbosity level (0-4) */

void printIntArray(int* x, int n)
{
  int i;
  for(i=0;i<n;i++)
    printf("%i:",x[i]);
}

void printDoubleArray(double* x, int n)
{
  int i;
  for(i=0;i<n;i++)
    printf("%f:",x[i]);
}

void printWordArray(SWORD* x)
{
  int i=0;
  for(;x[i].wnum!=0;i++)
    if(x[i].weight != 0)
      printf(" %i:%.2f ",(int)x[i].wnum,x[i].weight);
}

void printW(double *w, long sizePhi, long n,double C)
{
  int i;
  printf("---- w ----\n");
  for(i=0;i<sizePhi;i++)
    {
      printf("%f  ",w[i]);
    }
  printf("\n----- xi ----\n");
  for(;i<sizePhi+2*n;i++)
    {
      printf("%f ",1/sqrt(2*C)*w[i]);
    }
  printf("\n");

}

void save_const_set(CONSTSET cset, const char *fname) {
  FILE *fout = fopen(fname, "w");
  int i, j;
  SVECTOR *v;

  fprintf(fout, "%d  # num constraints\n", cset.m);
  for(i = 0; i < cset.m; i++) 
    fprintf(stderr, "%lf ", cset.rhs[i]);
  fprintf(fout, "  # rhs\n");
  for(i = 0; i < cset.m; i++) {
    for(v=cset.lhs[i]->fvec; v; v=v->next) {
      fprintf(fout,"%.32g ", v->factor);
      for (j=0; (v->words[j]).wnum; j++) {
	fprintf(fout,"%ld:%.8g ", (long)(v->words[j]).wnum, (double)(v->words[j]).weight);
      }
      if(v->userdefined)
	fprintf(fout,"#%s\n",v->userdefined);
      else
	fprintf(fout,"#\n");
    }
  }

  fclose(fout);
}

CONSTSET load_const_set(const char *fname) {
  CONSTSET cset;
  long max_sv, max_words, ll, i, queryid, slackid, wpos;
  double alpha, costfactor;
  char *line,*comment;
  SWORD *words;

  nol_ll(fname,&max_sv,&max_words,&ll); /* scan size of model file */
  max_words+=2;
  ll+=2;
  words = (SWORD *)my_malloc(sizeof(SWORD)*(max_words+10));
  line = (char *)my_malloc(sizeof(char)*ll);

  FILE *fin = fopen(fname, "r");
  fscanf(fin,"%d%*[^\n]\n", &cset.m);
  cset.lhs = (DOC **)my_malloc(sizeof(DOC *)*cset.m);
  cset.rhs = (double*)my_malloc(sizeof(double)*cset.m);

  for(i = 0; i < cset.m; i++) 
    fscanf(fin,"%lf ", &cset.rhs[i]);
  fscanf(fin,"%*[^\n]\n");
  
  for(i = 0; i < cset.m; i++)  {
    fgets(line, 500000, fin);
    if(!parse_document(line,words,&alpha,&queryid,&slackid,
		       &costfactor,&wpos,max_words,&comment)) {
      printf("\nParsing constraint set %s line %s", fname, line);
      exit(1);
    }
    cset.lhs[i] = create_example(-1, 0,0,0.0, create_svector(words,comment,1.0));
  }
  fclose(fin);

  return cset;
}

/**** end print methods ****/


void train_main(int argc, const char** argv, STRUCT_LEARN_PARM *struct_parm, STRUCTMODEL *structmodel, SVMStructMethod *m, StructuredSVMOnlineLearner **learner_ptr) {  
  SAMPLE sample;  /* training sample */
  LEARN_PARM learn_parm;
  KERNEL_PARM kernel_parm;
  int alg_type;
  char modelfile[1000], modelfile_in[1000], trainfile[1000];
  
  strcpy(modelfile, "");

  m = read_input_parameters_train(argc,argv,trainfile,modelfile,&verbosity,
				  &struct_verbosity,struct_parm,&learn_parm,
				  &kernel_parm,&alg_type, m, modelfile_in);
  m->svm_struct_learn_api_init(argc,argv);

  if(strlen(modelfile)) {
	  // Save command line arguments
	  char fname[1000];
	  sprintf(fname, "%s.sh", modelfile);
	  FILE *fout = fopen(fname, "w");
	  assert(fout);
	  for(int i = 0; i < argc; i++) {
		  if(i) fprintf(fout, " ");
		  fprintf(fout, "%s", argv[i]);
	  }
	  fclose(fout);
  }

  char learner[1000];  
  sprintf(learner, "%s.learner", strlen(modelfile_in) ? modelfile_in : modelfile);
  if(strlen(modelfile_in) && FileExists(modelfile_in)) {
    StructuredSVMOnlineLearner *learner = new StructuredSVMOnlineLearner(m, modelfile_in);
	if(*learner_ptr) *learner_ptr = learner;
    learner->Train();
    *structmodel = *learner->GetStructModel();
    delete learner;
	return;
  } else {
  if(struct_verbosity>=1) {
    printf("Reading training examples..."); fflush(stdout);
  }
  /* read the training examples */
  sample = m->read_struct_examples(trainfile,struct_parm);
  if(struct_verbosity>=1) {
    printf("done\n"); fflush(stdout);
  }
  
  /* Do the learning and return structmodel. */
  if(struct_parm->method != SPO_CUTTING_PLANE) {
    StructuredSVMOnlineLearner *learner = new StructuredSVMOnlineLearner(m, &learn_parm, &kernel_parm, struct_parm, &sample, trainfile, modelfile);
	if(*learner_ptr) *learner_ptr = learner;
    learner->Train();
    *structmodel = *learner->GetStructModel();
    delete learner;
  } else {
  if(alg_type == 0)
    svm_learn_struct(sample,struct_parm,&learn_parm,&kernel_parm,structmodel,NSLACK_ALG,m);
  else if(alg_type == 1)
    svm_learn_struct(sample,struct_parm,&learn_parm,&kernel_parm,structmodel,NSLACK_SHRINK_ALG,m);
  else if(alg_type == 2)
    svm_learn_struct_joint(sample,struct_parm,&learn_parm,&kernel_parm,structmodel,ONESLACK_PRIMAL_ALG,m);
  else if(alg_type == 3)
    svm_learn_struct_joint(sample,struct_parm,&learn_parm,&kernel_parm,structmodel,ONESLACK_DUAL_ALG,m);
  else if(alg_type == 4)
    svm_learn_struct_joint(sample,struct_parm,&learn_parm,&kernel_parm,structmodel,ONESLACK_DUAL_CACHE_ALG,m);
  else if(alg_type == 9)
    svm_learn_struct_joint_custom(sample,struct_parm,&learn_parm,&kernel_parm,structmodel,m);
  else
    exit(1);
  }

  /* Warning: The model contains references to the original data 'docs'.
     If you want to free the original data, and only keep the model, you 
     have to make a deep copy of 'model'. */
  if(struct_verbosity>=1) {
    printf("Writing learned model...");fflush(stdout);
  }
  if(strlen(modelfile)) 
    m->write_struct_model(modelfile,structmodel,struct_parm);
  if(struct_verbosity>=1) {
    printf("done\n");fflush(stdout);
  }

  m->free_struct_sample(sample);

  m->svm_struct_learn_api_exit();
  }
}
/*
SVMStructMethod *instantiate_SVMStructMethod(char *classtype)
{
#define SVM_FLY_BEHAVIOR_SEQUENCE_STRING "robiea"
#define SVM_BLOB_BEHAVIOR_SEQUENCE_STRING "zlatic"

  if(!classtype) {
	printf("SVMMultiClass\n");
	return new SVMMultiClass();
  }
  else if(strstr(classtype, SVM_FLY_BEHAVIOR_SEQUENCE_STRING)) {
  	printf("SVMFlyBehaviorSequence\n");
        return new SVMFlyBehaviorSequence(feat_name, behaviors, -1);
  }
  else if(strstr(classtype, SVM_BLOB_BEHAVIOR_SEQUENCE_STRING)) {
        printf("SVMBlobBehaviorSequence\n");
//	return new SVMBlobBehaviorSequence( FitParams *p    , behaviors, -1,    SVMFeatureParams *sparams);
  }
  else {
        printf("Unsupported classType\n");
	return 0;
  }
}
*/

void print_help_train(SVMStructMethod *m);

SVMStructMethod *read_input_parameters_train(int argc,const char *argv[],char *trainfile,
			   char *modelfile,
			   long *verbosity,long *struct_verbosity, 
			   STRUCT_LEARN_PARM *struct_parm,
			   LEARN_PARM *learn_parm, KERNEL_PARM *kernel_parm,
				       int *alg_type, SVMStructMethod *m, char *modelfile_in)
{
  long i;
  char type[100], cname[400], bname[400], feat_name[400], load_constraints_fname[400], save_constraints_fname[400];
  FitParams params;
  BehaviorGroups *behaviors;
  int beh;

  /* set default */
  (*alg_type)=DEFAULT_ALG_TYPE;
  struct_parm->C=-0.01;
  struct_parm->slack_norm=1;
  struct_parm->epsilon=DEFAULT_EPS;
  struct_parm->custom_argc=0;
  struct_parm->loss_function=DEFAULT_LOSS_FCT;
  struct_parm->loss_type=DEFAULT_RESCALING;
  struct_parm->newconstretrain=100;
  struct_parm->ccache_size=5;
  struct_parm->batch_size=100;

  strcpy (modelfile, "svm_struct_model");
  if(modelfile_in) strcpy (modelfile_in, "");
  strcpy (learn_parm->predfile, "trans_predictions");
  strcpy (learn_parm->alphafile, "");
  (*verbosity)=0;/*verbosity for svm_light*/
  (*struct_verbosity)=1; /*verbosity for struct learning portion*/
  learn_parm->biased_hyperplane=1;
  learn_parm->remove_inconsistent=0;
  learn_parm->skip_final_opt_check=0;
  learn_parm->svm_maxqpsize=10;
  learn_parm->svm_newvarsinqp=0;
  learn_parm->svm_iter_to_shrink=-9999;
  learn_parm->maxiter=100000;
  learn_parm->kernel_cache_size=40;
  learn_parm->svm_c=99999999;  /* overridden by struct_parm->C */
  learn_parm->eps=0.001;       /* overridden by struct_parm->epsilon */
  learn_parm->transduction_posratio=-1.0;
  learn_parm->svm_costratio=1.0;
  learn_parm->svm_costratio_unlab=1.0;
  learn_parm->svm_unlabbound=1E-5;
  learn_parm->epsilon_crit=0.001;
  learn_parm->epsilon_a=1E-10;  /* changed from 1e-15 */
  learn_parm->compute_loo=0;
  learn_parm->rho=1.0;
  learn_parm->xa_depth=0;
  kernel_parm->kernel_type=0;
  kernel_parm->poly_degree=3;
  kernel_parm->rbf_gamma=1.0;
  kernel_parm->coef_lin=1;
  kernel_parm->coef_const=1;
  strcpy(kernel_parm->custom,"empty");
  strcpy(type,"c");
  strcpy(save_constraints_fname, "");
  strcpy(load_constraints_fname, "");

  sprintf(bname, "%s/behaviors.txt", DATA_DIR);
  sprintf(cname, "%s/classifier", DATA_DIR);
  strcpy(feat_name, "");

  for(i=1;(i<argc) && ((argv[i])[0] == '-');i++) {
    switch ((argv[i])[1]) 
      { 
      case '?': print_help_train(m); exit(0);
      case 'a': i++; strcpy(learn_parm->alphafile,argv[i]); break;
      case 'c': i++; struct_parm->C=atof(argv[i]); break;
      case 'p': i++; struct_parm->slack_norm=atol(argv[i]); break;
      case 'e': i++; struct_parm->epsilon=atof(argv[i]); break;
      case 'k': i++; struct_parm->newconstretrain=atol(argv[i]); break;
      case 'h': i++; learn_parm->svm_iter_to_shrink=atol(argv[i]); break;
      case '#': i++; learn_parm->maxiter=atol(argv[i]); break;
      case 'm': i++; learn_parm->kernel_cache_size=atol(argv[i]); break;
      case 'w': i++; (*alg_type)=atol(argv[i]); break;
      case 'o': i++; struct_parm->loss_type=atol(argv[i]); break;
      case 'n': i++; learn_parm->svm_newvarsinqp=atol(argv[i]); break;
      case 'q': i++; learn_parm->svm_maxqpsize=atol(argv[i]); break;
      case 'x': 
	i++;
	params = default_parameters();
	behaviors = load_behaviors(bname, cname);
	assert(behaviors);
	beh = atol(argv[i]);
	m = new SVMBlobBehaviorSequence(&params, behaviors, beh);
	break; 
      case 'B': i++; sprintf(bname, "%s", argv[i]); behaviors = load_behaviors(bname, cname); break;
      case 'F': i++; sprintf(feat_name, "%s", argv[i]); break;
      case 'i': i++; strcpy(save_constraints_fname, argv[i]); break;
      case 'j': i++; strcpy(load_constraints_fname, argv[i]); break;
      case 'l': i++; struct_parm->loss_function=atol(argv[i]); break;
      case 'f': i++; struct_parm->ccache_size=atol(argv[i]); break;
      case 'b': i++; struct_parm->batch_size=atof(argv[i]); break;
      case 't': i++; kernel_parm->kernel_type=atol(argv[i]); break;
      case 'd': i++; kernel_parm->poly_degree=atol(argv[i]); break;
      case 'g': i++; kernel_parm->rbf_gamma=atof(argv[i]); break;
      case 's': i++; kernel_parm->coef_lin=atof(argv[i]); break;
      case 'r': i++; kernel_parm->coef_const=atof(argv[i]); break;
      case 'u': i++; strcpy(kernel_parm->custom,argv[i]); break;
      case '-': strcpy(struct_parm->custom_argv[struct_parm->custom_argc++],argv[i]);i++; strcpy(struct_parm->custom_argv[struct_parm->custom_argc++],argv[i]);break; 
      case 'v': i++; (*struct_verbosity)=atol(argv[i]); break;
      case 'y': i++; (*verbosity)=atol(argv[i]); break;
      case 'D':
          struct_parm->debug_weights = strstr(argv[i], "w") != NULL;
          struct_parm->debug_features = strstr(argv[i], "f") != NULL;
          struct_parm->debug_predictions = strstr(argv[i], "p") != NULL;
          struct_parm->debug_model = strstr(argv[i], "m") != NULL;
		  strcpy (struct_parm->debugdir, argv[++i]); 
		  break;
      default: printf("\nUnrecognized option %s!\n\n",argv[i]);
	       print_help_train(m);
	       exit(0);
      }
  }
  
  if(!m) {
    if(strlen(feat_name)) {
      m = new SVMFlyBehaviorSequence(feat_name, behaviors, -1);
    } else 
      m = new SVMMulticlass;
  }

  m->set_constraint_files(load_constraints_fname, save_constraints_fname);

  if(i>=argc) {
    printf("\nNot enough input parameters!\n\n");
    wait_any_key();
    print_help_train(m);
    exit(0);
  }
  strcpy (trainfile, argv[i]);
  if((i+1)<argc) {
    strcpy (modelfile, argv[i+1]);
  }
  if((i+2)<argc && modelfile_in) {
    strcpy (modelfile_in, argv[i+2]);
  }
  if(learn_parm->svm_iter_to_shrink == -9999) {
    learn_parm->svm_iter_to_shrink=100;
  }

  if((learn_parm->skip_final_opt_check) 
     && (kernel_parm->kernel_type == K_LINEAR)) {
    printf("\nIt does not make sense to skip the final optimality check for linear kernels.\n\n");
    learn_parm->skip_final_opt_check=0;
  }    
  if((learn_parm->skip_final_opt_check) 
     && (learn_parm->remove_inconsistent)) {
    printf("\nIt is necessary to do the final optimality check when removing inconsistent \nexamples.\n");
    wait_any_key();
    print_help_train(m);
    exit(0);
  }    
  if((learn_parm->svm_maxqpsize<2)) {
    printf("\nMaximum size of QP-subproblems not in valid range: %ld [2..]\n",learn_parm->svm_maxqpsize); 
    wait_any_key();
    print_help_train(m);
    exit(0);
  }
  if((learn_parm->svm_maxqpsize<learn_parm->svm_newvarsinqp)) {
    printf("\nMaximum size of QP-subproblems [%ld] must be larger than the number of\n",learn_parm->svm_maxqpsize); 
    printf("new variables [%ld] entering the working set in each iteration.\n",learn_parm->svm_newvarsinqp); 
    wait_any_key();
    print_help_train(m);
    exit(0);
  }
  if(learn_parm->svm_iter_to_shrink<1) {
    printf("\nMaximum number of iterations for shrinking not in valid range: %ld [1,..]\n",learn_parm->svm_iter_to_shrink);
    wait_any_key();
    print_help_train(m);
    exit(0);
  }
  if(struct_parm->C<0) {
    printf("\nYou have to specify a value for the parameter '-c' (C>0)!\n\n");
    wait_any_key();
    print_help_train(m);
    exit(0);
  }
  if(((*alg_type) < 0) || (((*alg_type) > 5) && ((*alg_type) != 9))) {
    printf("\nAlgorithm type must be either '0', '1', '2', '3', '4', or '9'!\n\n");
    wait_any_key();
    print_help_train(m);
    exit(0);
  }
  if(learn_parm->transduction_posratio>1) {
    printf("\nThe fraction of unlabeled examples to classify as positives must\n");
    printf("be less than 1.0 !!!\n\n");
    wait_any_key();
    print_help_train(m);
    exit(0);
  }
  if(learn_parm->svm_costratio<=0) {
    printf("\nThe COSTRATIO parameter must be greater than zero!\n\n");
    wait_any_key();
    print_help_train(m);
    exit(0);
  }
  if(struct_parm->epsilon<=0) {
    printf("\nThe epsilon parameter must be greater than zero!\n\n");
    wait_any_key();
    print_help_train(m);
    exit(0);
  }
  if((struct_parm->ccache_size<=0) && ((*alg_type) == 4)) {
    printf("\nThe cache size must be at least 1!\n\n");
    wait_any_key();
    print_help_train(m);
    exit(0);
  }
  if(((struct_parm->batch_size<=0) || (struct_parm->batch_size>100))  
     && ((*alg_type) == 4)) {
    printf("\nThe batch size must be in the interval ]0,100]!\n\n");
    wait_any_key();
    print_help_train(m);
    exit(0);
  }
  if((struct_parm->slack_norm<1) || (struct_parm->slack_norm>2)) {
    printf("\nThe norm of the slacks must be either 1 (L1-norm) or 2 (L2-norm)!\n\n");
    wait_any_key();
    print_help_train(m);
    exit(0);
  }
  if((struct_parm->loss_type != SLACK_RESCALING) 
     && (struct_parm->loss_type != MARGIN_RESCALING)) {
    printf("\nThe loss type must be either 1 (slack rescaling) or 2 (margin rescaling)!\n\n");
    wait_any_key();
    print_help_train(m);
    exit(0);
  }
  if(learn_parm->rho<0) {
    printf("\nThe parameter rho for xi/alpha-estimates and leave-one-out pruning must\n");
    printf("be greater than zero (typically 1.0 or 2.0, see T. Joachims, Estimating the\n");
    printf("Generalization Performance of an SVM Efficiently, ICML, 2000.)!\n\n");
    wait_any_key();
    print_help_train(m);
    exit(0);
  }
  if((learn_parm->xa_depth<0) || (learn_parm->xa_depth>100)) {
    printf("\nThe parameter depth for ext. xi/alpha-estimates must be in [0..100] (zero\n");
    printf("for switching to the conventional xa/estimates described in T. Joachims,\n");
    printf("Estimating the Generalization Performance of an SVM Efficiently, ICML, 2000.)\n");
    wait_any_key();
    print_help_train(m);
    exit(0);
  }

  m->parse_struct_parameters(struct_parm);

  return m;
}

void wait_any_key()
{
  printf("\n(more)\n");
  (void)getc(stdin);
}

void print_help_train(SVMStructMethod *m)
{
  printf("\nSVM-struct learning module: %s, %s, %s\n",INST_NAME,INST_VERSION,INST_VERSION_DATE);
  printf("   includes SVM-struct %s for learning complex outputs, %s\n",STRUCT_VERSION,STRUCT_VERSION_DATE);
  printf("   includes SVM-light %s quadratic optimizer, %s\n",VERSION,VERSION_DATE);
  copyright_notice();
  printf("   usage: svm_struct_learn [options] example_file model_file\n\n");
  printf("Arguments:\n");
  printf("         example_file-> file with training data\n");
  printf("         model_file  -> file to store learned decision rule in\n");

  printf("General Options:\n");
  printf("         -?          -> this help\n");
  printf("         -v [0..3]   -> verbosity level (default 1)\n");
  printf("         -y [0..3]   -> verbosity level for svm_light (default 0)\n");
  printf("Learning Options:\n");
  printf("         -c float    -> C: trade-off between training error\n");
  printf("                        and margin (default 0.01)\n");
  printf("         -x [1,2,3]  -> behavior: Train a behavior bout sequence segmenter for behavior group <behavior>\n");
  printf("         -p [1,2]    -> L-norm to use for slack variables. Use 1 for L1-norm,\n");
  printf("                        use 2 for squared slacks. (default 1)\n");
  printf("         -o [1,2]    -> Rescaling method to use for loss.\n");
  printf("                        1: slack rescaling\n");
  printf("                        2: margin rescaling\n");
  printf("                        (default %d)\n",DEFAULT_RESCALING);
  printf("         -l [0..]    -> Loss function to use.\n");
  printf("                        0: zero/one loss\n");
  printf("                        ?: see below in application specific options\n");
  printf("                        (default %d)\n",DEFAULT_LOSS_FCT);
  printf("Optimization Options (see [2][5]):\n");
  printf("         -w [0,..,9] -> choice of structural learning algorithm (default %d):\n",(int)DEFAULT_ALG_TYPE);
  printf("                        0: n-slack algorithm described in [2]\n");
  printf("                        1: n-slack algorithm with shrinking heuristic\n");
  printf("                        2: 1-slack algorithm (primal) described in [5]\n");
  printf("                        3: 1-slack algorithm (dual) described in [5]\n");
  printf("                        4: 1-slack algorithm (dual) with constraint cache [5]\n");
  printf("                        9: custom algorithm in svm_struct_learn_custom.c\n");
  printf("         -e float    -> epsilon: allow that tolerance for termination\n");
  printf("                        criterion (default %f)\n",DEFAULT_EPS);
  printf("         -k [1..]    -> number of new constraints to accumulate before\n"); 
  printf("                        recomputing the QP solution (default 100) (-w 0 and 1 only)\n");
  printf("         -f [5..]    -> number of constraints to cache for each example\n");
  printf("                        (default 5) (used with -w 4)\n");
  printf("         -b [1..100] -> percentage of training set for which to refresh cache\n");
  printf("                        when no epsilon violated constraint can be constructed\n");
  printf("                        from current cache (default 100%%) (used with -w 4)\n");
  printf("SVM-light Options for Solving QP Subproblems (see [3]):\n");
  printf("         -n [2..q]   -> number of new variables entering the working set\n");
  printf("                        in each svm-light iteration (default n = q). \n");
  printf("                        Set n < q to prevent zig-zagging.\n");
  printf("         -m [5..]    -> size of svm-light cache for kernel evaluations in MB\n");
  printf("                        (default 40) (used only for -w 1 with kernels)\n");
  printf("         -h [5..]    -> number of svm-light iterations a variable needs to be\n"); 
  printf("                        optimal before considered for shrinking (default 100)\n");
  printf("         -# int      -> terminate svm-light QP subproblem optimization, if no\n");
  printf("                        progress after this number of iterations.\n");
  printf("                        (default 100000)\n");
  printf("Kernel Options:\n");
  printf("         -t int      -> type of kernel function:\n");
  printf("                        0: linear (default)\n");
  printf("                        1: polynomial (s a*b+c)^d\n");
  printf("                        2: radial basis function exp(-gamma ||a-b||^2)\n");
  printf("                        3: sigmoid tanh(s a*b + c)\n");
  printf("                        4: user defined kernel from kernel.h\n");
  printf("         -d int      -> parameter d in polynomial kernel\n");
  printf("         -g float    -> parameter gamma in rbf kernel\n");
  printf("         -s float    -> parameter s in sigmoid/poly kernel\n");
  printf("         -r float    -> parameter c in sigmoid/poly kernel\n");
  printf("         -u string   -> parameter of user defined kernel\n");
  printf("Output Options:\n");
  printf("         -a string   -> write all alphas to this file after learning\n");
  printf("                        (in the same order as in the training set)\n");
  printf("Application-Specific Options:\n");
  m->print_struct_help();
  wait_any_key();

  printf("\nMore details in:\n");
  printf("[1] T. Joachims, Learning to Align Sequences: A Maximum Margin Aproach.\n");
  printf("    Technical Report, September, 2003.\n");
  printf("[2] I. Tsochantaridis, T. Joachims, T. Hofmann, and Y. Altun, Large Margin\n");
  printf("    Methods for Structured and Interdependent Output Variables, Journal\n");
  printf("    of Machine Learning Research (JMLR), Vol. 6(Sep):1453-1484, 2005.\n");
  printf("[3] T. Joachims, Making Large-Scale SVM Learning Practical. Advances in\n");
  printf("    Kernel Methods - Support Vector Learning, B. Schölkopf and C. Burges and\n");
  printf("    A. Smola (ed.), MIT Press, 1999.\n");
  printf("[4] T. Joachims, Learning to Classify Text Using Support Vector\n");
  printf("    Machines: Methods, Theory, and Algorithms. Dissertation, Kluwer,\n");
  printf("    2002.\n");
  printf("[5] T. Joachims, T. Finley, Chun-Nam Yu, Cutting-Plane Training of Structural\n");
  printf("    SVMs, Machine Learning Journal, to appear.\n");
}

SVMStructMethod *read_input_parameters_test(int argc,const char *argv[],char *testfile,
					    char *modelfile,char *predictionsfile, char *outfilelist, char *outdir,
			   STRUCT_LEARN_PARM *struct_parm, long *verbosity,long *struct_verbosity, SVMStructMethod *m);
void print_help_test(SVMStructMethod *);

int test_main(int argc, const char* argv[], STRUCT_LEARN_PARM *struct_parm, STRUCTMODEL *structmodel, SVMStructMethod *m)
{
  long correct=0,incorrect=0,no_accuracy=0;
  double runtime=0;
  double avgloss=0,l;
  FILE *predfl;
  STRUCT_TEST_STATS teststats;
  SAMPLE testsample;
  char testfile[200];
  char modelfile[200];
  char predictionsfile[200];
  char outfilelist[200];
  char outdir[200];
  STRUCTMODEL structmodel2;

  m = read_input_parameters_test(argc,argv,testfile,modelfile,predictionsfile,outfilelist,outdir,struct_parm,
				 &verbosity,&struct_verbosity, m);
  m->svm_struct_classify_api_init(argc,argv);
  
  if(strlen(outdir)) {
	  strcpy(struct_parm->debugdir, outdir);
	  // Save command line arguments
	  CreateDirectoryIfNecessary(outdir);
	  char fname[1000]; 
	  sprintf(fname, "%s/run.sh", outdir);
	  FILE *fout = fopen(fname, "w");
	  assert(fout);
	  for(int i = 0; i < argc; i++) {
		  if(i) fprintf(fout, " ");
		  fprintf(fout, "%s", argv[i]);
	  }
	  fclose(fout);
      sprintf(fname, "%s/index.html", outdir);
	  fout = fopen(fname, "w");
	  if(fout) fclose(fout);
  }

  if(struct_verbosity>=1) {
    printf("Reading model..."); fflush(stdout);
  }
  
  if(!structmodel) {
    structmodel2=m->read_struct_model(modelfile,struct_parm);
    structmodel = &structmodel2;
  }

  if(struct_verbosity>=1) {
    fprintf(stdout, "done.\n");
  }

  if(structmodel->svm_model->kernel_parm.kernel_type == K_LINEAR) { /* linear kernel */
    /* compute weight vector */
    add_weight_vector_to_linear_model(structmodel->svm_model);
    structmodel->w=structmodel->svm_model->lin_weights;
  }
  
  if(struct_verbosity>=1) {
    printf("Reading test examples..."); fflush(stdout);
  }
  testsample=m->read_struct_examples(testfile,struct_parm);
  if(struct_verbosity>=1) {
    printf("done.\n"); fflush(stdout);
  }
 
  if(struct_verbosity>=1) {
    printf("Classifying test examples..."); fflush(stdout);
  }

  if ((predfl = fopen (predictionsfile, "w")) == NULL)
  { perror (predictionsfile); exit (1); }

  int num;
  char **test_list = strlen(outfilelist) ? m->load_examples(!strcmp(outfilelist, "*") ? testfile : outfilelist, &num) : NULL;
  if(test_list) assert(num == testsample.n);

#ifdef USE_OPENMP
omp_lock_t lock;
omp_init_lock (&lock);
#pragma omp parallel for
#endif
  for(long i=0;i<testsample.n;i++) {
    double t1=get_runtime();
	fprintf(stderr, "Classifying sequence %d... file %s\n", (int)i, testsample.examples[i].labelname);
#if DEBUG > 0
	g_currFile = testsample.examples[i].labelname;
//	for(i=0; i<MAX_FEATURES; i++)
		//if(i < m->
#endif
	double score;
    LABEL y=m->classify_struct_example(&testsample.examples[i].x,structmodel,struct_parm, &score);
    fprintf(stderr, "  sequence %d done\n", (int)i);
    
#ifdef USE_OPENMP
    omp_set_lock(&lock);
#endif
    runtime+=(get_runtime()-t1);

    m->write_label(predfl,y);
    if(!strcmp(outfilelist, "*")) {
      char tmp[1000]; strcpy(tmp, test_list[i]);
      strcat(tmp, ".label");
	  if(strlen(outdir)) {
		  char folder[1000], fname[1000];
		  ExtractFolderAndFileName(tmp, folder, fname);
		  sprintf(tmp, "%s/%s", outdir, fname);
	  }
      m->save_example(((BehaviorBoutFeatures*)testsample.examples[i].x.data)->data, (BehaviorBoutSequence*)y.data, tmp);
    } else if(test_list)
      m->save_example(((BehaviorBoutFeatures*)testsample.examples[i].x.data)->data, (BehaviorBoutSequence*)y.data, test_list[i]);

    l=m->loss(testsample.examples[i].y,y,struct_parm);
    avgloss+=l;
    if(l == 0) 
      correct++;
    else
      incorrect++;
    m->eval_prediction(i,testsample.examples[i],y,structmodel,struct_parm,&teststats);

	  char ename[1000];
	  sprintf(ename, "%d", (int)i);
	  m->on_finished_find_most_violated_constraint(&y, &testsample.examples[i].y, -1, struct_parm, ename);


    if(m->empty_label(testsample.examples[i].y)) 
      { no_accuracy=1; } /* test data is not labeled */
    if(struct_verbosity>=2) {
      if((i+1) % 100 == 0) {
	printf("%ld..",i+1); fflush(stdout);
      }
    }
    m->free_label(y);
#ifdef USE_OPENMP
    omp_unset_lock(&lock);
#endif	    
  }  

  avgloss/=testsample.n;
  fclose(predfl);

  if(struct_verbosity>=1) {
    printf("done\n");
    printf("Runtime (without IO) in cpu-seconds: %.2f\n",
	   (float)(runtime/100.0));    
  }
  if((!no_accuracy) && (struct_verbosity>=1)) {
    printf("Average loss on test set: %.4f\n",(float)avgloss);
    printf("Zero/one-error on test set: %.2f%% (%ld correct, %ld incorrect, %d total)\n",(float)100.0*incorrect/testsample.n,correct,incorrect,testsample.n);
  }
  m->print_struct_testing_stats(testsample,structmodel,struct_parm,&teststats);
  m->free_struct_sample(testsample);
  if(&structmodel2 == structmodel)
    m->free_struct_model(structmodel2);

  m->svm_struct_classify_api_exit();

  return(0);
}


SVMStructMethod *read_input_parameters_test(int argc,const char *argv[],char *testfile,
					    char *modelfile,char *predictionsfile, char *outfilelist, char *outdir,
			   STRUCT_LEARN_PARM *struct_parm,
			   long *verbosity,long *struct_verbosity, SVMStructMethod *m)
{
  long i;
  FitParams params;
  BehaviorGroups *behaviors;
  int beh;
  char bname[400], cname[400], feat_name[400];
  
  /* set default */
  strcpy (modelfile, "svm_model");
  strcpy (predictionsfile, "svm_predictions"); 
  strcpy (outfilelist, ""); 
  strcpy (outdir, ""); 
  (*verbosity)=0;/*verbosity for svm_light*/
  (*struct_verbosity)=1; /*verbosity for struct learning portion*/
  struct_parm->custom_argc=0;

  sprintf(bname, "%s/behaviors.txt", DATA_DIR);
  sprintf(cname, "%s/classifier", DATA_DIR);
  strcpy(feat_name, "");

  for(i=1;(i<argc) && ((argv[i])[0] == '-');i++) {
    switch ((argv[i])[1]) 
      { 
      case 'h': print_help_test(m); exit(0);
      case '?': print_help_test(m); exit(0);
      case '-': strcpy(struct_parm->custom_argv[struct_parm->custom_argc++],argv[i]);i++; strcpy(struct_parm->custom_argv[struct_parm->custom_argc++],argv[i]);break; 
      case 'v': i++; (*struct_verbosity)=atol(argv[i]); break;
      case 'y': i++; (*verbosity)=atol(argv[i]); break;
      case 'B': i++; sprintf(bname, "%s", argv[i]); behaviors = load_behaviors(bname, cname); break;
      case 'F': i++; sprintf(feat_name, "%s", argv[i]); break;
      case 'x': 
	i++;
	params = default_parameters();
	sprintf(bname, "%s/behaviors.txt", DATA_DIR);
	sprintf(cname, "%s/classifier", DATA_DIR);
	behaviors = load_behaviors(bname, cname);
	assert(behaviors);
	beh = atol(argv[i]);
	m = new SVMBlobBehaviorSequence(&params, behaviors, beh);
	break; 
      case 'o': strcpy (outfilelist, argv[++i]); break;
      case 'O': strcpy (outfilelist, "*"); break;
      case 'D':
          struct_parm->debug_weights = strstr(argv[i], "w") != NULL;
          struct_parm->debug_features = strstr(argv[i], "f") != NULL;
          struct_parm->debug_predictions = strstr(argv[i], "p") != NULL;
          struct_parm->debug_model = strstr(argv[i], "m") != NULL;
		  strcpy (outdir, argv[++i]); 
		  break;
      default: printf("\nUnrecognized option %s!\n\n",argv[i]);
	       print_help_test(m);
	       exit(0);
      }
  }
  if(!m) {
    if(strlen(feat_name)) 
      m = new SVMFlyBehaviorSequence(feat_name, behaviors, -1);
    else
      m = new SVMMulticlass;
  }

  if((i+1)>=argc) {
    printf("\nNot enough input parameters!\n\n");
    print_help_test(m);
    exit(0);
  }
  strcpy (testfile, argv[i]);
  strcpy (modelfile, argv[i+1]);
  if((i+2)<argc) {
    strcpy (predictionsfile, argv[i+2]);
  }

  m->parse_struct_parameters_classify(struct_parm);
  return m;
}

void print_help_test(SVMStructMethod *m)
{
  printf("\nSVM-struct classification module: %s, %s, %s\n",INST_NAME,INST_VERSION,INST_VERSION_DATE);
  printf("   includes SVM-struct %s for learning complex outputs, %s\n",STRUCT_VERSION,STRUCT_VERSION_DATE);
  printf("   includes SVM-light %s quadratic optimizer, %s\n",VERSION,VERSION_DATE);
  copyright_notice();
  printf("   usage: svm_struct_classify [options] example_file model_file output_file\n\n");
  printf("options: -h         -> this help\n");
  printf("         -v [0..3]  -> verbosity level (default 2)\n\n");

  m->print_struct_help_classify();
}
