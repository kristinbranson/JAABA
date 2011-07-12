/***********************************************************************/
/*                                                                     */
/*   svm_struct_api.c                                                  */
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

#include <stdio.h>
#include <string.h>
#include "svm_struct/svm_struct_common.h"
#include "svm_struct_api_multiclass.h"
#include "svm_struct/svm_struct_learn.h"

void        SVMMulticlass::svm_struct_learn_api_init(int argc, const char* argv[])
{
  /* Called in learning part before anything else is done to allow
     any initializations that might be necessary. */
}

void        SVMMulticlass::svm_struct_learn_api_exit()
{
  /* Called in learning part at the very end to allow any clean-up
     that might be necessary. */
}

void        SVMMulticlass::svm_struct_classify_api_init(int argc, const char* argv[])
{
  /* Called in prediction part before anything else is done to allow
     any initializations that might be necessary. */
}

void        SVMMulticlass::svm_struct_classify_api_exit()
{
  /* Called in prediction part at the very end to allow any clean-up
     that might be necessary. */
}

SAMPLE      SVMMulticlass::read_struct_examples(char *file, STRUCT_LEARN_PARM *sparm)
{
  /* Reads training examples and returns them in sample. The number of
     examples must be written into sample.n */
  SAMPLE   sample;  /* sample */
  EXAMPLE  *examples;
  long     n;       /* number of examples */
  DOC **docs;       /* examples in original SVM-light format */ 
  double *target;
  long totwords,i,num_classes=0;
  CLASS_LABEL *y;

  /* Using the read_documents function from SVM-light */
  read_documents(file,&docs,&target,&totwords,&n);
  examples=(EXAMPLE *)my_malloc(sizeof(EXAMPLE)*n);
  for(i=0;i<n;i++)     /* find highest class label */
    if(num_classes < (target[i]+0.1)) 
      num_classes=target[i]+0.1;
  for(i=0;i<n;i++)     /* make sure all class labels are positive */
    if(target[i]<1) {
      printf("\nERROR: The class label '%lf' of example number %ld is not greater than '1'!\n",target[i],i+1);
      exit(1);
    } 
  for(i=0;i<n;i++) {          /* copy docs over into new datastructure */
    examples[i].x.data=docs[i];
    examples[i].y.data = y = (CLASS_LABEL*)malloc(sizeof(CLASS_LABEL));
    y->class_label=target[i]+0.1;
    y->scores=NULL;
    y->num_classes=num_classes;
  }
  free(target);
  free(docs);
  sample.n=n;
  sample.examples=examples;

  if(struct_verbosity>=0)
    printf(" (%d examples) ",sample.n);
  return(sample);
}


void        SVMMulticlass::init_struct_model(SAMPLE sample, STRUCTMODEL *sm, 
			      STRUCT_LEARN_PARM *sparm, LEARN_PARM *lparm, 
			      KERNEL_PARM *kparm)
{
  /* Initialize structmodel sm. The weight vector w does not need to be
     initialized, but you need to provide the maximum size of the
     feature space in sizePsi. This is the maximum number of different
     weights that can be learned. Later, the weight vector w will
     contain the learned weights for the model. */
  long i,totwords=0;
  SWORD *w;

  this->num_classes=1;
  for(i=0;i<sample.n;i++)     /* find highest class label */
    if(this->num_classes < (((CLASS_LABEL*)sample.examples[i].y.data)->class_label+0.1)) 
      this->num_classes=((CLASS_LABEL*)sample.examples[i].y.data)->class_label+0.1;
  for(i=0;i<sample.n;i++)     /* find highest feature number */
    for(w=((DOC*)sample.examples[i].x.data)->fvec->words;w->wnum;w++) 
      if(totwords < w->wnum) 
	totwords=w->wnum;
  this->num_features=totwords;
  if(struct_verbosity>=0)
    printf("Training set properties: %d features, %d classes\n",
	   this->num_features,this->num_classes);
  sm->sizePsi=this->num_features*this->num_classes;
  if(struct_verbosity>=2)
    printf("Size of Phi: %ld\n",sm->sizePsi);
}

CONSTSET    SVMMulticlass::init_struct_constraints(SAMPLE sample, STRUCTMODEL *sm, 
				    STRUCT_LEARN_PARM *sparm)
{
  /* Initializes the optimization problem. Typically, you do not need
     to change this function, since you want to start with an empty
     set of constraints. However, if for example you have constraints
     that certain weights need to be positive, you might put that in
     here. The constraints are represented as lhs[i]*w >= rhs[i]. lhs
     is an array of feature vectors, rhs is an array of doubles. m is
     the number of constraints. The function returns the initial
     set of constraints. */
  CONSTSET c;
  long     sizePsi=sm->sizePsi;
  long     i;
  SWORD     words[2];

  if(1) { /* normal case: start with empty set of constraints */
    c.lhs=NULL;
    c.rhs=NULL;
    c.m=0;
  }
  else { /* add constraints so that all learned weights are
            positive. WARNING: Currently, they are positive only up to
            precision epsilon set by -e. */
    c.lhs=(DOC**)my_malloc(sizeof(DOC *)*sizePsi);
    c.rhs=(double*)my_malloc(sizeof(double)*sizePsi);
    for(i=0; i<sizePsi; i++) {
      words[0].wnum=i+1;
      words[0].weight=1.0;
      words[1].wnum=0;
      /* the following slackid is a hack. we will run into problems,
         if we have move than 1000000 slack sets (ie examples) */
      c.lhs[i]=create_example(i,0,1000000+i,1,create_svector(words,NULL,1.0));
      c.rhs[i]=0.0;
    }
  }
  return(c);
}

LABEL       SVMMulticlass::classify_struct_example(SPATTERN x, STRUCTMODEL *sm, 
				    STRUCT_LEARN_PARM *sparm)
{
  /* Finds the label yhat for pattern x that scores the highest
     according to the linear evaluation function in sm, especially the
     weights sm.w. The returned label is taken as the prediction of sm
     for the pattern x. The weights correspond to the features defined
     by psi() and range from index 1 to index sm->sizePsi. If the
     function cannot find a label, it shall return an empty label as
     recognized by the function empty_label(y). */
  LABEL yy;
  CLASS_LABEL *y;
  DOC doc;
  long class_label,bestclass=-1,first=1,j;
  double score,bestscore=-1;
  SWORD *words;

  doc=*((DOC*)x.data);
  y = (CLASS_LABEL*)my_malloc(sizeof(CLASS_LABEL)+sizeof(double)*(this->num_classes+1));
  y->scores=(double *)(y+1);
  y->num_classes=this->num_classes;
  yy.data = y;

  words=doc.fvec->words;
  for(j=0;(words[j]).wnum != 0;j++) {       /* Check if feature numbers   */
    if((words[j]).wnum>this->num_features) /* are not larger than in     */
      (words[j]).wnum=0;                    /* model. Remove feature if   */
  }                                         /* necessary.                 */
  for(class_label=1;class_label<=this->num_classes;class_label++) {
    y->class_label=class_label;
    doc.fvec=psi(x,yy,sm,sparm);
    score=classify_example(sm->svm_model,&doc);
    free_svector(doc.fvec);
    y->scores[class_label]=score;
    if((bestscore<score)  || (first)) {
      bestscore=score;
      bestclass=class_label;
      first=0;
    }
  }
  y->class_label=bestclass;

  return(yy);
}

LABEL       SVMMulticlass::find_most_violated_constraint_slackrescaling(SPATTERN x, LABEL yy, 
						     STRUCTMODEL *sm, 
						     STRUCT_LEARN_PARM *sparm)
{
  /* Finds the label ybar for pattern x that that is responsible for
     the most violated constraint for the slack rescaling
     formulation. It has to take into account the scoring function in
     sm, especially the weights sm.w, as well as the loss
     function. The weights in sm.w correspond to the features defined
     by psi() and range from index 1 to index sm->sizePsi. Most simple
     is the case of the zero/one loss function. For the zero/one loss,
     this function should return the highest scoring label ybar, if
     ybar is unequal y; if it is equal to the correct label y, then
     the function shall return the second highest scoring label. If
     the function cannot find a label, it shall return an empty label
     as recognized by the function empty_label(y). */
  LABEL yybar;
  CLASS_LABEL *ybar = (CLASS_LABEL*)malloc(sizeof(CLASS_LABEL));
  DOC doc;
  long class_label,bestclass=-1,first=1;
  double score,score_y,score_ybar,bestscore=-1;

  /* NOTE: This function could be made much more efficient by not
     always computing a new PSI vector. */
  doc=*((DOC*)x.data);
  doc.fvec=psi(x,yy,sm,sparm);
  score_y=classify_example(sm->svm_model,&doc);
  free_svector(doc.fvec);

  yybar.data = ybar;
  ybar->scores=NULL;
  ybar->num_classes=this->num_classes;
  for(class_label=1;class_label<=this->num_classes;class_label++) {
    ybar->class_label=class_label;
    doc.fvec=psi(x,yybar,sm,sparm);
    score_ybar=classify_example(sm->svm_model,&doc);
    free_svector(doc.fvec);
    score=loss(yy,yybar,sparm)*(1.0-score_y+score_ybar);
    if((bestscore<score)  || (first)) {
      bestscore=score;
      bestclass=class_label;
      first=0;
    }
  }
  if(bestclass == -1) 
    printf("ERROR: Only one class\n");
  ybar->class_label=bestclass;
  if(struct_verbosity>=3)
    printf("[%ld:%.2f] ",bestclass,bestscore);
  return(yybar);
}

LABEL       SVMMulticlass::find_most_violated_constraint_marginrescaling(SPATTERN x, LABEL yy, 
						     STRUCTMODEL *sm, 
						     STRUCT_LEARN_PARM *sparm)
{
  /* Finds the label ybar for pattern x that that is responsible for
     the most violated constraint for the margin rescaling
     formulation. It has to take into account the scoring function in
     sm, especially the weights sm.w, as well as the loss
     function. The weights in sm.w correspond to the features defined
     by psi() and range from index 1 to index sm->sizePsi. Most simple
     is the case of the zero/one loss function. For the zero/one loss,
     this function should return the highest scoring label ybar, if
     ybar is unequal y; if it is equal to the correct label y, then
     the function shall return the second highest scoring label. If
     the function cannot find a label, it shall return an empty label
     as recognized by the function empty_label(y). */
  LABEL yybar;
  CLASS_LABEL *ybar = (CLASS_LABEL*)malloc(sizeof(CLASS_LABEL));
  DOC doc;
  long class_label,bestclass=-1,first=1;
  double score,bestscore=-1;

  /* NOTE: This function could be made much more efficient by not
     always computing a new PSI vector. */
  doc=*((DOC*)x.data);
  yybar.data = ybar;
  ybar->scores=NULL;
  ybar->num_classes=this->num_classes;
  for(class_label=1;class_label<=this->num_classes;class_label++) {
    ybar->class_label=class_label;
    doc.fvec=psi(x,yybar,sm,sparm);         
    score=classify_example(sm->svm_model,&doc);
    free_svector(doc.fvec);
    score+=loss(yy,yybar,sparm);
    if((bestscore<score)  || (first)) {
      bestscore=score;
      bestclass=class_label;
      first=0;
    }
  }
  if(bestclass == -1) 
    printf("ERROR: Only one class\n");
  ybar->class_label=bestclass;
  if(struct_verbosity>=3)
    printf("[%ld:%.2f] ",bestclass,bestscore);
  return(yybar);
}

int         SVMMulticlass::empty_label(LABEL y)
{
  /* Returns true, if y is an empty label. An empty label might be
     returned by find_most_violated_constraint_???(x, y, sm) if there
     is no incorrect label that can be found for x, or if it is unable
     to label x at all */
  return(((CLASS_LABEL*)y.data)->class_label<0.9);
}

SVECTOR     *SVMMulticlass::psi(SPATTERN x, LABEL yy, STRUCTMODEL *sm,
		 STRUCT_LEARN_PARM *sparm)
{
  /* Returns a feature vector describing the match between pattern x and
     label y. The feature vector is returned as an SVECTOR
     (i.e. pairs <featurenumber:featurevalue>), where the last pair has
     featurenumber 0 as a terminator. Featurenumbers start with 1 and end with
     sizePsi. This feature vector determines the linear evaluation
     function that is used to score labels. There will be one weight in
     sm.w for each feature. Note that psi has to match
     find_most_violated_constraint_???(x, y, sm) and vice versa. In
     particular, find_most_violated_constraint_???(x, y, sm) finds that
     ybar!=y that maximizes psi(x,ybar,sm)*sm.w (where * is the inner
     vector product) and the appropriate function of the loss.  */
  SVECTOR *fvec;
  CLASS_LABEL *y = (CLASS_LABEL*)yy.data;

  /* shift the feature numbers to the position of weight vector of class y */
  fvec=shift_s(((DOC*)x.data)->fvec,(y->class_label-1)*this->num_features);

  /* The following makes sure that the weight vectors for each class
     are treated separately when kernels are used . */
  fvec->kernel_id=y->class_label;

  return(fvec);
}

double      SVMMulticlass::loss(LABEL yy, LABEL yybar, STRUCT_LEARN_PARM *sparm)
{
  CLASS_LABEL *y = (CLASS_LABEL*)yy.data;
  CLASS_LABEL *ybar = (CLASS_LABEL*)yybar.data;

  /* loss for correct label y and predicted label ybar. The loss for
     y==ybar has to be zero. sparm->loss_function is set with the -l option. */
  if(sparm->loss_function == 0) { /* type 0 loss: 0/1 loss */
    if(y->class_label == ybar->class_label)     /* return 0, if y==ybar. return 100 else */
      return(0);
    else
      return(100);
  }
  if(sparm->loss_function == 1) { /* type 1 loss: squared difference */
    return((y->class_label-ybar->class_label)*(y->class_label-ybar->class_label));
  }
  else {
    /* Put your code for different loss functions here. But then
       find_most_violated_constraint_???(x, y, sm) has to return the
       highest scoring label with the largest loss. */
    printf("Unkown loss function\n");
    exit(1);
  }
}

int         SVMMulticlass::finalize_iteration(double ceps, int cached_constraint,
			       SAMPLE sample, STRUCTMODEL *sm,
			       CONSTSET cset, double *alpha, 
			       STRUCT_LEARN_PARM *sparm)
{
  /* This function is called just before the end of each cutting plane iteration. ceps is the amount by which the most violated constraint found in the current iteration was violated. cached_constraint is true if the added constraint was constructed from the cache. If the return value is FALSE, then the algorithm is allowed to terminate. If it is TRUE, the algorithm will keep iterating even if the desired precision sparm->epsilon is already reached. */
  return(0);
}

void        SVMMulticlass::print_struct_learning_stats(SAMPLE sample, STRUCTMODEL *sm,
					CONSTSET cset, double *alpha, 
					STRUCT_LEARN_PARM *sparm)
{
  /* This function is called after training and allows final touches to
     the model sm. But primarly it allows computing and printing any
     kind of statistic (e.g. training error) you might want. */

  /* Replace SV with single weight vector */
  MODEL *model=sm->svm_model;
  if(model->kernel_parm.kernel_type == K_LINEAR) {
    if(struct_verbosity>=1) {
      printf("Compacting linear model..."); fflush(stdout);
    }
    sm->svm_model=compact_linear_model(model);
    sm->w=sm->svm_model->lin_weights; /* short cut to weight vector */
    free_model(model,1);
    if(struct_verbosity>=1) {
      printf("done\n"); fflush(stdout);
    }
  }  
}

void        SVMMulticlass::write_struct_model(const char *file, STRUCTMODEL *sm, 
			       STRUCT_LEARN_PARM *sparm)
{
  /* Writes structural model sm to file file. */
  FILE *modelfl;
  long j,i,sv_num;
  MODEL *model=sm->svm_model;
  SVECTOR *v;

  if ((modelfl = fopen (file, "w")) == NULL)
  { perror (file); exit (1); }
  fprintf(modelfl,"SVM-multiclass Version %s\n",INST_VERSION);
  fprintf(modelfl,"%d # number of classes\n",
	  this->num_classes);
  fprintf(modelfl,"%d # number of base features\n",
	  this->num_features);
  fprintf(modelfl,"%d # loss function\n",
	  sparm->loss_function);
  fprintf(modelfl,"%ld # kernel type\n",
	  model->kernel_parm.kernel_type);
  fprintf(modelfl,"%ld # kernel parameter -d \n",
	  model->kernel_parm.poly_degree);
  fprintf(modelfl,"%.8g # kernel parameter -g \n",
	  model->kernel_parm.rbf_gamma);
  fprintf(modelfl,"%.8g # kernel parameter -s \n",
	  model->kernel_parm.coef_lin);
  fprintf(modelfl,"%.8g # kernel parameter -r \n",
	  model->kernel_parm.coef_const);
  fprintf(modelfl,"%s# kernel parameter -u \n",model->kernel_parm.custom);
  fprintf(modelfl,"%ld # highest feature index \n",model->totwords);
  fprintf(modelfl,"%ld # number of training documents \n",model->totdoc);
 
  sv_num=1;
  for(i=1;i<model->sv_num;i++) {
   for(v=model->supvec[i]->fvec;v;v=v->next) 
      sv_num++;
  }
  fprintf(modelfl,"%ld # number of support vectors plus 1 \n",sv_num);
  fprintf(modelfl,"%.8g # threshold b, each following line is a SV (starting with alpha*y)\n",model->b);

  for(i=1;i<model->sv_num;i++) {
    for(v=model->supvec[i]->fvec;v;v=v->next) {
      fprintf(modelfl,"%.32g ",model->alpha[i]*v->factor);
      fprintf(modelfl,"qid:%ld ",v->kernel_id);
      for (j=0; (v->words[j]).wnum; j++) {
	fprintf(modelfl,"%ld:%.8g ",
		(long)(v->words[j]).wnum,
		(double)(v->words[j]).weight);
      }
      if(v->userdefined)
	fprintf(modelfl,"#%s\n",v->userdefined);
      else
	fprintf(modelfl,"#\n");
    /* NOTE: this could be made more efficient by summing the
       alpha's of identical vectors before writing them to the
       file. */
    }
  }
  fclose(modelfl);
}

void        SVMMulticlass::print_struct_testing_stats(SAMPLE sample, STRUCTMODEL *sm,
				       STRUCT_LEARN_PARM *sparm, 
				       STRUCT_TEST_STATS *teststats)
{
  /* This function is called after making all test predictions in
     svm_struct_classify and allows computing and printing any kind of
     evaluation (e.g. precision/recall) you might want. You can use
     the function eval_prediction to accumulate the necessary
     statistics for each prediction. */
}

void        SVMMulticlass::eval_prediction(long exnum, EXAMPLE ex, LABEL ypred, 
			    STRUCTMODEL *sm, STRUCT_LEARN_PARM *sparm, 
			    STRUCT_TEST_STATS *teststats)
{
  /* This function allows you to accumlate statistic for how well the
     predicition matches the labeled example. It is called from
     svm_struct_classify. See also the function
     print_struct_testing_stats. */
  if(exnum == 0) { /* this is the first time the function is
		      called. So initialize the teststats */
  }
}

STRUCTMODEL SVMMulticlass::read_struct_model(const char *file, STRUCT_LEARN_PARM *sparm)
{
  /* Reads structural model sm from file file. This function is used
     only in the prediction module, not in the learning module. */
  FILE *modelfl;
  STRUCTMODEL sm;
  long i,queryid,slackid;
  double costfactor;
  long max_sv,max_words,ll,wpos;
  char *line,*comment;
  SWORD *words;
  char version_buffer[100];
  MODEL *model;

  nol_ll(file,&max_sv,&max_words,&ll); /* scan size of model file */
  max_words+=2;
  ll+=2;

  words = (SWORD *)my_malloc(sizeof(SWORD)*(max_words+10));
  line = (char *)my_malloc(sizeof(char)*ll);
  model = (MODEL *)my_malloc(sizeof(MODEL));

  if ((modelfl = fopen (file, "r")) == NULL)
  { perror (file); exit (1); }

  fscanf(modelfl,"SVM-multiclass Version %s\n",version_buffer);
  if(strcmp(version_buffer,INST_VERSION)) {
    perror ("Version of model-file does not match version of svm_struct_classify!"); 
    exit (1); 
  }
  fscanf(modelfl,"%d%*[^\n]\n", &this->num_classes);  
  fscanf(modelfl,"%d%*[^\n]\n", &this->num_features);  
  fscanf(modelfl,"%d%*[^\n]\n", &sparm->loss_function);  
  fscanf(modelfl,"%ld%*[^\n]\n", &model->kernel_parm.kernel_type);  
  fscanf(modelfl,"%ld%*[^\n]\n", &model->kernel_parm.poly_degree);
  fscanf(modelfl,"%lf%*[^\n]\n", &model->kernel_parm.rbf_gamma);
  fscanf(modelfl,"%lf%*[^\n]\n", &model->kernel_parm.coef_lin);
  fscanf(modelfl,"%lf%*[^\n]\n", &model->kernel_parm.coef_const);
  fscanf(modelfl,"%[^#]%*[^\n]\n", model->kernel_parm.custom);

  fscanf(modelfl,"%ld%*[^\n]\n", &model->totwords);
  fscanf(modelfl,"%ld%*[^\n]\n", &model->totdoc);
  fscanf(modelfl,"%ld%*[^\n]\n", &model->sv_num);
  fscanf(modelfl,"%lf%*[^\n]\n", &model->b);

  model->supvec = (DOC **)my_malloc(sizeof(DOC *)*model->sv_num);
  model->alpha = (double *)my_malloc(sizeof(double)*model->sv_num);
  model->index=NULL;
  model->lin_weights=NULL;

  for(i=1;i<model->sv_num;i++) {
    fgets(line,(int)ll,modelfl);
    if(!parse_document(line,words,&(model->alpha[i]),&queryid,&slackid,
		       &costfactor,&wpos,max_words,&comment)) {
      printf("\nParsing error while reading model file in SV %ld!\n%s",
	     i,line);
      exit(1);
    }
    model->supvec[i] = create_example(-1,0,0,0.0,
				      create_svector(words,comment,1.0));
    model->supvec[i]->fvec->kernel_id=queryid;
  }
  fclose(modelfl);
  free(line);
  free(words);
  if(verbosity>=1) {
    fprintf(stdout, " (%d support vectors read) ",(int)(model->sv_num-1));
  }
  sm.svm_model=model;
  sm.sizePsi=model->totwords;
  sm.w=NULL;
  return(sm);
}

void        SVMMulticlass::write_label(FILE *fp, LABEL yy)
{
  CLASS_LABEL *y = (CLASS_LABEL*)yy.data;

  /* Writes label y to file handle fp. */
  int i;
  fprintf(fp,"%d",y->class_label);
  if(y->scores) 
    for(i=1;i<=y->num_classes;i++)
      fprintf(fp," %f",y->scores[i]);
  fprintf(fp,"\n");
} 

void        SVMMulticlass::free_pattern(SPATTERN x) {
  /* Frees the memory of x. */
  free_example((DOC*)x.data,1);
}

void        SVMMulticlass::free_label(LABEL y) {
  /* Frees the memory of y. */
  free(y.data);
}

void        SVMMulticlass::free_struct_model(STRUCTMODEL sm) 
{
  /* Frees the memory of model. */
  /* if(sm.w) free(sm.w); */ /* this is free'd in free_model */
  if(sm.svm_model) free_model(sm.svm_model,1);
  /* add free calls for user defined data here */
}

void        SVMMulticlass::free_struct_sample(SAMPLE s)
{
  /* Frees the memory of sample s. */
  int i;
  for(i=0;i<s.n;i++) { 
    free_pattern(s.examples[i].x);
    free_label(s.examples[i].y);
  }
  free(s.examples);
}

void        SVMMulticlass::print_struct_help()
{
  /* Prints a help text that is appended to the common help text of
     svm_struct_learn. */

  printf("          none\n\n");
  printf("Based on multi-class SVM formulation described in:\n");
  printf("          K. Crammer and Y. Singer. On the Algorithmic Implementation of\n");
  printf("          Multi-class SVMs, JMLR, 2001.\n");
}

void         SVMMulticlass::parse_struct_parameters(STRUCT_LEARN_PARM *sparm)
{
  /* Parses the command line parameters that start with -- */
  int i;

  for(i=0;(i<sparm->custom_argc) && ((sparm->custom_argv[i])[0] == '-');i++) {
    switch ((sparm->custom_argv[i])[2]) 
      { 
      case 'a': i++; /* strcpy(learn_parm->alphafile,argv[i]); */ break;
      case 'e': i++; /* sparm->epsilon=atof(sparm->custom_argv[i]); */ break;
      case 'k': i++; /* sparm->newconstretrain=atol(sparm->custom_argv[i]); */ break;
      }
  }
}

void        SVMMulticlass::print_struct_help_classify()
{
  /* Prints a help text that is appended to the common help text of
     svm_struct_classify. */
}

void        SVMMulticlass::parse_struct_parameters_classify(STRUCT_LEARN_PARM *sparm)
{
  /* Parses the command line parameters that start with -- for the
     classification module */
  int i;

  for(i=0;(i<sparm->custom_argc) && ((sparm->custom_argv[i])[0] == '-');i++) {
    switch ((sparm->custom_argv[i])[2]) 
      { 
      /* case 'x': i++; strcpy(xvalue,sparm->custom_argv[i]); break; */
      default: printf("\nUnrecognized option %s!\n\n",sparm->custom_argv[i]);
	       exit(0);
      }
  }
}

SAMPLE SVMMulticlass::struct_examples_from_array(float *feat, int *classes, int num_examples, int num_features, int num_classes, unsigned char *mask) {
  SAMPLE   sample;
  EXAMPLE *examples=(EXAMPLE *)my_malloc(sizeof(EXAMPLE)*num_examples);
  int i, j;
  float *ptr = feat;
  SWORD *words;
  CLASS_LABEL *y;

  sample.n = 0;
  sample.examples = examples;
  for(i = 0; i < num_examples; i++, ptr += num_features) { 
    if(mask && !mask[i]) 
      continue; 

    words = (SWORD*)my_malloc(sizeof(SWORD)*(num_features+10));
    for(j = 0; j < num_features; j++) {
      words[j].wnum = j+1;
      words[j].weight = ptr[j];
    }
    words[num_features].wnum = 0;
    words[num_features].weight = 0;

    examples[sample.n].x.data = create_example(i, 0, 0, 1, create_svector(words,NULL,1.0));
    examples[sample.n].y.data = y = (CLASS_LABEL*)malloc(sizeof(CLASS_LABEL));
    y->scores=NULL;
    y->num_classes=num_classes;    
    y->class_label = classes ? classes[i] : 0;
    sample.n++;
  }

  return sample;
}

void SVMMulticlass::svm_struct_classify_from_array(float *feat, int *preds, float *scores, int num_examples, STRUCTMODEL *model, STRUCT_LEARN_PARM *sparm) {
  SAMPLE test = struct_examples_from_array(feat, NULL, num_examples, this->num_features, this->num_classes, NULL);
  LABEL y;
  int i;

  for(i = 0; i < num_examples; i++) {
    y=classify_struct_example(test.examples[i].x, model, sparm);
    preds[i] = ((CLASS_LABEL*)y.data)->class_label;
    if(scores) {
      scores[0] = 0;
      memcpy(scores+1, ((CLASS_LABEL*)y.data)->scores, this->num_classes*sizeof(float));
    }
  }
  free_struct_sample(test);
}


void SVMMulticlass::svm_struct_learn_from_array(STRUCT_LEARN_PARM *struct_parm, STRUCTMODEL *structmodel, 
				 float *feat, int *classes, unsigned char *mask, int num_examples, int num_features, int num_classes, int argc, const char *argv[]) {
  SAMPLE sample;  /* training sample */
  LEARN_PARM learn_parm;
  KERNEL_PARM kernel_parm;
  int alg_type;

  svm_struct_learn_api_init(argc,argv);
  struct_learn_read_input_parameters(argc,argv,&verbosity,
				     &struct_verbosity,struct_parm,&learn_parm,
				     &kernel_parm,&alg_type);

  /* read the training examples */
  sample=struct_examples_from_array(feat, classes, num_examples, num_features, num_classes, mask);
  if(struct_verbosity>=1) {
    printf("done\n"); fflush(stdout);
  }

  /* Do the learning and return structmodel. */
  if(alg_type == 0)
    svm_learn_struct(sample,struct_parm,&learn_parm,&kernel_parm,structmodel,NSLACK_ALG,this);
  else if(alg_type == 1)
    svm_learn_struct(sample,struct_parm,&learn_parm,&kernel_parm,structmodel,NSLACK_SHRINK_ALG,this);
  else if(alg_type == 2)
    svm_learn_struct_joint(sample,struct_parm,&learn_parm,&kernel_parm,structmodel,ONESLACK_PRIMAL_ALG,this);
  else if(alg_type == 3)
    svm_learn_struct_joint(sample,struct_parm,&learn_parm,&kernel_parm,structmodel,ONESLACK_DUAL_ALG,this);
  else if(alg_type == 4)
    svm_learn_struct_joint(sample,struct_parm,&learn_parm,&kernel_parm,structmodel,ONESLACK_DUAL_CACHE_ALG,this);
  else if(alg_type == 9)
    svm_learn_struct_joint_custom(sample,struct_parm,&learn_parm,&kernel_parm,structmodel,this);
  else
    exit(1);

  /* Warning: The model contains references to the original data 'docs'.
     If you want to free the original data, and only keep the model, you 
     have to make a deep copy of 'model'. */
  if(struct_verbosity>=1) {
    printf("Writing learned model...");fflush(stdout);
  }
  if(struct_verbosity>=1) {
    printf("done\n");fflush(stdout);
  }

  free_struct_sample(sample);
  //free_struct_model(structmodel);

  svm_struct_learn_api_exit();
}

void SVMMulticlass::struct_learn_read_input_parameters(int argc,const char *argv[],
			   long *verbosity,long *struct_verbosity, 
			   STRUCT_LEARN_PARM *struct_parm,
			   LEARN_PARM *learn_parm, KERNEL_PARM *kernel_parm,
			   int *alg_type)
{
  long i;
  char type[100];
  
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

  for(i=1;(i<argc) && ((argv[i])[0] == '-');i++) {
    switch ((argv[i])[1]) 
      { 
      case '?': struct_learn_print_help(); exit(0);
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
      default: printf("\nUnrecognized option %s!\n\n",argv[i]);
	       struct_learn_print_help();
	       exit(0);
      }
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
    struct_learn_wait_any_key();
    struct_learn_print_help();
    exit(0);
  }    
  if((learn_parm->svm_maxqpsize<2)) {
    printf("\nMaximum size of QP-subproblems not in valid range: %ld [2..]\n",learn_parm->svm_maxqpsize); 
    struct_learn_wait_any_key();
    struct_learn_print_help();
    exit(0);
  }
  if((learn_parm->svm_maxqpsize<learn_parm->svm_newvarsinqp)) {
    printf("\nMaximum size of QP-subproblems [%ld] must be larger than the number of\n",learn_parm->svm_maxqpsize); 
    printf("new variables [%ld] entering the working set in each iteration.\n",learn_parm->svm_newvarsinqp); 
    struct_learn_wait_any_key();
    struct_learn_print_help();
    exit(0);
  }
  if(learn_parm->svm_iter_to_shrink<1) {
    printf("\nMaximum number of iterations for shrinking not in valid range: %ld [1,..]\n",learn_parm->svm_iter_to_shrink);
    struct_learn_wait_any_key();
    struct_learn_print_help();
    exit(0);
  }
  if(struct_parm->C<0) {
    printf("\nYou have to specify a value for the parameter '-c' (C>0)!\n\n");
    struct_learn_wait_any_key();
    struct_learn_print_help();
    exit(0);
  }
  if(((*alg_type) < 0) || (((*alg_type) > 5) && ((*alg_type) != 9))) {
    printf("\nAlgorithm type must be either '0', '1', '2', '3', '4', or '9'!\n\n");
    struct_learn_wait_any_key();
    struct_learn_print_help();
    exit(0);
  }
  if(learn_parm->transduction_posratio>1) {
    printf("\nThe fraction of unlabeled examples to classify as positives must\n");
    printf("be less than 1.0 !!!\n\n");
    struct_learn_wait_any_key();
    struct_learn_print_help();
    exit(0);
  }
  if(learn_parm->svm_costratio<=0) {
    printf("\nThe COSTRATIO parameter must be greater than zero!\n\n");
    struct_learn_wait_any_key();
    struct_learn_print_help();
    exit(0);
  }
  if(struct_parm->epsilon<=0) {
    printf("\nThe epsilon parameter must be greater than zero!\n\n");
    struct_learn_wait_any_key();
    struct_learn_print_help();
    exit(0);
  }
  if((struct_parm->ccache_size<=0) && ((*alg_type) == 4)) {
    printf("\nThe cache size must be at least 1!\n\n");
    struct_learn_wait_any_key();
    struct_learn_print_help();
    exit(0);
  }
  if(((struct_parm->batch_size<=0) || (struct_parm->batch_size>100))  
     && ((*alg_type) == 4)) {
    printf("\nThe batch size must be in the interval ]0,100]!\n\n");
    struct_learn_wait_any_key();
    struct_learn_print_help();
    exit(0);
  }
  if((struct_parm->slack_norm<1) || (struct_parm->slack_norm>2)) {
    printf("\nThe norm of the slacks must be either 1 (L1-norm) or 2 (L2-norm)!\n\n");
    struct_learn_wait_any_key();
    struct_learn_print_help();
    exit(0);
  }
  if((struct_parm->loss_type != SLACK_RESCALING) 
     && (struct_parm->loss_type != MARGIN_RESCALING)) {
    printf("\nThe loss type must be either 1 (slack rescaling) or 2 (margin rescaling)!\n\n");
    struct_learn_wait_any_key();
    struct_learn_print_help();
    exit(0);
  }
  if(learn_parm->rho<0) {
    printf("\nThe parameter rho for xi/alpha-estimates and leave-one-out pruning must\n");
    printf("be greater than zero (typically 1.0 or 2.0, see T. Joachims, Estimating the\n");
    printf("Generalization Performance of an SVM Efficiently, ICML, 2000.)!\n\n");
    struct_learn_wait_any_key();
    struct_learn_print_help();
    exit(0);
  }
  if((learn_parm->xa_depth<0) || (learn_parm->xa_depth>100)) {
    printf("\nThe parameter depth for ext. xi/alpha-estimates must be in [0..100] (zero\n");
    printf("for switching to the conventional xa/estimates described in T. Joachims,\n");
    printf("Estimating the Generalization Performance of an SVM Efficiently, ICML, 2000.)\n");
    struct_learn_wait_any_key();
    struct_learn_print_help();
    exit(0);
  }

  parse_struct_parameters(struct_parm);
}

void SVMMulticlass::struct_learn_wait_any_key()
{
  printf("\n(more)\n");
  (void)getc(stdin);
}

void SVMMulticlass::struct_learn_print_help()
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
  print_struct_help();
  struct_learn_wait_any_key();

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

