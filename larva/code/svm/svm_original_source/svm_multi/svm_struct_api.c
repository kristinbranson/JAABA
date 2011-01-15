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
#include "svm_struct_api.h"

void        svm_struct_learn_api_init(int argc, char* argv[])
{
  /* Called in learning part before anything else is done to allow
     any initializations that might be necessary. */
}

void        svm_struct_learn_api_exit()
{
  /* Called in learning part at the very end to allow any clean-up
     that might be necessary. */
}

void        svm_struct_classify_api_init(int argc, char* argv[])
{
  /* Called in prediction part before anything else is done to allow
     any initializations that might be necessary. */
}

void        svm_struct_classify_api_exit()
{
  /* Called in prediction part at the very end to allow any clean-up
     that might be necessary. */
}

SAMPLE      read_struct_examples(char *file, STRUCT_LEARN_PARM *sparm)
{
  /* Reads training examples and returns them in sample. The number of
     examples must be written into sample.n */
  SAMPLE   sample;  /* sample */
  EXAMPLE  *examples;
  long     n;       /* number of examples */
  DOC **docs;       /* examples in original SVM-light format */ 
  double *target;
  long totwords,i,num_classes=0;

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
    examples[i].x.doc=docs[i];
    examples[i].y.class=target[i]+0.1;
    examples[i].y.scores=NULL;
    examples[i].y.num_classes=num_classes;
  }
  free(target);
  free(docs);
  sample.n=n;
  sample.examples=examples;

  if(struct_verbosity>=0)
    printf(" (%d examples) ",sample.n);
  return(sample);
}

void        init_struct_model(SAMPLE sample, STRUCTMODEL *sm, 
			      STRUCT_LEARN_PARM *sparm, LEARN_PARM *lparm, 
			      KERNEL_PARM *kparm)
{
  /* Initialize structmodel sm. The weight vector w does not need to be
     initialized, but you need to provide the maximum size of the
     feature space in sizePsi. This is the maximum number of different
     weights that can be learned. Later, the weight vector w will
     contain the learned weights for the model. */
  long i,totwords=0;
  WORD *w;

  sparm->num_classes=1;
  for(i=0;i<sample.n;i++)     /* find highest class label */
    if(sparm->num_classes < (sample.examples[i].y.class+0.1)) 
      sparm->num_classes=sample.examples[i].y.class+0.1;
  for(i=0;i<sample.n;i++)     /* find highest feature number */
    for(w=sample.examples[i].x.doc->fvec->words;w->wnum;w++) 
      if(totwords < w->wnum) 
	totwords=w->wnum;
  sparm->num_features=totwords;
  if(struct_verbosity>=0)
    printf("Training set properties: %d features, %d classes\n",
	   sparm->num_features,sparm->num_classes);
  sm->sizePsi=sparm->num_features*sparm->num_classes;
  if(struct_verbosity>=2)
    printf("Size of Phi: %ld\n",sm->sizePsi);
}

CONSTSET    init_struct_constraints(SAMPLE sample, STRUCTMODEL *sm, 
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
  WORD     words[2];

  if(1) { /* normal case: start with empty set of constraints */
    c.lhs=NULL;
    c.rhs=NULL;
    c.m=0;
  }
  else { /* add constraints so that all learned weights are
            positive. WARNING: Currently, they are positive only up to
            precision epsilon set by -e. */
    c.lhs=my_malloc(sizeof(DOC *)*sizePsi);
    c.rhs=my_malloc(sizeof(double)*sizePsi);
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

LABEL       classify_struct_example(PATTERN x, STRUCTMODEL *sm, 
				    STRUCT_LEARN_PARM *sparm)
{
  /* Finds the label yhat for pattern x that scores the highest
     according to the linear evaluation function in sm, especially the
     weights sm.w. The returned label is taken as the prediction of sm
     for the pattern x. The weights correspond to the features defined
     by psi() and range from index 1 to index sm->sizePsi. If the
     function cannot find a label, it shall return an empty label as
     recognized by the function empty_label(y). */
  LABEL y;
  DOC doc;
  long class,bestclass=-1,first=1,j;
  double score,bestscore=-1;
  WORD *words;

  doc=*(x.doc);
  y.scores=(double *)my_malloc(sizeof(double)*(sparm->num_classes+1));
  y.num_classes=sparm->num_classes;
  words=doc.fvec->words;
  for(j=0;(words[j]).wnum != 0;j++) {       /* Check if feature numbers   */
    if((words[j]).wnum>sparm->num_features) /* are not larger than in     */
      (words[j]).wnum=0;                    /* model. Remove feature if   */
  }                                         /* necessary.                 */
  for(class=1;class<=sparm->num_classes;class++) {
    y.class=class;
    doc.fvec=psi(x,y,sm,sparm);
    score=classify_example(sm->svm_model,&doc);
    free_svector(doc.fvec);
    y.scores[class]=score;
    if((bestscore<score)  || (first)) {
      bestscore=score;
      bestclass=class;
      first=0;
    }
  }
  y.class=bestclass;
  return(y);
}

LABEL       find_most_violated_constraint_slackrescaling(PATTERN x, LABEL y, 
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
  LABEL ybar;
  DOC doc;
  long class,bestclass=-1,first=1;
  double score,score_y,score_ybar,bestscore=-1;

  /* NOTE: This function could be made much more efficient by not
     always computing a new PSI vector. */
  doc=*(x.doc);
  doc.fvec=psi(x,y,sm,sparm);
  score_y=classify_example(sm->svm_model,&doc);
  free_svector(doc.fvec);

  ybar.scores=NULL;
  ybar.num_classes=sparm->num_classes;
  for(class=1;class<=sparm->num_classes;class++) {
    ybar.class=class;
    doc.fvec=psi(x,ybar,sm,sparm);
    score_ybar=classify_example(sm->svm_model,&doc);
    free_svector(doc.fvec);
    score=loss(y,ybar,sparm)*(1.0-score_y+score_ybar);
    if((bestscore<score)  || (first)) {
      bestscore=score;
      bestclass=class;
      first=0;
    }
  }
  if(bestclass == -1) 
    printf("ERROR: Only one class\n");
  ybar.class=bestclass;
  if(struct_verbosity>=3)
    printf("[%ld:%.2f] ",bestclass,bestscore);
  return(ybar);
}

LABEL       find_most_violated_constraint_marginrescaling(PATTERN x, LABEL y, 
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
  LABEL ybar;
  DOC doc;
  long class,bestclass=-1,first=1;
  double score,bestscore=-1;

  /* NOTE: This function could be made much more efficient by not
     always computing a new PSI vector. */
  doc=*(x.doc);
  ybar.scores=NULL;
  ybar.num_classes=sparm->num_classes;
  for(class=1;class<=sparm->num_classes;class++) {
    ybar.class=class;
    doc.fvec=psi(x,ybar,sm,sparm);         
    score=classify_example(sm->svm_model,&doc);
    free_svector(doc.fvec);
    score+=loss(y,ybar,sparm);
    if((bestscore<score)  || (first)) {
      bestscore=score;
      bestclass=class;
      first=0;
    }
  }
  if(bestclass == -1) 
    printf("ERROR: Only one class\n");
  ybar.class=bestclass;
  if(struct_verbosity>=3)
    printf("[%ld:%.2f] ",bestclass,bestscore);
  return(ybar);
}

int         empty_label(LABEL y)
{
  /* Returns true, if y is an empty label. An empty label might be
     returned by find_most_violated_constraint_???(x, y, sm) if there
     is no incorrect label that can be found for x, or if it is unable
     to label x at all */
  return(y.class<0.9);
}

SVECTOR     *psi(PATTERN x, LABEL y, STRUCTMODEL *sm,
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

  /* shift the feature numbers to the position of weight vector of class y */
  fvec=shift_s(x.doc->fvec,(y.class-1)*sparm->num_features);

  /* The following makes sure that the weight vectors for each class
     are treated separately when kernels are used . */
  fvec->kernel_id=y.class;

  return(fvec);
}

double      loss(LABEL y, LABEL ybar, STRUCT_LEARN_PARM *sparm)
{
  /* loss for correct label y and predicted label ybar. The loss for
     y==ybar has to be zero. sparm->loss_function is set with the -l option. */
  if(sparm->loss_function == 0) { /* type 0 loss: 0/1 loss */
    if(y.class == ybar.class)     /* return 0, if y==ybar. return 100 else */
      return(0);
    else
      return(100);
  }
  if(sparm->loss_function == 1) { /* type 1 loss: squared difference */
    return((y.class-ybar.class)*(y.class-ybar.class));
  }
  else {
    /* Put your code for different loss functions here. But then
       find_most_violated_constraint_???(x, y, sm) has to return the
       highest scoring label with the largest loss. */
    printf("Unkown loss function\n");
    exit(1);
  }
}

int         finalize_iteration(double ceps, int cached_constraint,
			       SAMPLE sample, STRUCTMODEL *sm,
			       CONSTSET cset, double *alpha, 
			       STRUCT_LEARN_PARM *sparm)
{
  /* This function is called just before the end of each cutting plane iteration. ceps is the amount by which the most violated constraint found in the current iteration was violated. cached_constraint is true if the added constraint was constructed from the cache. If the return value is FALSE, then the algorithm is allowed to terminate. If it is TRUE, the algorithm will keep iterating even if the desired precision sparm->epsilon is already reached. */
  return(0);
}

void        print_struct_learning_stats(SAMPLE sample, STRUCTMODEL *sm,
					CONSTSET cset, double *alpha, 
					STRUCT_LEARN_PARM *sparm)
{
  /* This function is called after training and allows final touches to
     the model sm. But primarly it allows computing and printing any
     kind of statistic (e.g. training error) you might want. */

  /* Replace SV with single weight vector */
  MODEL *model=sm->svm_model;
  if(model->kernel_parm.kernel_type == LINEAR) {
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

void        write_struct_model(char *file, STRUCTMODEL *sm, 
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
	  sparm->num_classes);
  fprintf(modelfl,"%d # number of base features\n",
	  sparm->num_features);
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

void        print_struct_testing_stats(SAMPLE sample, STRUCTMODEL *sm,
				       STRUCT_LEARN_PARM *sparm, 
				       STRUCT_TEST_STATS *teststats)
{
  /* This function is called after making all test predictions in
     svm_struct_classify and allows computing and printing any kind of
     evaluation (e.g. precision/recall) you might want. You can use
     the function eval_prediction to accumulate the necessary
     statistics for each prediction. */
}

void        eval_prediction(long exnum, EXAMPLE ex, LABEL ypred, 
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

STRUCTMODEL read_struct_model(char *file, STRUCT_LEARN_PARM *sparm)
{
  /* Reads structural model sm from file file. This function is used
     only in the prediction module, not in the learning module. */
  FILE *modelfl;
  STRUCTMODEL sm;
  long i,queryid,slackid;
  double costfactor;
  long max_sv,max_words,ll,wpos;
  char *line,*comment;
  WORD *words;
  char version_buffer[100];
  MODEL *model;

  nol_ll(file,&max_sv,&max_words,&ll); /* scan size of model file */
  max_words+=2;
  ll+=2;

  words = (WORD *)my_malloc(sizeof(WORD)*(max_words+10));
  line = (char *)my_malloc(sizeof(char)*ll);
  model = (MODEL *)my_malloc(sizeof(MODEL));

  if ((modelfl = fopen (file, "r")) == NULL)
  { perror (file); exit (1); }

  fscanf(modelfl,"SVM-multiclass Version %s\n",version_buffer);
  if(strcmp(version_buffer,INST_VERSION)) {
    perror ("Version of model-file does not match version of svm_struct_classify!"); 
    exit (1); 
  }
  fscanf(modelfl,"%d%*[^\n]\n", &sparm->num_classes);  
  fscanf(modelfl,"%d%*[^\n]\n", &sparm->num_features);  
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

void        write_label(FILE *fp, LABEL y)
{
  /* Writes label y to file handle fp. */
  int i;
  fprintf(fp,"%d",y.class);
  if(y.scores) 
    for(i=1;i<=y.num_classes;i++)
      fprintf(fp," %f",y.scores[i]);
  fprintf(fp,"\n");
} 

void        free_pattern(PATTERN x) {
  /* Frees the memory of x. */
  free_example(x.doc,1);
}

void        free_label(LABEL y) {
  /* Frees the memory of y. */
  if(y.scores) free(y.scores);
}

void        free_struct_model(STRUCTMODEL sm) 
{
  /* Frees the memory of model. */
  /* if(sm.w) free(sm.w); */ /* this is free'd in free_model */
  if(sm.svm_model) free_model(sm.svm_model,1);
  /* add free calls for user defined data here */
}

void        free_struct_sample(SAMPLE s)
{
  /* Frees the memory of sample s. */
  int i;
  for(i=0;i<s.n;i++) { 
    free_pattern(s.examples[i].x);
    free_label(s.examples[i].y);
  }
  free(s.examples);
}

void        print_struct_help()
{
  /* Prints a help text that is appended to the common help text of
     svm_struct_learn. */

  printf("          none\n\n");
  printf("Based on multi-class SVM formulation described in:\n");
  printf("          K. Crammer and Y. Singer. On the Algorithmic Implementation of\n");
  printf("          Multi-class SVMs, JMLR, 2001.\n");
}

void         parse_struct_parameters(STRUCT_LEARN_PARM *sparm)
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

void        print_struct_help_classify()
{
  /* Prints a help text that is appended to the common help text of
     svm_struct_classify. */
}

void        parse_struct_parameters_classify(STRUCT_LEARN_PARM *sparm)
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

