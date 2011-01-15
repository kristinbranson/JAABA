/************************************************************************/
/*                                                                      */
/*   svm_struct_latent_cccp.c                                           */
/*                                                                      */
/*   Main Optimization Code for Latent SVM^struct                       */
/*                                                                      */
/*   Author: Chun-Nam Yu                                                */
/*   Date: 30.Sep.08                                                    */
/*                                                                      */
/*   This software is available for non-commercial use only. It must    */
/*   not be modified and distributed without prior permission of the    */
/*   author. The author is not responsible for implications from the    */
/*   use of this software.                                              */
/*                                                                      */
/************************************************************************/

#include <stdio.h>
#include <assert.h>
//#include "./svm_light/svm_common.h"
//#include "svm_struct_latent_api_types.h"
#include "svm_struct_latent_api.h"
#include "./svm_light/svm_learn.h"


#define ALPHA_THRESHOLD 1E-14
#define IDLE_ITER 20
#define CLEANUP_CHECK 100
#define STOP_PREC 1E-2
#define UPDATE_BOUND 3

#define MAX_OUTER_ITER 400

#define MAX(x,y) ((x) < (y) ? (y) : (x))
#define MIN(x,y) ((x) > (y) ? (y) : (x))

#define DEBUG_LEVEL 0

void my_read_input_parameters(int argc, char* argv[], char *trainfile, char *modelfile,
			      LEARN_PARM *learn_parm, KERNEL_PARM *kernel_parm, STRUCT_LEARN_PARM *struct_parm);

void my_wait_any_key();

int resize_cleanup(int size_active, int *idle, double *alpha, double *delta, double *gammaG0, double *proximal_rhs, double **G, DOC **dXc, double *cut_error);

double sprod_nn(double *a, double *b, long n) {
  double ans=0.0;
  long i;
  for (i=1;i<n+1;i++) {
    ans+=a[i]*b[i];
  }
  return(ans);
}

void add_vector_nn(double *w, double *dense_x, long n, double factor) {
  long i;
  for (i=1;i<n+1;i++) {
    w[i]+=factor*dense_x[i];
  }
}

double* add_list_nn(SVECTOR *a, long totwords) 
     /* computes the linear combination of the SVECTOR list weighted
	by the factor of each SVECTOR. assumes that the number of
	features is small compared to the number of elements in the
	list */
{
    SVECTOR *f;
    long i;
    double *sum;

    sum=create_nvector(totwords);

    for(i=0;i<=totwords;i++) 
      sum[i]=0;

    for(f=a;f;f=f->next)  
      add_vector_ns(sum,f,f->factor);

    return(sum);
}


SVECTOR* find_cutting_plane(EXAMPLE *ex, SVECTOR **fycache, double *margin, long m, STRUCTMODEL *sm, STRUCT_LEARN_PARM *sparm) {

  long i;
  SVECTOR *f, *fy, *fybar, *lhs;
  LABEL       ybar;
  LATENT_VAR hbar;
  double lossval;
  double *new_constraint;

  long l,k;
  SVECTOR *fvec;
  WORD *words;  

  /* find cutting plane */
  lhs = NULL;
  *margin = 0;
  for (i=0;i<m;i++) {
    find_most_violated_constraint_marginrescaling(ex[i].x, ex[i].y, &ybar, &hbar, sm, sparm);
    /* get difference vector */
    fy = copy_svector(fycache[i]);
    fybar = psi(ex[i].x,ybar,hbar,sm,sparm);
    lossval = loss(ex[i].y,ybar,hbar,sparm);
    free_label(ybar);
    free_latent_var(hbar);

    /* scale difference vector */
    for (f=fy;f;f=f->next) {
      f->factor*=1.0/m;
      //f->factor*=ex[i].x.example_cost/m;
    }
    for (f=fybar;f;f=f->next) {
      f->factor*=-1.0/m;
      //f->factor*=-ex[i].x.example_cost/m;
    }
    /* add ybar to constraint */
    append_svector_list(fy,lhs);
    append_svector_list(fybar,fy);
    lhs = fybar;
    *margin+=lossval/m;
    //*margin+=lossval*ex[i].x.example_cost/m;
  }
  /* compact the linear representation */
  new_constraint = add_list_nn(lhs, sm->sizePsi);
  free_svector(lhs);

  l=0;
  for (i=1;i<sm->sizePsi+1;i++) {
    if (fabs(new_constraint[i])>1E-10) l++; // non-zero
  }
  words = (WORD*)my_malloc(sizeof(WORD)*(l+1)); 
  assert(words!=NULL);
  k=0;
  for (i=1;i<sm->sizePsi+1;i++) {
    if (fabs(new_constraint[i])>1E-10) {
      words[k].wnum = i;
      words[k].weight = new_constraint[i]; 
      k++;
    }
  }
  words[k].wnum = 0;
  words[k].weight = 0.0;
  fvec = create_svector(words,"",1);

  free(words);
  free(new_constraint);

  return(fvec); 

}


double cutting_plane_algorithm(double *w, long m, int MAX_ITER, double C, double epsilon, SVECTOR **fycache, EXAMPLE *ex, STRUCTMODEL *sm, STRUCT_LEARN_PARM *sparm) {
  long i,j;
  double xi;
  double *alpha;
  double **G; /* Gram matrix */
  DOC **dXc; /* constraint matrix */
  double *delta; /* rhs of constraints */
  SVECTOR *new_constraint;
  double dual_obj, alphasum;
  int iter, size_active; 
  double value;
  int r;
  int *idle; /* for cleaning up */
  double margin;
  double primal_obj;
  double *proximal_rhs;
  double *gammaG0=NULL;
  double min_rho = 0.001;
  double max_rho;
  double serious_counter=0;
  double rho = 1.0; /* temporarily set it to 1 first */

  double expected_descent, primal_obj_b=-1, reg_master_obj;
  int null_step=1;
  double *w_b;
  double kappa=0.1;
  double temp_var;
  double proximal_term, primal_lower_bound;

  double v_k; 
  double obj_difference; 
  double *cut_error; // cut_error[i] = alpha_{k,i} at current center x_k
  double sigma_k; 
  double m2 = 0.2;
  double m3 = 0.9;
  double gTd; 
  double last_sigma_k=0; 

  double initial_primal_obj;
  int suff_decrease_cond=0;
  double decrease_proportion = 0.2; // start from 0.2 first 

  double z_k_norm;
  double last_z_k_norm=0;


  /* set parameters for hideo solver */
  LEARN_PARM lparm;
  KERNEL_PARM kparm;
  MODEL *svmModel=NULL;
  lparm.biased_hyperplane = 0;
  lparm.epsilon_crit = MIN(epsilon,0.001);
  lparm.svm_c = C;
  lparm.sharedslack = 1;
  kparm.kernel_type = LINEAR;

  lparm.remove_inconsistent=0;
  lparm.skip_final_opt_check=0;
  lparm.svm_maxqpsize=10;
  lparm.svm_newvarsinqp=0;
  lparm.svm_iter_to_shrink=-9999;
  lparm.maxiter=100000;
  lparm.kernel_cache_size=40;
  lparm.eps = epsilon; 
  lparm.transduction_posratio=-1.0;
  lparm.svm_costratio=1.0;
  lparm.svm_costratio_unlab=1.0;
  lparm.svm_unlabbound=1E-5;
  lparm.epsilon_a=1E-10;  /* changed from 1e-15 */
  lparm.compute_loo=0;
  lparm.rho=1.0;
  lparm.xa_depth=0;
  strcpy(lparm.alphafile,"");
  kparm.poly_degree=3;
  kparm.rbf_gamma=1.0;
  kparm.coef_lin=1;
  kparm.coef_const=1;
  strcpy(kparm.custom,"empty");
 


  w_b = create_nvector(sm->sizePsi);
  clear_nvector(w_b,sm->sizePsi);
  /* warm start */
  for (i=1;i<sm->sizePsi+1;i++) {
    w_b[i] = w[i];
  }

  iter = 0;
  size_active = 0;
  xi = 0.0;
  alpha = NULL;
  G = NULL;
  dXc = NULL;
  delta = NULL;
  idle = NULL;

  proximal_rhs = NULL;
  cut_error = NULL; 

  new_constraint = find_cutting_plane(ex, fycache, &margin, m, sm, sparm);
  value = margin - sprod_ns(w, new_constraint);
	
  primal_obj_b = 0.5*sprod_nn(w_b,w_b,sm->sizePsi)+C*value;
  primal_obj = 0.5*sprod_nn(w,w,sm->sizePsi)+C*value;
  primal_lower_bound = 0;
  expected_descent = -primal_obj_b;
  initial_primal_obj = primal_obj_b; 

  max_rho = C; 

  printf("Running CCCP inner loop solver: "); fflush(stdout); 

  while ((!suff_decrease_cond)&&(expected_descent<-epsilon)&&(iter<MAX_ITER)) { 
    iter+=1;
    size_active+=1;

#if (DEBUG_LEVEL>0)
    printf("ITER %d\n", iter); 
#endif
    printf("."); fflush(stdout); 

    /* add  constraint */
    dXc = (DOC**)realloc(dXc, sizeof(DOC*)*size_active);
    assert(dXc!=NULL);
    dXc[size_active-1] = (DOC*)malloc(sizeof(DOC));
    dXc[size_active-1]->fvec = new_constraint; 
    dXc[size_active-1]->slackid = 1; // only one common slackid (one-slack)
    dXc[size_active-1]->costfactor = 1.0;

    delta = (double*)realloc(delta, sizeof(double)*size_active);
    assert(delta!=NULL);
    delta[size_active-1] = margin;
    alpha = (double*)realloc(alpha, sizeof(double)*size_active);
    assert(alpha!=NULL);
    alpha[size_active-1] = 0.0;
    idle = (int*)realloc(idle, sizeof(int)*size_active);
    assert(idle!=NULL); 
    idle[size_active-1] = 0;
    /* proximal point */
    proximal_rhs = (double*)realloc(proximal_rhs, sizeof(double)*size_active);
    assert(proximal_rhs!=NULL); 
    cut_error = (double*)realloc(cut_error, sizeof(double)*size_active); 
    assert(cut_error!=NULL); 
    // note g_i = - new_constraint
    cut_error[size_active-1] = C*(sprod_ns(w_b, new_constraint) - sprod_ns(w, new_constraint)); 
    cut_error[size_active-1] += (primal_obj_b - 0.5*sprod_nn(w_b,w_b,sm->sizePsi)); 
    cut_error[size_active-1] -= (primal_obj - 0.5*sprod_nn(w,w,sm->sizePsi)); 

    gammaG0 = (double*)realloc(gammaG0, sizeof(double)*size_active);
    assert(gammaG0!=NULL);
      
    /* update Gram matrix */
    G = (double**)realloc(G, sizeof(double*)*size_active);
    assert(G!=NULL);
    G[size_active-1] = NULL;
    for (j=0;j<size_active;j++) {
      G[j] = (double*)realloc(G[j], sizeof(double)*size_active);
      assert(G[j]!=NULL);
    }
    for (j=0;j<size_active-1;j++) {
      G[size_active-1][j] = sprod_ss(dXc[size_active-1]->fvec, dXc[j]->fvec);
      G[j][size_active-1] = G[size_active-1][j];
    }
    G[size_active-1][size_active-1] = sprod_ss(dXc[size_active-1]->fvec,dXc[size_active-1]->fvec);

	
    /* update gammaG0 */
    if (null_step==1) {
      gammaG0[size_active-1] = sprod_ns(w_b, dXc[size_active-1]->fvec);
    } else {
      for (i=0;i<size_active;i++) {
	gammaG0[i] = sprod_ns(w_b, dXc[i]->fvec); 
      }
    }

     /* update proximal_rhs */
    for (i=0;i<size_active;i++) {
      proximal_rhs[i] = (1+rho)*delta[i] - rho*gammaG0[i];
    }


    /* solve QP to update alpha */
    //dual_obj = 0; 
    //r = mosek_qp_optimize(G, proximal_rhs, alpha, (long) size_active, C, &dual_obj,rho);
    if (size_active>1) {
      if (svmModel!=NULL) free_model(svmModel,0);
      svmModel = (MODEL*)my_malloc(sizeof(MODEL)); 
      svm_learn_optimization(dXc,proximal_rhs,size_active,sm->sizePsi,&lparm,&kparm,NULL,svmModel,alpha); 
    } else {
      assert(size_active==1); 
      alpha[0] = C; 
    }
    /* DEBUG */
    //printf("r: %d\n", r); fflush(stdout);
    /* END DEBUG */

    clear_nvector(w,sm->sizePsi);
    for (j=0;j<size_active;j++) {
      if (alpha[j]>C*ALPHA_THRESHOLD) {
	add_vector_ns(w,dXc[j]->fvec,alpha[j]/(1+rho));
      }
    }

    /* compute dual obj */
    dual_obj = +0.5*(1+rho)*sprod_nn(w,w,sm->sizePsi);
    for (j=0;j<size_active;j++) {
      dual_obj -= proximal_rhs[j]/(1+rho)*alpha[j];
    }

    z_k_norm = sqrt(sprod_nn(w,w,sm->sizePsi)); 

    add_vector_nn(w, w_b, sm->sizePsi, rho/(1+rho));

    
    /* detect if step size too small */
    sigma_k = 0; 
    alphasum = 0; 
    for (j=0;j<size_active;j++) {
      sigma_k += alpha[j]*cut_error[j]; 
      alphasum+=alpha[j]; 
    }
    sigma_k/=C; 
    gTd = -C*(sprod_ns(w,new_constraint) - sprod_ns(w_b,new_constraint));

#if (DEBUG_LEVEL>0)
    for (j=0;j<size_active;j++) {
      printf("alpha[%d]: %.8g, cut_error[%d]: %.8g\n", j, alpha[j], j, cut_error[j]);
    }
    printf("sigma_k: %.8g\n", sigma_k); 
    printf("alphasum: %.8g\n", alphasum);
    printf("g^T d: %.8g\n", gTd); 
    fflush(stdout); 
#endif


    /* update cleanup information */
    for (j=0;j<size_active;j++) {
      if (alpha[j]<ALPHA_THRESHOLD*C) {
	idle[j]++;
      } else {
        idle[j]=0;
      }
    }

    new_constraint = find_cutting_plane(ex, fycache, &margin, m, sm, sparm);
    value = margin - sprod_ns(w, new_constraint);

    /* print primal objective */
    primal_obj = 0.5*sprod_nn(w,w,sm->sizePsi)+C*value;
     
#if (DEBUG_LEVEL>0)
    printf("ITER PRIMAL_OBJ %.4f\n", primal_obj); fflush(stdout);
#endif
    
 
    temp_var = sprod_nn(w_b,w_b,sm->sizePsi); 
    proximal_term = 0.0;
    for (i=1;i<sm->sizePsi+1;i++) {
      proximal_term += (w[i]-w_b[i])*(w[i]-w_b[i]);
    }
    
    reg_master_obj = -dual_obj+0.5*rho*temp_var/(1+rho);
    expected_descent = reg_master_obj - primal_obj_b;

    v_k = (reg_master_obj - proximal_term*rho/2) - primal_obj_b; 

    primal_lower_bound = MAX(primal_lower_bound, reg_master_obj - 0.5*rho*(1+rho)*proximal_term);

#if (DEBUG_LEVEL>0)
    printf("ITER REG_MASTER_OBJ: %.4f\n", reg_master_obj);
    printf("ITER EXPECTED_DESCENT: %.4f\n", expected_descent);
    printf("ITER PRIMLA_OBJ_B: %.4f\n", primal_obj_b);
    printf("ITER RHO: %.4f\n", rho);
    printf("ITER ||w-w_b||^2: %.4f\n", proximal_term);
    printf("ITER PRIMAL_LOWER_BOUND: %.4f\n", primal_lower_bound);
    printf("ITER V_K: %.4f\n", v_k); 
#endif
    obj_difference = primal_obj - primal_obj_b; 


    if (primal_obj<primal_obj_b+kappa*expected_descent) {
      /* extra condition to be met */
      if ((gTd>m2*v_k)||(rho<min_rho+1E-8)) {
#if (DEBUG_LEVEL>0)
	printf("SERIOUS STEP\n");
#endif
	/* update cut_error */
	for (i=0;i<size_active;i++) {
	  cut_error[i] -= (primal_obj_b - 0.5*sprod_nn(w_b,w_b,sm->sizePsi)); 
	  cut_error[i] -= C*sprod_ns(w_b, dXc[i]->fvec); 
	  cut_error[i] += (primal_obj - 0.5*sprod_nn(w,w,sm->sizePsi));
	  cut_error[i] += C*sprod_ns(w, dXc[i]->fvec); 
	}
	primal_obj_b = primal_obj;
	for (i=1;i<sm->sizePsi+1;i++) {
	  w_b[i] = w[i];
	}
	null_step = 0;
	serious_counter++;	
      } else {
	/* increase step size */
#if (DEBUG_LEVEL>0)
	printf("NULL STEP: SS(ii) FAILS.\n");
#endif
	serious_counter--; 
	rho = MAX(rho/10,min_rho);
      }
    } else { /* no sufficient decrease */
      serious_counter--; 
      if ((cut_error[size_active-1]>m3*last_sigma_k)&&(fabs(obj_difference)>last_z_k_norm+last_sigma_k)) {
#if (DEBUG_LEVEL>0)
	printf("NULL STEP: NS(ii) FAILS.\n");
#endif
	rho = MIN(10*rho,max_rho);
      } 
#if (DEBUG_LEVEL>0)
      else printf("NULL STEP\n");
#endif
    }
    /* update last_sigma_k */
    last_sigma_k = sigma_k; 
    last_z_k_norm = z_k_norm; 


    /* break away from while loop if more than certain proportioal decrease in primal objective */
    if (primal_obj_b/initial_primal_obj<1-decrease_proportion) {
      suff_decrease_cond = 1; 
    }

    /* clean up */
    if (iter % CLEANUP_CHECK == 0) {
      size_active = resize_cleanup(size_active, idle, alpha, delta, gammaG0, proximal_rhs, G, dXc, cut_error);
    }

	  
  } // end cutting plane while loop 

  printf(" Inner loop optimization finished.\n"); fflush(stdout); 
      
  /* free memory */
  for (j=0;j<size_active;j++) {
    free(G[j]);
    free_example(dXc[j],0);	
  }
  free(G);
  free(dXc);
  free(alpha);
  free(delta);
  free_svector(new_constraint);
  free(idle);
  free(gammaG0);
  free(proximal_rhs);
  free(cut_error); 

  /* copy and free */
  for (i=1;i<sm->sizePsi+1;i++) {
    w[i] = w_b[i];
  }
  free(w_b);

  return(primal_obj_b);

}



int main(int argc, char* argv[]) {

  double *w; /* weight vector */
  int outer_iter;
  long m, i;
  double C, epsilon;
  LEARN_PARM learn_parm;
  KERNEL_PARM kernel_parm;
  char trainfile[1024];
  char modelfile[1024];
  int MAX_ITER;
  /* new struct variables */
  SVECTOR **fycache, *diff, *fy;
  EXAMPLE *ex;
  SAMPLE sample;
  STRUCT_LEARN_PARM sparm;
  STRUCTMODEL sm;
  
  double decrement;
  double primal_obj, last_primal_obj;
  double cooling_eps; 
  double stop_crit; 
 

  /* read input parameters */
  my_read_input_parameters(argc, argv, trainfile, modelfile, &learn_parm, &kernel_parm, &sparm); 

  epsilon = learn_parm.eps;
  C = learn_parm.svm_c;
  MAX_ITER = learn_parm.maxiter;

  /* read in examples */
  sample = read_struct_examples(trainfile,&sparm); 
  ex = sample.examples;
  m = sample.n;
  
  /* initialization */
  init_struct_model(sample,&sm,&sparm,&learn_parm,&kernel_parm); 

  w = create_nvector(sm.sizePsi);
  clear_nvector(w, sm.sizePsi);
  sm.w = w; /* establish link to w, as long as w does not change pointer */

  /* some training information */
  printf("C: %.8g\n", C);
  printf("epsilon: %.8g\n", epsilon);
  printf("sample.n: %ld\n", sample.n); 
  printf("sm.sizePsi: %ld\n", sm.sizePsi); fflush(stdout);
  

  /* impute latent variable for first iteration */
  init_latent_variables(&sample,&learn_parm,&sm,&sparm);
  
  /* prepare feature vector cache for correct labels with imputed latent variables */
  fycache = (SVECTOR**)malloc(m*sizeof(SVECTOR*));
  for (i=0;i<m;i++) {
    fy = psi(ex[i].x, ex[i].y, ex[i].h, &sm, &sparm);
    diff = add_list_ss(fy);
    free_svector(fy);
    fy = diff;
    fycache[i] = fy;
  }
    
  /* outer loop: latent variable imputation */
  outer_iter = 0;
  last_primal_obj = 0;
  decrement = 0;
  cooling_eps = 0.5*C*epsilon; 
  while ((outer_iter<2)||((!stop_crit)&&(outer_iter<MAX_OUTER_ITER))) { 
    printf("OUTER ITER %d\n", outer_iter); 
    /* cutting plane algorithm */
    primal_obj = cutting_plane_algorithm(w, m, MAX_ITER, C, cooling_eps, fycache, ex, &sm, &sparm); 
    
    /* compute decrement in objective in this outer iteration */
    decrement = last_primal_obj - primal_obj;
    last_primal_obj = primal_obj;
    printf("primal objective: %.4f\n", primal_obj);
    printf("decrement: %.4f\n", decrement); fflush(stdout);
    
    stop_crit = (decrement<C*epsilon)&&(cooling_eps<0.5*C*epsilon+1E-8);

    cooling_eps = -decrement*0.01;
    cooling_eps = MAX(cooling_eps, 0.5*C*epsilon);
    printf("cooling_eps: %.8g\n", cooling_eps); 

  
    /* impute latent variable using updated weight vector */
    for (i=0;i<m;i++) {
      free_latent_var(ex[i].h);
      ex[i].h = infer_latent_variables(ex[i].x, ex[i].y, &sm, &sparm);
    }
    /* re-compute feature vector cache */
    for (i=0;i<m;i++) {
      free_svector(fycache[i]);
      fy = psi(ex[i].x, ex[i].y, ex[i].h, &sm, &sparm);
      diff = add_list_ss(fy);
      free_svector(fy);
      fy = diff;
      fycache[i] = fy;
    }


    outer_iter++;  
  } // end outer loop
  

  /* write structural model */
  write_struct_model(modelfile, &sm, &sparm);
  // skip testing for the moment  

  /* free memory */
  free_struct_sample(sample);
  free_struct_model(sm, &sparm);
  for(i=0;i<m;i++) {
    free_svector(fycache[i]);
  }
  free(fycache);
   
  return(0); 
  
}



void my_read_input_parameters(int argc, char *argv[], char *trainfile, char* modelfile,
			      LEARN_PARM *learn_parm, KERNEL_PARM *kernel_parm, STRUCT_LEARN_PARM *struct_parm) {
  
  long i;

  /* set default */
  learn_parm->maxiter=20000;
  learn_parm->svm_maxqpsize=100;
  learn_parm->svm_c=100.0;
  learn_parm->eps=0.001;
  learn_parm->biased_hyperplane=12345; /* store random seed */
  learn_parm->remove_inconsistent=10; 
  kernel_parm->kernel_type=0;
  kernel_parm->rbf_gamma=0.05;
  kernel_parm->coef_lin=1;
  kernel_parm->coef_const=1;
  kernel_parm->poly_degree=3;

  struct_parm->custom_argc=0;

  for(i=1;(i<argc) && ((argv[i])[0] == '-');i++) {
    switch ((argv[i])[1]) {
    case 'c': i++; learn_parm->svm_c=atof(argv[i]); break;
    case 'e': i++; learn_parm->eps=atof(argv[i]); break;
    case 's': i++; learn_parm->svm_maxqpsize=atol(argv[i]); break; 
    case 'g': i++; kernel_parm->rbf_gamma=atof(argv[i]); break;
    case 'd': i++; kernel_parm->poly_degree=atol(argv[i]); break;
    case 'r': i++; learn_parm->biased_hyperplane=atol(argv[i]); break; 
    case 't': i++; kernel_parm->kernel_type=atol(argv[i]); break;
    case 'n': i++; learn_parm->maxiter=atol(argv[i]); break;
    case 'p': i++; learn_parm->remove_inconsistent=atol(argv[i]); break; 
    case '-': strcpy(struct_parm->custom_argv[struct_parm->custom_argc++],argv[i]);i++; strcpy(struct_parm->custom_argv[struct_parm->custom_argc++],argv[i]);break; 
    default: printf("\nUnrecognized option %s!\n\n",argv[i]);
      exit(0);
    }

  }

  if(i>=argc) {
    printf("\nNot enough input parameters!\n\n");
    my_wait_any_key();
    exit(0);
  }
  strcpy (trainfile, argv[i]);

  if((i+1)<argc) {
    strcpy (modelfile, argv[i+1]);
  }
  
  parse_struct_parameters(struct_parm);

}



void my_wait_any_key()
{
  printf("\n(more)\n");
  (void)getc(stdin);
}



int resize_cleanup(int size_active, int *idle, double *alpha, double *delta, double *gammaG0, double *proximal_rhs, double **G, DOC **dXc, double *cut_error) {
  int i,j, new_size_active;
  long k;

  i=0;
  while ((i<size_active)&&(idle[i]<IDLE_ITER)) i++;
  j=i;
  while((j<size_active)&&(idle[j]>=IDLE_ITER)) j++;

  while (j<size_active) {
    /* copying */
    alpha[i] = alpha[j];
    delta[i] = delta[j];
    gammaG0[i] = gammaG0[j];
    cut_error[i] = cut_error[j]; 
    
    free(G[i]);
    G[i] = G[j]; 
    G[j] = NULL;
    free_example(dXc[i],0);
    dXc[i] = dXc[j];
    dXc[j] = NULL;

    i++;
    j++;
    while((j<size_active)&&(idle[j]>=IDLE_ITER)) j++;
  }
  for (k=i;k<size_active;k++) {
    if (G[k]!=NULL) free(G[k]);
    if (dXc[k]!=NULL) free_example(dXc[k],0);
  }
  new_size_active = i;
  alpha = (double*)realloc(alpha, sizeof(double)*new_size_active);
  delta = (double*)realloc(delta, sizeof(double)*new_size_active);
  gammaG0 = (double*)realloc(gammaG0, sizeof(double)*new_size_active);
  proximal_rhs = (double*)realloc(proximal_rhs, sizeof(double)*new_size_active);
  G = (double**)realloc(G, sizeof(double*)*new_size_active);
  dXc = (DOC**)realloc(dXc, sizeof(DOC*)*new_size_active);
  cut_error = (double*)realloc(cut_error, sizeof(double)*new_size_active); 
  
  /* resize G and idle */
  i=0;
  while ((i<size_active)&&(idle[i]<IDLE_ITER)) i++;
  j=i;
  while((j<size_active)&&(idle[j]>=IDLE_ITER)) j++;

  while (j<size_active) {
    idle[i] = idle[j];
    for (k=0;k<new_size_active;k++) {
      G[k][i] = G[k][j];
    }
    i++;
    j++;
    while((j<size_active)&&(idle[j]>=IDLE_ITER)) j++;
  }  
  idle = (int*)realloc(idle, sizeof(int)*new_size_active);
  for (k=0;k<new_size_active;k++) {
    G[k] = (double*)realloc(G[k], sizeof(double)*new_size_active);
  }
  return(new_size_active);

}
