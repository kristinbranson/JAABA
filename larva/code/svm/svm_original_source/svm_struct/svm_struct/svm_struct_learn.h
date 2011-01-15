/***********************************************************************/
/*                                                                     */
/*   svm_struct_learn.h                                                */
/*                                                                     */
/*   Basic algorithm for learning structured outputs (e.g. parses,     */
/*   sequences, multi-label classification) with a Support Vector      */ 
/*   Machine.                                                          */
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

#ifndef SVM_STRUCT_LEARN
#define SVM_STRUCT_LEARN

#ifdef __cplusplus
extern "C" {
#endif
#include "../svm_light/svm_common.h"
#include "../svm_light/svm_learn.h"
#ifdef __cplusplus
}
#endif
#include "svm_struct_common.h" 
#include "../svm_struct_api_types.h" 

#define  SLACK_RESCALING    1
#define  MARGIN_RESCALING   2

#define  NSLACK_ALG               0
#define  NSLACK_SHRINK_ALG        1
#define  ONESLACK_PRIMAL_ALG      2
#define  ONESLACK_DUAL_ALG        3
#define  ONESLACK_DUAL_CACHE_ALG  4

typedef struct ccacheelem {
  SVECTOR *fydelta; /* left hand side of constraint */
  double  rhs;      /* right hand side of constraint */
  double  viol;     /* violation score under current model */
  struct ccacheelem *next; /* next in linked list */
} CCACHEELEM;

typedef struct ccache {
  int        n;              /* number of examples */
  CCACHEELEM **constlist;    /* array of pointers to constraint lists
				- one list per example. The first
				element of the list always points to
				the most violated constraint under the
				current model for each example. */
  STRUCTMODEL *sm;           /* pointer to model */
  double  *avg_viol_gain; /* array of average values by which
			     violation of globally most violated
			     constraint exceeds that of most violated
			     constraint in cache */
  int     *changed;       /* array of boolean indicating whether the
			     most violated ybar change compared to
			     last iter? */
} CCACHE;

void find_most_violated_constraint(SVECTOR **fydelta, double *lossval, 
				   EXAMPLE *ex, SVECTOR *fycached, long n, 
				   STRUCTMODEL *sm,STRUCT_LEARN_PARM *sparm,
				   double *rt_viol, double *rt_psi, 
				   long *argmax_count);
CCACHE *create_constraint_cache(SAMPLE sample, STRUCT_LEARN_PARM *sparm, 
				STRUCTMODEL *sm);
void free_constraint_cache(CCACHE *ccache);
double add_constraint_to_constraint_cache(CCACHE *ccache, MODEL *svmModel, 
	  				  int exnum, SVECTOR *fydelta, 
					  double rhs, double gainthresh,
					  int maxconst, double *rt_cachesum);
void update_constraint_cache_for_model(CCACHE *ccache, MODEL *svmModel);
double compute_violation_of_constraint_in_cache(CCACHE *ccache, double thresh);
double find_most_violated_joint_constraint_in_cache(CCACHE *ccache, 
  		     double thresh, double *lhs_n, SVECTOR **lhs, double *rhs);
void svm_learn_struct(SAMPLE sample, STRUCT_LEARN_PARM *sparm,
		      LEARN_PARM *lparm, KERNEL_PARM *kparm, 
		      STRUCTMODEL *sm, int alg_type);
void svm_learn_struct_joint(SAMPLE sample, STRUCT_LEARN_PARM *sparm,
		      LEARN_PARM *lparm, KERNEL_PARM *kparm, 
		      STRUCTMODEL *sm, int alg_type);
void svm_learn_struct_joint_custom(SAMPLE sample, STRUCT_LEARN_PARM *sparm,
		      LEARN_PARM *lparm, KERNEL_PARM *kparm, 
		      STRUCTMODEL *sm);
void remove_inactive_constraints(CONSTSET *cset, double *alpha, 
			         long i, long *alphahist, long mininactive);
MATRIX *init_kernel_matrix(CONSTSET *cset, KERNEL_PARM *kparm); 
MATRIX *update_kernel_matrix(MATRIX *matrix, int newpos, CONSTSET *cset,
			     KERNEL_PARM *kparm);
 
#endif


