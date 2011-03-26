/***********************************************************************/
/*                                                                     */
/*   svm_struct_main.c                                                 */
/*                                                                     */
/*   Command line interface to the alignment learning module of the    */
/*   Support Vector Machine.                                           */
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



/* the following enables you to use svm-learn out of C++ */

#include "../svm_light/svm_common.h"
#include "../svm_light/svm_learn.h"

# include "svm_struct_learn.h"
# include "svm_struct_common.h"

#include <stdio.h>
#include <string.h>
#include <assert.h>

char trainfile[200];           /* file with training examples */
char modelfile[200];           /* file for resulting classifier */
 

int main (int argc, const char* argv[])
{  
  STRUCT_LEARN_PARM struct_parm;
  STRUCTMODEL structmodel;

  memset(&structmodel, 0, sizeof(STRUCTMODEL));
  memset(&struct_parm, 0, sizeof(STRUCT_LEARN_PARM));
  train_main(argc, argv, &struct_parm, &structmodel, NULL);

  return 0;
}

/*---------------------------------------------------------------------------*/


