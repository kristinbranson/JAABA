/***********************************************************************/
/*                                                                     */
/*   svm_struct_classify.c                                             */
/*                                                                     */
/*   Classification module of SVM-struct.                              */
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
/************************************************************************/

#include <stdio.h>

#include "../svm_light/svm_common.h"

#include "../svm_struct_api.h"
#include "svm_struct_common.h"
#include "../svm_struct_api_multiclass.h"



int main(int argc, const char* argv[])
{
  STRUCT_LEARN_PARM struct_parm;

  memset(&struct_parm, 0, sizeof(STRUCT_LEARN_PARM));
  return test_main(argc, argv, &struct_parm, NULL, NULL);
}
  
