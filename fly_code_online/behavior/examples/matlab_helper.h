#ifndef MATLAB_HELPER_H
#define MATLAB_HELPER_H

#include "mat.h"
#include "matrix.h"


#define MAT_GET_VARIABLE(var, name) \
  var = matGetVariable(pmat, name); \
  if(!var) { \
    fprintf(stderr, "Error finding variable '%s' in %s\n", name, pname); \
    return false; \
  }
#define MAT_GET_DOUBLE_FIELD(arr, var, name) \
  tmp = mxGetField(arr, id, name); \
  if(!tmp) { \
    fprintf(stderr, "Error finding field '%s' in %s\n", name, pname); \
    return false; \
  } \
  var = *(double*)mxGetData(tmp); \
  /*mxDestroyArray(tmp);*/
#define MAT_GET_STRING(var, name) \
  if(!tmp) { \
    fprintf(stderr, "Error finding field '%s' in %s\n", name, pname); \
    return false; \
  } \
  mxGetString(tmp, var, mxGetM(tmp)*mxGetN(tmp)+1); \
  /*mxDestroyArray(tmp);*/
#define MAT_GET_STRING_FIELD(arr, var, name) \
  tmp = mxGetField(arr, id, name); \
  MAT_GET_STRING(var, name);
#define MAT_GET_CELL_STRING(arr, var, i, name) \
  tmp = mxGetCell(arr, i); \
  MAT_GET_STRING(var, name)
#define MAT_GET_DOUBLE_ARRAY(var, name, len)	\
  if(!tmp) { \
    fprintf(stderr, "Error finding field '%s' in %s\n", name, pname); \
    return false; \
  } \
  if(len && (int)(mxGetM(tmp)*mxGetN(tmp)) != len) {			\
    fprintf(stderr, "Error reading %s, data length was %d, expected %d\n", name, (int)(mxGetM(tmp)*mxGetN(tmp)), len); \
    return false; \
  } else if(!len) {		   \
    len = mxGetM(tmp)*mxGetN(tmp); \
    var = (double*)malloc(sizeof(double)*len); \
  } \
  memcpy(var, mxGetData(tmp), len*sizeof(double));	\
  /*mxDestroyArray(tmp);*/
#define MAT_GET_DOUBLE_ARRAY_FIELD(arr, var, name, len)	\
  tmp = mxGetField(arr, id, name); \
  MAT_GET_DOUBLE_ARRAY(var, name, len);
#define MAT_GET_CELL_DOUBLE_ARRAY_VARIABLE(var, name, len, id)	\
  tmp2 = matGetVariable(pmat, name); \
  if(tmp2) tmp = mxGetCell(tmp2, id);	\
  MAT_GET_DOUBLE_ARRAY(var, name, len); \
  mxDestroyArray(tmp2);
#define MAT_READ_DATA_FILE(var, fname, len)		\
  pmat = matOpen(fname, "r"); \
  if (pmat == NULL) { \
    fprintf(stderr, "Error opening mat file %s\n", fname); \
    return false; \
  } \
  MAT_GET_VARIABLE(tmp, "data"); \
  tmp2 = mxGetFieldByNumber(tmp, id, 0); \
  if(!tmp2) { \
    fprintf(stderr, "Error finding field 'data' in %s\n", pname); \
    return false; \
  } \
  if((int)(mxGetM(tmp2)*mxGetN(tmp2)) != len) {				\
    fprintf(stderr, "Error in %s, data length was %d, expected %d\n", fname, (int)(mxGetM(tmp2)*mxGetN(tmp2)), len); \
    return false; \
  } \
  memcpy(var, mxGetData(tmp2), sizeof(double)*len); \
  matClose(pmat); \
  mxDestroyArray(tmp);		 \
  /*mxDestroyArray(tmp2);*/

#endif
