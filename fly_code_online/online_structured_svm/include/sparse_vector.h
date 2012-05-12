#ifndef __SPARSE_VECTOR_H
#define __SPARSE_VECTOR_H

#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>



///@cond
#define VFLOAT double

#ifndef INFINITY
#define INFINITY HUGE_VAL
#endif

#define REALLOC_SCALE 1.1
#define REALLOC_ADD 128
///@endcond


/**
 * @file sparse_vector.h
 * @brief Routines for a vector that is stored sparsely in memory
 */

#include "util.h"

int int_compare (const void * a, const void * b);


/// @cond
typedef struct _SparseVectorElement {
  int ind;
  VFLOAT val;
} SparseVectorElement;
/// @endcond

/**
 * @class SparseVector
 * @brief A vector that is stored sparsely in memory
 */
class SparseVector {
  int numNonZero;
  int maxIndex;
  int numAlloc;
  SparseVectorElement *elements;
  VFLOAT *non_sparse;
  int *used_inds;
  int numUsed;
  double multiplier;
  int shiftAmount;

public:

  /**
   * @brief Create an empty vector (e.g., all entries are 0)
   */
  SparseVector() {
    numNonZero = 0;
    numAlloc = 0;
    maxIndex = -1;
    elements = NULL;
    non_sparse = NULL;
    used_inds = NULL;
    numUsed = 0;
    shiftAmount = -1;
    multiplier = 1;
  }

  /**
   * @brief Create sparse vector from a non-sparse vector
   * @param v a non-sparse vector, assumed to be an array of length n
   * @param n The length of v
   */
  SparseVector(double *v, int n, bool isNonSparse=false) {
    numNonZero = 0;
    maxIndex = n-1;
    numAlloc = 0;
    elements = NULL;
    non_sparse = NULL;
    used_inds = NULL;
    numUsed = 0;
    shiftAmount = -1;
    multiplier = 1;
    if(isNonSparse) {
      non_sparse = (double*)malloc(sizeof(double)*n);
      memcpy(non_sparse, v, sizeof(double)*n);
    } else {
      for(int i = 0; i < n; i++) {
	if(v[i]) {
	  if((numNonZero+2) > numAlloc) {
	    numAlloc = (int)(numAlloc*REALLOC_SCALE)+REALLOC_ADD;
	    elements = (SparseVectorElement*)realloc(elements, numAlloc*sizeof(SparseVectorElement));
	  }
	  elements[numNonZero].ind = i;
	  elements[numNonZero].val = v[i];
	  numNonZero++;
	}
      }
    }
  }

  
  /**
   * @brief Create a copy of another vector
   * @param v The vector to copy
   */
  SparseVector(const SparseVector &v) {
    non_sparse = NULL;
    elements = NULL;
    used_inds = NULL;
    numAlloc = numNonZero = v.numNonZero;
    maxIndex = v.maxIndex;
    numUsed = v.numUsed;
    shiftAmount = v.shiftAmount;
    multiplier = v.multiplier;
    if(v.elements) {
      if(shiftAmount >= 0) 
	elements = v.elements;
      else {
	elements = (SparseVectorElement*)malloc(sizeof(SparseVectorElement)*numAlloc);
	memcpy(elements, v.elements, sizeof(SparseVectorElement)*numAlloc); 
      }
    } else {
      assert(v.shiftAmount < 0);
      non_sparse = (double*)malloc(sizeof(double)*(v.maxIndex+1));
      memcpy(non_sparse, v.non_sparse, sizeof(double)*(v.maxIndex+1));
      if(v.used_inds) {
	used_inds = (int*)malloc(sizeof(int)*(v.maxIndex+1));
	memcpy(used_inds, v.used_inds, sizeof(int)*(v.maxIndex+1));
      }
    }
  }

  /**
   * @brief Create a copy of another vector
   * @param v The vector to copy
   * @param takeover if true, effectively destroys vector v.  This is used for converting an object v on the stack to the heap
   */
  SparseVector(SparseVector &v, bool takeover) {
    elements = NULL;
    non_sparse = NULL;
    used_inds = NULL;
    numNonZero = v.numNonZero;
    numAlloc = v.numAlloc;
    maxIndex = v.maxIndex;
    if(takeover) {
      elements = v.elements;
      v.elements = NULL;
      non_sparse = v.non_sparse;
      v.non_sparse = NULL;
      used_inds = v.used_inds;
      v.used_inds = NULL;
      numUsed = v.numUsed;
      v.numUsed = 0;
      shiftAmount = v.shiftAmount;
      v.shiftAmount = -1;
      multiplier = v.multiplier;
      v.multiplier = -1;
    } else {
      assert(v.shiftAmount < 0);
      if(v.elements) {
	elements = (SparseVectorElement*)malloc(sizeof(SparseVectorElement)*numAlloc);
	memcpy(elements, v.elements, sizeof(SparseVectorElement)*numAlloc); 
      } 
      if(v.non_sparse) {
	non_sparse = (double*)malloc(sizeof(double)*v.Length());
	memcpy(non_sparse, v.non_sparse, sizeof(double)*v.Length());
	if(v.used_inds) {
	  used_inds = (int*)malloc(sizeof(int)*(v.Length()));
	  memcpy(used_inds, v.used_inds, sizeof(int)*v.numUsed);
	}
      }
    }
  }
  
  ~SparseVector() {
    if(elements && shiftAmount < 0) free(elements);
    if(non_sparse) free(non_sparse);
    if(used_inds) free(used_inds);
  }

  SparseVectorElement *get_elements(int *num) {
    *num = numNonZero;
    return elements;
  }

  /**
   * @brief Get a non-sparse array of this sparse vector.  It should be free'd using delete []
   */
  template <class T>
  T *get_non_sparse(int sz, T *non_sparse_mem=NULL) { 
    assert(sz >= Length());
    assert(shiftAmount < 0);
    T *v = non_sparse_mem ? non_sparse_mem : (T*)malloc(sizeof(T)*sz);
    if(non_sparse) {
      int j;
      for(j = 0; j < my_min(sz,maxIndex+1); j++)
	v[j] = (T)non_sparse[j];
      for(; j < sz; j++)
	v[j] = 0;
    } else {
      if(!non_sparse_mem)
	for(int i = 0; i < sz; i++)
	  v[i] = 0;
      for(int j = 0; j < numNonZero; j++)
	v[elements[j].ind] = (T)elements[j].val;
    }
    return v;
  }
    
  void zero() {
    assert(shiftAmount < 0);
    if(non_sparse) {
      if(used_inds) {
	for(int i = 0; i < numUsed; i++) non_sparse[used_inds[i]] = 0;
      } else {
	for(int i = 0; i < maxIndex+1; i++) non_sparse[i] = 0;
      }
      numUsed = 0;
    }
    numNonZero = 0;
  }

  void make_non_sparse(bool non_s, int sz=-1, bool track_used_inds=false, double *non_sparse_mem=NULL) {
    assert(shiftAmount < 0);
    if(non_s && !non_sparse) {
      non_sparse = get_non_sparse<double>(sz == -1 ? Length() : sz, non_sparse_mem);
      numUsed = 0;
      if(track_used_inds) {
	used_inds = (int*)malloc(sizeof(int)*my_max(sz,Length()));
	for(int i = 0; i < numNonZero; i++) 
	  if(elements[i].val)
	    used_inds[numUsed++] = elements[i].ind;
      }
      maxIndex = (sz == -1 ? Length() : sz)-1;
      numAlloc = 0;
      if(elements) free(elements);
      elements = NULL;
    } else if(!non_s && non_sparse && used_inds) {
      numNonZero = 0;
      qsort(used_inds, numUsed, sizeof(int), int_compare);
      numAlloc = numUsed+1;
      elements = (SparseVectorElement*)realloc(elements, numAlloc*sizeof(SparseVectorElement));
      for(int i = 0; i < numUsed; i++) {
	if(non_sparse[used_inds[i]] && (!i || used_inds[i] != used_inds[i-1])) {
	  elements[numNonZero].ind = used_inds[i];
	  elements[numNonZero++].val = non_sparse[used_inds[i]];
	}
	if(non_sparse_mem) non_sparse_mem[used_inds[i]] = 0;
      }
      if(non_sparse != non_sparse_mem) 
	free(non_sparse);  
      non_sparse = NULL;
      free(used_inds);   used_inds = NULL;   numUsed = 0;
    } else if(!non_s && non_sparse) {
      numNonZero = 0;
      for(int i = 0; i < maxIndex+1; i++) {
	if(non_sparse[i]) {
	  if((numNonZero+2) > numAlloc) {
	    numAlloc = (int)(numAlloc*REALLOC_SCALE)+REALLOC_ADD;
	    elements = (SparseVectorElement*)realloc(elements, numAlloc*sizeof(SparseVectorElement));
	  }
	  elements[numNonZero].ind = i;
	  elements[numNonZero++].val = non_sparse[i];
	}
      }
      free(non_sparse);  non_sparse = NULL;
    }
  }


  /**
   * @brief Get the length (the 1 plus the maximum non-zero index) of this sparse vector
   */
  int Length() { return maxIndex+1; }

  /**
   * @brief Move this object from the stack to the heap
   */
  SparseVector *ptr() {  
    return new SparseVector(*this, true);
  } 

  /**
   * @brief load this vector from a string encoding str, in similar format to SVM lite files
   * @return A pointer into str after the last character parsed in str
   */
  const char *from_string(const char *str) {
    numAlloc = numNonZero = 0;
    maxIndex = -1;    
    elements = NULL;
    shiftAmount = -1;
    double w;
    int n;

    while(*str == ' ') str++;

    while(sscanf(str, "%d:%lg", &n, &w) == 2) {
      if(n < 0) {
        fprintf(stderr, "Error bad vector index %d\n", n);
        return NULL;
      }
      if((numNonZero+2) > numAlloc) {
        numAlloc = (int)(numAlloc*REALLOC_SCALE)+REALLOC_ADD;
        elements = (SparseVectorElement*)realloc(elements, numAlloc*sizeof(SparseVectorElement));
      }
      elements[numNonZero].ind = n;
      elements[numNonZero++].val = (VFLOAT)w;
      maxIndex = my_max(maxIndex, n);

      while(isdigit(*str) || *str == ':' || *str == '.' || *str == 'e' || *str == '-') str++;
      while(*str == ' ') str++;
    }
    
    return str;
  }

  /**
   * @brief Create a JSON encoding of this vector
   * @return a JSON encoding of this vector
   */
  Json::Value save() {
    assert(shiftAmount < 0);
    char *str = to_string();
    Json::Value a = str;
    free(str);
    return a;
  }

  /**
   * @brief Initialize this vector from a JSON encoding
   * @param v a JSON encoding of this vector
   */
  bool load(const Json::Value &v) {
    int len = strlen(v.asString().c_str())+1;
    char *str = new char[len];
    strcpy(str, v.asString().c_str());
    bool b = from_string(str) != NULL;
    delete [] str;
    return b;
  }

  /**
   * @brief Save this vector to a string encoding str, in similar format to SVM lite files
   * @return A pointer into str after the last character written
   */
  char *to_string() {
    assert(shiftAmount < 0);
    int alloc = 100;
    char *str = (char*)malloc(sizeof(char)*alloc);
    int len = 0;
    strcpy(str, "");
    if(non_sparse) {
      int n = 0;
      for(int i = 0; i < maxIndex+1; i++) {
	if(len+100 > alloc) { alloc = (int)(alloc*1.1 + 200); str = (char*)realloc(str, sizeof(char)*alloc); }
	if(n) {
	  strcpy(str+len, " ");
	  len++;
	}
	if(non_sparse[i]) {
	  sprintf(str+len, "%d:%lg", (int)i, non_sparse[i]);
	  len += strlen(str+len);
	  n++;
	}
      }
    } else {
      for(int i = 0; i < numNonZero; i++) {
	if(len+100 > alloc) { alloc = (int)(alloc*1.1 + 200); str = (char*)realloc(str, sizeof(char)*alloc); }
	if(i) {
	  strcpy(str+len, " ");
	  len++;
	}
	sprintf(str+len, "%d:%lg", (int)elements[i].ind, elements[i].val);
	len += strlen(str+len);
      }
    }
    return str;
  }

  /**
   * @brief load this vector from a file in binary format
   */
  bool read(FILE *fin) {
    shiftAmount = -1;
    if(!fread(&numNonZero, sizeof(int), 1, fin) ||
       !fread(&maxIndex, sizeof(int), 1, fin))
      return false;
    numAlloc = numNonZero+2;
    elements = (SparseVectorElement*)realloc(elements, numAlloc*sizeof(SparseVectorElement));
    return fread(elements, 1, sizeof(SparseVectorElement)*numNonZero, fin) ? true : false;
  }

  /**
   * @brief save this vector from a file in binary format
   */
  bool write(FILE *fout) {
    assert(shiftAmount < 0);
    make_non_sparse(false);
    return fwrite(&numNonZero, sizeof(int), 1, fout) &&
      fwrite(&maxIndex, sizeof(int), 1, fout) &&
      fwrite(elements, 1, sizeof(SparseVectorElement)*numNonZero, fout);
  }

  /**
   * @brief Set this vector equal to v
   */
  SparseVector &operator=(SparseVector const& v) {
    assert(v.shiftAmount < 0);
    numAlloc = numNonZero = v.numNonZero;
    maxIndex = v.maxIndex;
    elements = NULL;
    non_sparse = NULL;
    used_inds = NULL;
    numUsed = v.numUsed;
    if(v.elements) {
      elements = (SparseVectorElement*)malloc(sizeof(SparseVectorElement)*numAlloc);
      memcpy(elements, v.elements, sizeof(SparseVectorElement)*numAlloc);
    }
    if(v.non_sparse) {
      non_sparse = (double*)malloc(sizeof(double)*(v.maxIndex+1));
      memcpy(non_sparse, v.non_sparse, sizeof(double)*(v.maxIndex+1));
      if(v.used_inds) {
	used_inds = (int*)malloc(sizeof(int)*(v.maxIndex+1));
	memcpy(used_inds, v.used_inds, sizeof(double)*v.numUsed);
      }
    }
    return *this;
  }

  /**
   * @brief Return a copy of this vector with all indices of this vector shifted by shiftAmount
   */
  SparseVector shift(int shiftAmount) {
    assert(this->shiftAmount < 0);
    assert(!non_sparse);
    SparseVector retval;
    retval.elements = elements;
    retval.shiftAmount = shiftAmount;
    retval.numNonZero = numNonZero;
    retval.maxIndex = maxIndex;

    /*
    SparseVector retval(*this);
    for(int i = 0; i < numNonZero; i++)
      retval.elements[i].ind += shiftAmount;
    retval.maxIndex += shiftAmount;
    */
    return retval;
  }

  /**
   * @brief Return a copy of this vector with all indices of this vector shifted by shiftAmount
   */
  SparseVector extract_subset(int start, int end) {
    assert(shiftAmount < 0);
    assert(!non_sparse);
    SparseVector retval;
    for(int i = 0; i < numNonZero; i++) {
      if(elements[i].ind >= start && elements[i].ind < end) {
	if((retval.numNonZero+1) > retval.numAlloc) {
	  retval.numAlloc = (int)(retval.numAlloc*REALLOC_SCALE)+REALLOC_ADD;
	  retval.elements = (SparseVectorElement*)realloc(retval.elements, retval.numAlloc*sizeof(SparseVectorElement));
	}
	retval.maxIndex = retval.elements[retval.numNonZero].ind = elements[i].ind-start;
	retval.elements[retval.numNonZero++].val = elements[i].val;
      }
    }
    return retval;
  }

  /**
   * @brief Helper function that multiplies different sections of this vector by different scalars
   * @param region_inds a num_regions array, where the ith section is assumed to apply to vector indices region_inds[i-1]...region_inds[i]-1
   * @param region_multipliers a num_regions array, where each entry of this vector in the ith section is multiplied by the scalar  region_multipliers[i]
   * @param num_regions the number of sections
   */
  void multiply_regions(int *region_inds, VFLOAT *region_multipliers, int num_regions) {
    assert(shiftAmount < 0);
    assert(!non_sparse);
    int region = 0;
    for(int i = 0; i < numNonZero && region < num_regions; i++) {
      while(region < num_regions && elements[i].ind >= region_inds[region]) region++;
      if(region < num_regions) elements[i].val *= region_multipliers[region];
    }
  }


  /**
   * @brief Compute the dot product of this vector and v
   */
  VFLOAT dot(SparseVector const& v, bool *use=NULL) {
    VFLOAT retval = 0;
    if(non_sparse) {
      if(v.non_sparse) {
	assert(v.shiftAmount < 0 && shiftAmount < 0);
	for(int i = 0; i < (maxIndex+1) && i < (v.maxIndex+1); i++)
	  if(!use || use[i])
	    retval += non_sparse[i]*v.non_sparse[i];
      } else {
	double *ptr = v.shiftAmount > 0 ? non_sparse+v.shiftAmount : non_sparse;
	for(int i = 0; i < v.numNonZero; i++)
	  if(!use || use[elements[i].ind])
	    retval += ptr[v.elements[i].ind]*v.elements[i].val;
      }
    } else {
      assert(v.multiplier == 1 && multiplier == 1);
      if(v.shiftAmount> 0 || shiftAmount > 0) {
	int i1 = 0, i2 = 0;
	int s1 = shiftAmount >= 0 ? shiftAmount : 0;
	int s2 = v.shiftAmount >= 0 ? v.shiftAmount : 0;
	while(i1 < numNonZero && i2 < v.numNonZero) {
	  if(i1 >= numNonZero || (i2 < v.numNonZero && elements[i1].ind+s1 > v.elements[i2].ind+s2))
	    i2++;
	  else if(i2 >= v.numNonZero || (i1 < numNonZero && v.elements[i2].ind+s2 > elements[i1].ind+s1))
	    i1++;
	  else {
	    if(!use || use[elements[i1].ind+s1]) retval += elements[i1].val * v.elements[i2].val;
	    i1++; i2++;
	  }
	}
      } else {
	int i1 = 0, i2 = 0;
	while(i1 < numNonZero && i2 < v.numNonZero) {
	  if(i1 >= numNonZero || (i2 < v.numNonZero && elements[i1].ind > v.elements[i2].ind))
	    i2++;
	  else if(i2 >= v.numNonZero || (i1 < numNonZero && v.elements[i2].ind > elements[i1].ind))
	    i1++;
	  else {
	    if(!use || use[elements[i1].ind]) retval += elements[i1].val * v.elements[i2].val;
	    i1++; i2++;
	  }
	}
      }
    }
    return retval;
  }

  /**
   * @brief Return a copy of this vector scaled by f
   */
  SparseVector operator*(VFLOAT f) {
    SparseVector retval;
    if(non_sparse) {
      assert(shiftAmount < 0);
      if(used_inds) {
	qsort(used_inds, numUsed, sizeof(int), int_compare);
	retval.numAlloc = numUsed;
	retval.numNonZero = 0;
	retval.maxIndex = maxIndex;
	retval.elements = (SparseVectorElement*)malloc(sizeof(SparseVectorElement)*retval.numAlloc);
	for(int i = 0; i < numUsed; i++) {
	  if(!i || used_inds[i] != used_inds[i-1]) {
	    retval.elements[retval.numNonZero].ind = used_inds[i];
	    retval.elements[retval.numNonZero++].val = non_sparse[used_inds[i]]*f;
	  }
	}
      } else {
	retval.non_sparse = get_non_sparse<double>(Length());
	retval.maxIndex = maxIndex;
	for(int i = 0; i < maxIndex+1; i++) 
	  retval.non_sparse[i] *= f;
      }
    } else {
      retval.numAlloc = retval.numNonZero = numNonZero;
      retval.maxIndex = maxIndex;
      if(shiftAmount >= 0) {
	retval.elements = elements;
	retval.multiplier = multiplier*f;
	retval.shiftAmount = shiftAmount;
      } else {
	retval.elements = (SparseVectorElement*)malloc(sizeof(SparseVectorElement)*retval.numAlloc);
	for(int i = 0; i < numNonZero; i++) {
	  retval.elements[i].val = elements[i].val*f;
	  retval.elements[i].ind = elements[i].ind;
	}
      }
    }
    return retval;
  }

  /**
   * @brief Return a copy of this vector scaled by f
   */
  SparseVector *mult_scalar(VFLOAT f, bool *use) {
    assert(shiftAmount < 0);
    SparseVector *retval = new SparseVector(*this);
    if(non_sparse) {
      for(int i = 0; i < maxIndex+1; i++)
	if(!use || use[i]) 
	  retval->non_sparse[i] *= f;
    } else {
      for(int i = 0; i < numNonZero; i++)
	if(!use || use[retval->elements[i].ind]) 
	  retval->elements[i].val *= f;
    }
    return retval;
  }

  /**
   * @brief Return a vector that is a component-wise multiplication of this vector and v
   */
  SparseVector operator*(SparseVector const& v) {
    assert(!v.non_sparse);
    assert(shiftAmount < 0 && v.shiftAmount < 0);
    SparseVector retval;
    if(non_sparse) {
      retval.numAlloc = v.numNonZero;
      retval.elements = (SparseVectorElement*)realloc(retval.elements, retval.numNonZero*sizeof(SparseVectorElement));
      for(int i = 0; i < v.numNonZero; i++) {
	retval.elements[i].ind = v.elements[i].ind;
	retval.elements[i].val = non_sparse[v.elements[i].ind]*v.elements[i].val;
      }
    } else {
      int i1 = 0, i2 = 0;
      retval.numAlloc = my_min(numAlloc, v.numAlloc);
      retval.elements = (SparseVectorElement*)realloc(retval.elements, retval.numAlloc*sizeof(SparseVectorElement));
      while(i1 < numNonZero && i2 < v.numNonZero) {
	if(i1 >= numNonZero || (i2 < v.numNonZero && elements[i1].ind > v.elements[i2].ind)) 
	  i2++;
	else if(i2 >= v.numNonZero || (i1 < numNonZero && v.elements[i2].ind > elements[i1].ind)) 
	  i1++;
	else {
	  if((retval.numNonZero+1) > retval.numAlloc) {
	    retval.numAlloc = (int)(retval.numAlloc*REALLOC_SCALE)+REALLOC_ADD;
	    retval.elements = (SparseVectorElement*)realloc(retval.elements, retval.numAlloc*sizeof(SparseVectorElement));
	  }
	  retval.maxIndex = retval.numNonZero;
	  retval.elements[retval.numNonZero].ind = elements[i1].ind;
	  retval.elements[retval.numNonZero].val = elements[i1++].val * v.elements[i2++].val;
	  if(retval.elements[retval.numNonZero].val)
	    retval.numNonZero++;
	}
      }
    }
    return retval;
  }

  /**
   * @brief Return a vector that is the sum of this vector and v
   */
  SparseVector operator+(SparseVector const& v) {
    assert(shiftAmount < 0 && v.shiftAmount < 0);
    SparseVector retval;
    if(non_sparse) {
      retval.non_sparse = get_non_sparse<double>(my_max(Length(),v.maxIndex+1));
      retval.used_inds = (int*)malloc(sizeof(int)*Length());
      if(v.non_sparse) {
	for(int i = 0; i < my_min(v.maxIndex,maxIndex)+1; i++) {
	  if(!retval.non_sparse[i] && v.non_sparse[i]) retval.used_inds[retval.numUsed++] = i;
	  retval.non_sparse[i] -= v.non_sparse[i];
	}
      } else {
	for(int i = 0; i < v.numNonZero; i++) {
	  if(!retval.non_sparse[v.elements[i].ind]) retval.used_inds[retval.numUsed++] = v.elements[i].ind;
	  retval.non_sparse[v.elements[i].ind] += v.elements[i].val;
	}
      }
    } else {
      int i1 = 0, i2 = 0;
      retval.numAlloc = my_max(numAlloc, v.numAlloc);
      retval.elements = (SparseVectorElement*)realloc(retval.elements, retval.numAlloc*sizeof(SparseVectorElement));
      while(i1 < numNonZero || i2 < v.numNonZero) {
	if((retval.numNonZero+1) > retval.numAlloc) {
	  retval.numAlloc = (int)(retval.numAlloc*REALLOC_SCALE)+REALLOC_ADD;
	  retval.elements = (SparseVectorElement*)realloc(retval.elements, retval.numAlloc*sizeof(SparseVectorElement));
	}
	if(i1 >= numNonZero || (i2 < v.numNonZero && elements[i1].ind > v.elements[i2].ind)) {
	  retval.elements[retval.numNonZero].ind = v.elements[i2].ind;
	  retval.elements[retval.numNonZero++].val = v.elements[i2++].val;
	} else if(i2 >= v.numNonZero || (i1 < numNonZero && v.elements[i2].ind > elements[i1].ind)) {
	  retval.elements[retval.numNonZero].ind = elements[i1].ind;
	  retval.elements[retval.numNonZero++].val = elements[i1++].val;
	} else {
	  retval.elements[retval.numNonZero].ind = elements[i1].ind;
	  retval.elements[retval.numNonZero].val = elements[i1++].val + v.elements[i2++].val;
	  if(retval.elements[retval.numNonZero].val)
	    retval.numNonZero++;
	}
      }
    }
    retval.maxIndex = my_max(maxIndex, v.maxIndex);
    return retval;
  }

  /**
   * @brief Return a vector that is the difference of this vector and v (this minus v)
   */
  SparseVector operator-(SparseVector const& v) {
    SparseVector retval;
    if(non_sparse) {
      assert(shiftAmount < 0 && v.shiftAmount < 0);
      retval.non_sparse = get_non_sparse<double>(my_max(Length(),v.maxIndex+1));
      retval.used_inds = (int*)malloc(sizeof(int)*Length());
      if(v.non_sparse) {
	for(int i = 0; i < my_min(v.maxIndex,maxIndex)+1; i++) {
	  if(!retval.non_sparse[i] && v.non_sparse[i]) retval.used_inds[retval.numUsed++] = i;
	  retval.non_sparse[i] -= v.non_sparse[i];
	}
      } else {
	for(int i = 0; i < v.numNonZero; i++) {
	  if(!retval.non_sparse[v.elements[i].ind]) retval.used_inds[retval.numUsed++] = v.elements[i].ind;
	  retval.non_sparse[v.elements[i].ind] -= v.elements[i].val;
	}
      }
    } else {
      assert(v.multiplier == 1 && multiplier == 1);
      if(v.shiftAmount> 0 || shiftAmount > 0) {
	int i1 = 0, i2 = 0;
	int s1 = shiftAmount >= 0 ? shiftAmount : 0;
	int s2 = v.shiftAmount >= 0 ? v.shiftAmount : 0;
	retval.numAlloc = my_max(numAlloc, v.numAlloc);
	retval.elements = (SparseVectorElement*)realloc(retval.elements, retval.numAlloc*sizeof(SparseVectorElement));
	while(i1 < numNonZero || i2 < v.numNonZero) {
	  if((retval.numNonZero+1) > retval.numAlloc) {
	    retval.numAlloc = (int)(retval.numAlloc*REALLOC_SCALE)+REALLOC_ADD;
	    retval.elements = (SparseVectorElement*)realloc(retval.elements, retval.numAlloc*sizeof(SparseVectorElement));
	  }
	  if(i1 >= numNonZero || (i2 < v.numNonZero && elements[i1].ind+s1 > v.elements[i2].ind+s2)) {
	    retval.elements[retval.numNonZero].ind = v.elements[i2].ind+s2;
	    retval.elements[retval.numNonZero++].val = -v.elements[i2++].val;
	  } else if(i2 >= v.numNonZero || (i1 < numNonZero && v.elements[i2].ind+s2 > elements[i1].ind+s1)) {
	    retval.elements[retval.numNonZero].ind = elements[i1].ind+s1;
	    retval.elements[retval.numNonZero++].val = elements[i1++].val;
	  } else {
	    retval.elements[retval.numNonZero].ind = elements[i1].ind+s1;
	    retval.elements[retval.numNonZero].val = elements[i1++].val - v.elements[i2++].val;
	    if(retval.elements[retval.numNonZero].val)
	      retval.numNonZero++;
	  }
	}
      } else {
	int i1 = 0, i2 = 0;
	retval.numAlloc = my_max(numAlloc, v.numAlloc);
	retval.elements = (SparseVectorElement*)realloc(retval.elements, retval.numAlloc*sizeof(SparseVectorElement));
	while(i1 < numNonZero || i2 < v.numNonZero) {
	  if((retval.numNonZero+1) > retval.numAlloc) {
	    retval.numAlloc = (int)(retval.numAlloc*REALLOC_SCALE)+REALLOC_ADD;
	    retval.elements = (SparseVectorElement*)realloc(retval.elements, retval.numAlloc*sizeof(SparseVectorElement));
	  }
	  if(i1 >= numNonZero || (i2 < v.numNonZero && elements[i1].ind > v.elements[i2].ind)) {
	    retval.elements[retval.numNonZero].ind = v.elements[i2].ind;
	    retval.elements[retval.numNonZero++].val = -v.elements[i2++].val;
	  } else if(i2 >= v.numNonZero || (i1 < numNonZero && v.elements[i2].ind > elements[i1].ind)) {
	    retval.elements[retval.numNonZero].ind = elements[i1].ind;
	    retval.elements[retval.numNonZero++].val = elements[i1++].val;
	  } else {
	    retval.elements[retval.numNonZero].ind = elements[i1].ind;
	    retval.elements[retval.numNonZero].val = elements[i1++].val - v.elements[i2++].val;
	    if(retval.elements[retval.numNonZero].val)
	      retval.numNonZero++;
	  }
	}
      }
    }
    retval.maxIndex = my_max(maxIndex, v.maxIndex);
    return retval;
  }

  
  /**
   * @brief Multiply this vector times f
   */
  void operator*=(VFLOAT f) {
    if(non_sparse) {
      if(used_inds) {
	for(int i = 0; i < numUsed; i++) 
	  non_sparse[used_inds[i]] *= f;
      } else {
	for(int i = 0; i < maxIndex+1; i++) 
	  non_sparse[i] *= f;
      }
    } else {
      for(int i = 0; i < numNonZero; i++)
	elements[i].val *= f;
    }
  }

  /**
   * @brief Do a component-wise multiply by vector times v
   */
  void operator*=(SparseVector const& v) {
    assert(shiftAmount < 0 && v.shiftAmount < 0);
    SparseVector retval = (*this) * v;
    takeover(retval);
  } 

  /**
   * @brief Add v to this vector
   */
  void operator+=(SparseVector const& v) {
    if(non_sparse && !v.non_sparse) {
      if(used_inds) {
	for(int i = 0; i < v.numNonZero; i++) {
	  if(!non_sparse[v.elements[i].ind] && v.elements[i].val) used_inds[numUsed++] = v.elements[i].ind;
	  non_sparse[v.elements[i].ind] += v.elements[i].val;
	}
      } else {
	double *ptr = v.shiftAmount > 0 ? non_sparse + v.shiftAmount : non_sparse;
	if(v.multiplier != 1) {
	  double *ptr = non_sparse + v.shiftAmount;
	  for(int i = 0; i < v.numNonZero; i++) 
	    ptr[v.elements[i].ind] += v.multiplier*v.elements[i].val;
	} else {
	  for(int i = 0; i < v.numNonZero; i++) 
	    ptr[v.elements[i].ind] += v.elements[i].val;
	}
      }
    } else {
      assert(shiftAmount < 0 && v.shiftAmount < 0);
      SparseVector retval = (*this) + v;
      takeover(retval);
    }
  }

  /**
   * @brief subtract v from this vector
   */
  void operator-=(SparseVector const& v) {
    if(non_sparse && !v.non_sparse) {
      if(used_inds) {
	for(int i = 0; i < v.numNonZero; i++) {
	  if(!non_sparse[v.elements[i].ind] && v.elements[i].val) used_inds[numUsed++] = v.elements[i].ind;
	  non_sparse[v.elements[i].ind] -= v.elements[i].val;
	}
      } else {
	double *ptr = v.shiftAmount > 0 ? non_sparse + v.shiftAmount : non_sparse;
	if(v.multiplier != 1) {
	  double *ptr = non_sparse + v.shiftAmount;
	  for(int i = 0; i < v.numNonZero; i++) 
	    ptr[v.elements[i].ind] -= v.multiplier*v.elements[i].val;
	} else {
	  for(int i = 0; i < v.numNonZero; i++) 
	    ptr[v.elements[i].ind] -= v.elements[i].val;
	}
      }
    } else {
      SparseVector retval = (*this) - v;
      takeover(retval);
    }
  }

private:
  /*
   * @brief helper function to takeover memory allocated in v
   */
  void takeover(SparseVector &v) {
    if(elements && shiftAmount < 0) free(elements);
    if(non_sparse) free(non_sparse);
    if(used_inds) free(used_inds);
    non_sparse = v.non_sparse;
    elements = v.elements;
    used_inds = v.used_inds;
    maxIndex = v.maxIndex;
    numNonZero = v.numNonZero;
    numAlloc = v.numAlloc;
    numUsed = v.numUsed;
    shiftAmount = v.shiftAmount;
    multiplier = v.multiplier;
    v.non_sparse = NULL;
    v.elements = NULL;
    v.numAlloc = 0;
    v.numUsed = 0;
    v.multiplier = 1;
    v.shiftAmount = -11;
  }
};




#endif

