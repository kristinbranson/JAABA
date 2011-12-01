#ifndef __SPARSE_VECTOR_H
#define __SPARSE_VECTOR_H

#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#define VFLOAT double

#ifndef INFINITY
#define INFINITY HUGE_VAL
#endif

#define REALLOC_SCALE 1.1
#define REALLOC_ADD 128


#include "util.h"


typedef struct _SparseVectorElement {
  int ind;
  VFLOAT val;
} SparseVectorElement;

/**
 * @class SparseVector
 * @brief A vector that is stored sparsely in memory
 */
class SparseVector {
  int numNonZero;
  int maxIndex;
  int numAlloc;
  SparseVectorElement *elements;

public:

  /**
   * @brief Create an empty vector (e.g., all entries are 0)
   */
  SparseVector() {
    numNonZero = 0;
    numAlloc = 0;
    maxIndex = -1;
    elements = NULL;
  }

  /**
   * @brief Create sparse vector from a non-sparse vector
   * @param v a non-sparse vector, assumed to be an array of length n
   * @param n The length of v
   */
  SparseVector(double *v, int n) {
    numNonZero = 0;
    maxIndex = n-1;
    numAlloc = 0;
    elements = NULL;
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

  
  /**
   * @brief Create a copy of another vector
   * @param v The vector to copy
   */
  SparseVector(const SparseVector &v) {
    numAlloc = numNonZero = v.numNonZero;
    maxIndex = v.maxIndex;
    elements = (SparseVectorElement*)malloc(sizeof(SparseVectorElement)*numAlloc);
    memcpy(elements, v.elements, sizeof(SparseVectorElement)*numAlloc); 
  }

  /**
   * @brief Create a copy of another vector
   * @param v The vector to copy
   * @param takeover if true, effectively destroys vector v.  This is used for converting an object v on the stack to the heap
   */
  SparseVector(SparseVector &v, bool takeover) {
    numNonZero = v.numNonZero;
    numAlloc = v.numAlloc;
    maxIndex = v.maxIndex;
    if(takeover) {
      elements = v.elements;
      v.elements = NULL;
    } else {
      elements = (SparseVectorElement*)malloc(sizeof(SparseVectorElement)*numAlloc);
      memcpy(elements, v.elements, sizeof(SparseVectorElement)*numAlloc); 
    }
  }
  
  ~SparseVector() {
    if(elements) free(elements);
  }

  /**
   * @brief Get a non-sparse array of this sparse vector.  It should be free'd using delete []
   */
  template <class T>
  T *non_sparse(int sz) { 
    assert(sz >= Length());
    T *v = new T[sz];
    for(int i = 0; i < sz; i++)
      v[i] = 0;
    for(int j = 0; j < numNonZero; j++)
      v[elements[j].ind] = (T)elements[j].val;
    return v;
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

  Json::Value save() {
    char *str = new char[3000000];
    to_string(str);
    Json::Value a = str;
    delete [] str;
    return a;
    /*
    Json::Value a;
    for(int i = 0; i < numNonZero; i++) {
      Json::Value v;  v["i"] = elements[i].ind;  v["v"] = elements[i].val;
      a[i] = v;
    }
    return a;*/
  }
  bool load(const Json::Value &v) {
    char *str = new char[3000000];
    strcpy(str, v.asString().c_str());
    bool b = from_string(str) != NULL;
    delete [] str;
    return b;
    /*
    maxIndex = -1;    
    if(!v.isArray()) { fprintf(stderr, "Error reading sparse vector\n"); return false; }
    numAlloc = numNonZero = (int)v.size();
    elements = (SparseVectorElement*)realloc(elements, numAlloc*sizeof(SparseVectorElement));
    for(int i = 0; i < numNonZero; i++) {
      if((elements[i].ind=v[i].get("i",-1).asInt()) < 0 || !v[i].isMember("v")) {
	free(elements); elements = NULL;
	numAlloc = numNonZero = 0;
	fprintf(stderr, "Error reading sparse vector\n"); return false; 
      } else
	elements[i].val = v[i]["v"].asDouble();
      maxIndex = my_max(maxIndex, elements[i].ind);
    }
    return true;
    */
  }

  /**
   * @brief Save this vector to a string encoding str, in similar format to SVM lite files
   * @return A pointer into str after the last character written
   */
  void to_string(char *str) {
    strcpy(str, "");
    for(int i = 0; i < numNonZero; i++) {
      if(i) {
        strcpy(str, " ");
        str++;
      }
      sprintf(str, "%d:%lg", (int)elements[i].ind, elements[i].val);
      str += strlen(str);
    }
  }

  /**
   * @brief load this vector from a file in binary format
   */
  bool read(FILE *fin) {
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
    return fwrite(&numNonZero, sizeof(int), 1, fout) &&
      fwrite(&maxIndex, sizeof(int), 1, fout) &&
      fwrite(elements, 1, sizeof(SparseVectorElement)*numNonZero, fout);
  }

  /**
   * @brief Set this vector equal to v
   */
  SparseVector &operator=(SparseVector const& v) {
    numAlloc = numNonZero = v.numNonZero;
    maxIndex = v.maxIndex;
    elements = (SparseVectorElement*)malloc(sizeof(SparseVectorElement)*numAlloc);
    memcpy(elements, v.elements, sizeof(SparseVectorElement)*numAlloc);
    return *this;
  }

  /**
   * @brief Return a copy of this vector with all indices of this vector shifted by shiftAmount
   */
  SparseVector shift(int shiftAmount) {
    SparseVector retval(*this);
    for(int i = 0; i < numNonZero; i++)
      retval.elements[i].ind += shiftAmount;
    retval.maxIndex += shiftAmount;
    return retval;
  }

  /**
   * @brief Helper function that multiplies different sections of this vector by different scalars
   * @param region_inds a num_regions array, where the ith section is assumed to apply to vector indices region_inds[i-1]...region_inds[i]-1
   * @param region_multipliers a num_regions array, where each entry of this vector in the ith section is multiplied by the scalar  region_multipliers[i]
   * @param num_regions the number of sections
   */
  void multiply_regions(int *region_inds, VFLOAT *region_multipliers, int num_regions) {
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
    int i1 = 0, i2 = 0;
    while(i1 < numNonZero || i2 < v.numNonZero) {
      if(i1 >= numNonZero || (i2 < v.numNonZero && elements[i1].ind > v.elements[i2].ind))
        i2++;
      else if(i2 >= v.numNonZero || (i1 < numNonZero && v.elements[i2].ind > elements[i1].ind))
        i1++;
      else {
        if(!use || use[elements[i1].ind]) retval += elements[i1].val * v.elements[i2].val;
        i1++; i2++;
      }
    }
    return retval;
  }

  /**
   * @brief Return a copy of this vector scaled by f
   */
  SparseVector operator*(VFLOAT f) {
    SparseVector retval = *this;
    for(int i = 0; i < numNonZero; i++)
      retval.elements[i].val *= f;
    return retval;
  }

  /**
   * @brief Return a copy of this vector scaled by f
   */
  SparseVector *mult_scalar(VFLOAT f, bool *use) {
    SparseVector *retval = new SparseVector(*this);
    for(int i = 0; i < numNonZero; i++)
      if(!use || use[retval->elements[i].ind]) 
	retval->elements[i].val *= f;
    return retval;
  }

  /**
   * @brief Return a vector that is a component-wise multiplication of this vector and v
   */
  SparseVector operator*(SparseVector const& v) {
    SparseVector retval;
    int i1 = 0, i2 = 0;
	retval.numAlloc = my_min(numAlloc, v.numAlloc);
    retval.elements = (SparseVectorElement*)realloc(retval.elements, retval.numAlloc*sizeof(SparseVectorElement));
    while(i1 < numNonZero || i2 < v.numNonZero) {
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
    return retval;
  }

  /**
   * @brief Return a vector that is the sum of this vector and v
   */
  SparseVector operator+(SparseVector const& v) {
    SparseVector retval;
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
    retval.maxIndex = my_max(maxIndex, v.maxIndex);
    return retval;
  }

  /**
   * @brief Return a vector that is the difference of this vector and v (this minus v)
   */
  SparseVector operator-(SparseVector const& v) {
    SparseVector retval;
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
    retval.maxIndex = my_max(maxIndex, v.maxIndex);
    return retval;
  }

  
  /**
   * @brief Multiply this vector times f
   */
  void operator*=(VFLOAT f) {
    for(int i = 0; i < numNonZero; i++)
      elements[i].val *= f;
  }

  /**
   * @brief Do a component-wise multiply by vector times v
   */
  void operator*=(SparseVector const& v) {
    SparseVector retval = (*this) * v;
    takeover(retval);
  } 

  /**
   * @brief Add v to this vector
   */
  void operator+=(SparseVector const& v) {
    SparseVector retval = (*this) + v;
    takeover(retval);
  }

  /**
   * @brief subtract v from this vector
   */
  void operator-=(SparseVector const& v) {
    SparseVector retval = (*this) - v;
    takeover(retval);
  }

private:
  /*
   * @brief helper function to takeover memory allocated in v
   */
  void takeover(SparseVector &v) {
    if(elements) free(elements);
    elements = v.elements;
    maxIndex = v.maxIndex;
    numNonZero = v.numNonZero;
    numAlloc = v.numAlloc;
    v.elements = NULL;
    v.numAlloc = 0;
  }
};




#endif

