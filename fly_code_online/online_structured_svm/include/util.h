#ifndef __UTIL_H
#define __UTIL_H

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <vector>
#include <map>


#include "json/json.h"

/**
 * @file util.h
 * @brief Miscellaneous helper routines 
 */

//#define DEBUG_MEMORY_LEAKS

#ifndef WIN32
#undef DEBUG_MEMORY_LEAKS
#endif

#ifdef DEBUG_MEMORY_LEAKS
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#endif

#ifndef M_E
#define M_E 2.7182818284590452354
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef INFINITY
#define INFINITY HUGE_VAL
#endif


#define my_max(x,y) ((x) > (y) ? (x) : (y))
#define my_min(x,y) ((x) < (y) ? (x) : (y))
#define my_abs(x) ((x) < 0 ? -(x) : (x))
#define my_round(x) (x < 0 ? ((int)(x-.5f)) : ((int)(x+.5f)))

#define SAFE(x, mi, ma) (my_min(ma,my_max(mi,x)))
#define SQR(x) ((x)*(x))
#define LOG2(x) (log((double)x)/log(2.0))
#define LOG_E(x) (log((double)x)/log(M_E))
#define LOG_B(x,b) (log((double)x)/log((double)b))

#define IN_RANGE(v, lower, upper) my_min(upper, my_max(v, lower))

#define RAND_FLOAT (((float)rand())/RAND_MAX)
#define RAND_DOUBLE (((double)rand())/RAND_MAX)
#define RAND_GAUSSIAN (sqrt(-2*LOG_E(RAND_DOUBLE))*cos(2*M_PI*RAND_DOUBLE))

#define SIGMOID2(t) (1.0f/(1+exp(-t)))

#ifndef isnan
#define isnan(x) ((x) != (x))
#endif

#ifdef WIN32
#define SLEEP(us) Sleep(us/1000)
#else
#define SLEEP(us) usleep(us)
#endif

#ifndef WIN32
  #define stricmp strcasecmp
  #define strnicmp strncasecmp
#endif


#ifdef WIN32
#include <winsock2.h>
#include <windows.h>
#include "jsonrpc/stdint.h"
#else
#include <stdint.h>
#endif

inline void chomp(char *str) { 
	int i = (int)strlen(str)-1; 
	while(i >= 0 && isspace(str[i])) {
		str[i--] = '\0';
	}
}

inline void pre_chomp(char *str) { 
  int i = 0, n = 0;
  while(isspace(str[i])) i++;
  if(i) { 
    while(str[i]) {
     str[n++] = str[i];
     i++;
    }
    str[n] = '\0';
  }
}


inline char *my_strtok(char *str, const char *delim, char **next) {
  char *str2 = str ? str : *next;
  if(!str2) return NULL;
  int len = (int)strlen(str2);
  for(int i = 0; i < len; i++) {
    for(int j = 0; j < (int)strlen(delim); j++) {
      if(str2[i] == delim[j]) {
        str2[i] = '\0';
        *next = i < len-1 ? str2 + i + 1 : NULL;
        return str2;
      }
    }
  }
  *next = NULL;
  return str2;
}


inline void tolower_str(char *str) { 
  for(int i = 0; i < (int)strlen(str); i++)
    str[i] = tolower(str[i]);
}

inline void toupper_str(char *str) { 
  for(int i = 0; i < (int)strlen(str); i++)
    str[i] = toupper(str[i]);
}

inline char *StringCopy(const char *str) { char *retval = (char*)malloc((strlen(str)+1)*sizeof(char)); strcpy(retval, str); return retval; }
inline void StringFree(char *str) { free(str); }

inline void StringToUpperCase(char *str) {
  for(int i = 0; i < (int)strlen(str); i++)
    str[i] = toupper(str[i]);
}

inline void StringToLowerCase(char *str) {
  for(int i = 0; i < (int)strlen(str); i++)
    str[i] = tolower(str[i]);
}

// Given a file path f1 from the current directory, get the relative path to 
// f1 if we are in the folder f2 instead of in the current directory
inline char *GetRelativePath(const char *f1, const char *f2, char *out) {
  int i;
  strcpy(out, "");

  if(f1[0] == '/' || f1[0] == '\\') {
    strcpy(out, f1); 
    return out; 
  }

  char rels[1000]; strcpy(rels, "");
  i = 0;
  bool isEq = true;
  int lastEq = -1, lastDirEq=0;
  while(f2[i]) {
    if(isEq && f1[i] != f2[i]) {
      isEq = false;
      lastEq = i-1;
    }
    if(f2[i] == '/') {
      if(!isEq) strcat(rels, "../");
      else lastDirEq = i;
    }
    i++;
  }

  while(f1[lastDirEq] == '/')
    lastDirEq++;
  strcat(out, rels);
  strcat(out, f1+lastDirEq);
  return out;
}

inline char *GetFullPath(const char *relativePath, const char *currDir, char *fullPath) {
  if(relativePath[0] == '/' || (relativePath[1] == ':' && 
				(relativePath[2] == '/' || relativePath[2] == '\\'))) {
    strcpy(fullPath, relativePath);   // path is already a full path
  } else {
    sprintf(fullPath, "%s/%s", currDir, relativePath);
  }
  return fullPath;
}


inline void ExtractFilename(const char *str, char *fname) {
  int i = (int)strlen(str)-1;
  while(i >= 0 && str[i] != '/' && str[i] != '\\')
    i--;
  strcpy(fname, str+i+1);
}


inline void ExtractFolderAndFileName(const char *str, char *folder, char *fname) {
  int i = (int)strlen(str)-1;
  while(i >= 0 && str[i] != '/' && str[i] != '\\')
    i--;
  strcpy(fname, str+i+1);
  
  if(i > 0) i--;
  int num = 0;
  while(i >= 0 && str[i] != '/' && str[i] != '\\') {
    i--;
    num++;
  }
  strcpy(folder, str+i+1);
  folder[num] = '\0';
}

inline void ExtractPathname(const char *str, char *fname) {
  int i = (int)strlen(str);
  while(i > 0 && str[i] != '/' && str[i] != '\\') {
    i--;
  }
  strcpy(fname, str);
  fname[i] = '\0';
}

inline void StripFileExtension(char *str) {
  int i = (int)strlen(str);
  while(i >= 0) {
    if(str[i] == '.') {
      strcpy(str+i, "");
      break;
    }
    i--;
  }
}

inline const char *GetFileExtension(const char *str) {
  int i = (int)strlen(str);
  while(i >= 0) {
    if(str[i] == '.') {
      return str+i+1;
      break;
    }
    i--;
  }
  return NULL;
}

inline void StringReplaceChar(char *str, char c, char c2) {
  for(int i = 0; i < (int)strlen(str); i++)
    if(str[i] == c)
      str[i] = c2;
}

inline int SplitString(char *str, char **toks, const char *d) {
  int num = 0;
  char *ptr;
  while((ptr=strtok(str, d)) != NULL) {
    toks[num++] = ptr;
    str = NULL;
  }
  return num;
}

inline bool IsSafeFileName(const char *fname) { 
  return !strstr(fname, "/") && !strstr(fname, "\\") && !strstr(fname, "~"); 
}

inline bool FileExists(const char *fname) {
  FILE *fin = fopen(fname, "r");
  if(fin) fclose(fin);
  return fin != NULL;
}

template <typename T>
inline T L2Norm(T *v1, T *v2, int sz) {
  T sum = 0, d;

  for(int i = 0; i < sz; i++) {
    d = v1[i]-v2[i];
    sum += SQR(d);
  }
  return sqrt(sum);
}

inline float mod_2pi(float x, float d=2*M_PI)
{
  while(x < 0) x += d;
  while(x > (float)(2*M_PI)) x -= d;
  return x;
}


inline char *sgets(char *buf, int n, char **source) {
  int k = 0;
  //fprintf(stderr, "sgets %s\n", *source);

  while (n-- > 1 && **source && **source != '\n')  {
    buf[k++] = *(*source)++;
  }
  buf[k] = 0;
  if (n < 1)
    return buf;
  if (!**source) {
    if (k <= 0)
      return NULL;
    buf[k++] = '\n';
    buf[k] = 0;
    //(*source)++;
    return buf;
  }
  buf[k++] = *(*source)++;
  buf[k] = 0;
  return buf;
}


inline char *LoadTextFile(const char *fname) {
  FILE *fin = fopen(fname, "r");
  if(!fin) return NULL;

  char line[1000];
  char *retval = (char*)malloc(4);
  int len = 50000; 
  *(int*)retval = 0;
  int num = 0;
  while(fgets(line, 999, fin)) {
    chomp(line);
    retval = (char*)realloc(retval, len+strlen(line)+100);
    ((int*)retval)[1+num] = len;
    strcpy(retval+len, line);
    len += (int)strlen(line)+1;
    *(int*)retval = num++;
  }
  
  fclose(fin);
  return retval;
}

inline bool HasString(const char *files, const char *str) {
  int num = *((int*)files);
  for(int i = 0; i < num; i++) {
    if(!strcmp(files + ((int*)files)[i+1], str))
      return true;
  }
  return false;
}

inline char *ReadStringFile(const char *fname) {
  int len = 0;
  char *retval = NULL;
  FILE *fin = fopen(fname, "r");
  if(!fin) return NULL;
  char line[30000];
  while(fgets(line, 29999, fin)) {
    retval = (char*)realloc(retval, len + strlen(line) + 1024);
    strcpy(retval + len, line);
    len += (int)strlen(line);
  }
  fclose(fin);
  return retval;
}

inline unsigned char *ReadBinaryFile(const char *fname, int *len) {
  *len = 0;
  FILE *fin = fopen(fname, "rb");
  if(!fin) return NULL;
  int num;
  unsigned char *retval = (unsigned char*)malloc(1024);
  while((num=(int)fread(retval+*len, 1, 1024, fin)) > 0) {
    retval = (unsigned char*)realloc(retval, *len + num + 1024);
    *len += num;
  }
  fclose(fin);
  return retval;
}

inline bool SaveText(const char *str, const char *fname) {
  FILE *fout = fopen(fname, "w");
  if(fout) {
    fprintf(fout, "%s", str);
    fclose(fout);
    return true;
  }
  return false;
}


char *http_post_request(const char* hostname, const char* api, const char* parameters, int *err);


//  gamma.cpp -- computation of gamma function.
//      Algorithms and coefficient values from "Computation of Special
//      Functions", Zhang and Jin, John Wiley and Sons, 1996.
//
//  (C) 2003, C. Bond. All rights reserved.
//
// Returns gamma function of argument 'x'.
//
// NOTE: Returns 1e308 if argument is a negative integer or 0,
//      or if argument exceeds 171.
//
inline double gamma(double x)
{
    int i,k,m;
    double ga,gr,r,z;

    static double g[] = {
        1.0,
        0.5772156649015329,
       -0.6558780715202538,
       -0.420026350340952e-1,
        0.1665386113822915,
       -0.421977345555443e-1,
       -0.9621971527877e-2,
        0.7218943246663e-2,
       -0.11651675918591e-2,
       -0.2152416741149e-3,
        0.1280502823882e-3,
       -0.201348547807e-4,
       -0.12504934821e-5,
        0.1133027232e-5,
       -0.2056338417e-6,
        0.6116095e-8,
        0.50020075e-8,
       -0.11812746e-8,
        0.1043427e-9,
        0.77823e-11,
       -0.36968e-11,
        0.51e-12,
       -0.206e-13,
       -0.54e-14,
        0.14e-14};

    if (x > 171.0) return 1e308;    // This value is an overflow flag.
    if (x == (int)x) {
        if (x > 0.0) {
            ga = 1.0;               // use factorial
            for (i=2;i<x;i++) {
               ga *= i;
            }
         }
         else
            ga = 1e308;
     }
     else {
        if (fabs(x) > 1.0) {
            z = fabs(x);
            m = (int)z;
            r = 1.0;
            for (k=1;k<=m;k++) {
                r *= (z-k);
            }
            z -= m;
        }
        else
            z = x;
        gr = g[24];
        for (k=23;k>=0;k--) {
            gr = gr*z+g[k];
        }
        ga = 1.0/(gr*z);
        if (fabs(x) > 1.0) {
            ga *= r;
            if (x < 0.0) {
                ga = -M_PI/(x*ga*sin(M_PI*x));
            }
        }
    }
    return ga;
}

// Returns an array of size num that contains a random permutation of the numbers 0:(num-1)
inline int *RandPerm(int num) {
  int i, ind, tmp;
  int *a = (int*)malloc(num*sizeof(int));
  for(i = 0; i < num; i++)
    a[i] = i;
  for(i = 0; i < num; i++) {
    ind = i + rand()%(num-i); 
    tmp = a[i];
    a[i] = a[ind];
    a[ind] = tmp;
  }
  return a;
}

inline int *RandSplit(int num0, int numTotal) {
  int *perm = RandPerm(numTotal);
  int *retval = (int*)malloc(numTotal*sizeof(int));
  int i;
  for(i = 0; i < num0; i++)
    retval[perm[i]] = 0;
  for(; i < numTotal; i++)
    retval[perm[i]] = 1;
  free(perm);
  return retval;
}


template <typename T>
inline T *Create1DArray(int r) {
  T *retval = (T*)malloc(sizeof(T)*r);
  for(int i = 0; i < r; i++) 
    retval[i] = 0;
  return retval;
}

template <typename T>
inline T **Create2DArray(int r, int c) {
  T **retval = (T**)malloc((sizeof(T*)+sizeof(T)*c)*r);
  T *ptr = (T*)(retval+r);
  for(int i = 0; i < r; i++, ptr += c) {
    retval[i] = ptr;
    for(int j = 0; j < c; j++) 
      retval[i][j] = 0;
  }
  return retval;
}

template <typename T>
inline T ***Create3DArray(int d1, int d2, int d3) {
  T ***retval = (T***)malloc((sizeof(T**)+(sizeof(T*)+sizeof(T)*d3)*d2)*d1);
  T **ptr = (T**)(retval+d1);
  T *ptr2 = (T*)(ptr+d1*d2);
  for(int i = 0; i < d1; i++, ptr += d2) {
    retval[i] = ptr;
    for(int j = 0; j < d2; j++, ptr2 += d3) {
      retval[i][j] = ptr2;
      for(int k = 0; k < d3; k++) 
	retval[i][j][k] = 0;
    }
  }
  return retval;
}

template <typename T>
inline T ****Create4DArray(int d1, int d2, int d3, int d4) {
  T ****retval = (T****)malloc((sizeof(T***)+(sizeof(T**)+(sizeof(T*)+sizeof(T)*d4)*d3)*d2)*d1);
  T ***ptr = (T***)(retval+d1);
  T **ptr2 = (T**)(ptr+d1*d2);
  T *ptr3 = (T*)(ptr+d1*d2*d3);
  for(int i = 0; i < d1; i++, ptr += d2) {
    retval[i] = ptr;
    for(int j = 0; j < d2; j++, ptr2 += d3) {
      retval[i][j] = ptr2;
      for(int k = 0; k < d3; k++, ptr3 += d4) {
	retval[i][j][k] = ptr3;
	for(int l = 0; l < d4; l++) 
	  retval[i][j][k][l] = 0;
      }
    }
  }
  return retval;
}

template <typename T>
inline bool Save1DArray(FILE *fout, T *m, int d1, const char *name, const char *header) {
  int n = 1, c = sizeof(T), l = strlen(name), l2 = strlen(header);
  return fwrite(&l, sizeof(int), 1, fout) &&
    fwrite(name, sizeof(char), l, fout) &&
    fwrite(&l2, sizeof(int), 1, fout) &&
    fwrite(header, sizeof(char), l2, fout) &&
    fwrite(&n, sizeof(int), 1, fout) &&
    fwrite(&d1, sizeof(int), 1, fout) &&
    fwrite(&c, sizeof(int), 1, fout) &&
    fwrite(m, sizeof(T), d1, fout);
}

template <typename T>
inline bool Save2DArray(FILE *fout, T **m, int d1, int d2, const char *name, const char *header) {
  int n = 2, c = sizeof(T), l = strlen(name), l2 = strlen(header);
  if(!(fwrite(&l, sizeof(int), 1, fout) &&
       fwrite(name, sizeof(char), l, fout) &&
       fwrite(&l2, sizeof(int), 1, fout) &&
       fwrite(header, sizeof(char), l2, fout) &&
       fwrite(&n, sizeof(int), 1, fout) &&
       fwrite(&d1, sizeof(int), 1, fout) &&
       fwrite(&d2, sizeof(int), 1, fout) &&
       fwrite(&c, sizeof(int), 1, fout)))
    return false;
  for(int i = 0; i < d1; i++)
    if(!fwrite(m[i], sizeof(T), d2, fout))
      return false;
  return true;
}

template <typename T>
inline bool Save3DArray(FILE *fout, T ***m, int d1, int d2, int d3, const char *name, const char *header) {
  int n = 3, c = sizeof(T), l = strlen(name), l2 = strlen(header);
  if(!(fwrite(&l, sizeof(int), 1, fout) &&
       fwrite(name, sizeof(char), l, fout) &&
       fwrite(&l2, sizeof(int), 1, fout) &&
       fwrite(header, sizeof(char), l2, fout) &&
       fwrite(&n, sizeof(int), 1, fout) &&
       fwrite(&d1, sizeof(int), 1, fout) &&
       fwrite(&d2, sizeof(int), 1, fout) &&
       fwrite(&d3, sizeof(int), 1, fout) &&
       fwrite(&c, sizeof(int), 1, fout)))
    return false;
  for(int i = 0; i < d1; i++)
    for(int j = 0; j < d2; j++)
      if(!fwrite(m[i][j], sizeof(T), d3, fout))
	return false;
  return true;
}


template <typename T>
inline bool Save4DArray(FILE *fout, T ****m, int d1, int d2, int d3, int d4, const char *name, const char *header) {
  int n = 4, c = sizeof(T), l = strlen(name), l2 = strlen(header);
  if(!(fwrite(&l, sizeof(int), 1, fout) &&
       fwrite(name, sizeof(char), l, fout) &&
       fwrite(&l2, sizeof(int), 1, fout) &&
       fwrite(header, sizeof(char), l2, fout) &&
       fwrite(&n, sizeof(int), 1, fout) &&
       fwrite(&d1, sizeof(int), 1, fout) &&
       fwrite(&d2, sizeof(int), 1, fout) &&
       fwrite(&d3, sizeof(int), 1, fout) &&
       fwrite(&d4, sizeof(int), 1, fout) &&
       fwrite(&c, sizeof(int), 1, fout)))
    return false;
  for(int i = 0; i < d1; i++)
    for(int j = 0; j < d2; j++)
      for(int k = 0; k < d3; k++)
	if(!fwrite(m[i][j][k], sizeof(T), d4, fout))
	  return false;
  return true;
}

inline void SaveMatlabImport(const char *fname) {
  char fname2[1000];
  ExtractFilename(fname, fname2);
  StripFileExtension(fname2);
  strcat(fname2, ".bin");

  FILE *fout = fopen(fname, "w");
  assert(fout);
  fprintf(fout, "fid=fopen('%s', 'r');\n", fname2);
  fprintf(fout, "num_mats=fread(fid,1,'int32');\n");
  fprintf(fout, "for t=1:num_mats\n");
  fprintf(fout, "  mat_l=fread(fid,1,'int32');\n");
  fprintf(fout, "  mat_name=fread(fid,mat_l,'schar');\n");
  fprintf(fout, "  mat_l2=fread(fid,1,'int32');\n");
  fprintf(fout, "  mat_header=fread(fid,mat_l2,'schar');\n");
  fprintf(fout, "  mat_n=fread(fid,1,'int32');\n");
  fprintf(fout, "  mat_d=fread(fid,mat_n,'int32');\n");
  fprintf(fout, "  mat_c=fread(fid,1,'int32');\n");
  fprintf(fout, "  mat=zeros(mat_d');\n");
  fprintf(fout, "  if mat_n==1, mat=fread(fid,mat_d(1),char(mat_header'));\n");
  fprintf(fout, "  elseif mat_n==2, for j=1:mat_d(1) mat(j,:)=fread(fid,mat_d(2),char(mat_header')); end\n");
  fprintf(fout, "  elseif mat_n==3, for j=1:mat_d(1) for k=1:mat_d(2) mat(j,k,:)=fread(fid,mat_d(3),char(mat_header')); end; end\n");
  fprintf(fout, "  elseif mat_n==4, for j=1:mat_d(1) for k=1:mat_d(2) for l=1:mat_d(3) mat(j,k,l,:)=fread(fid,mat_d(4),char(mat_header')); end; end; end\n");
  fprintf(fout, "  end\n");
  fprintf(fout, "  assignin('base', char(mat_name'), mat);\n");
  fprintf(fout, "end\n");
  fclose(fout);
}

inline void CreateDirectoryIfNecessary(const char *dirName, int permissions=0) {
#ifndef WIN32
  char str[400];
  sprintf(str, "mkdir %s", dirName);
  if(system(str) < 0) fprintf(stderr, "%s failed\n", str);
  if(permissions) {
    sprintf(str, "chmod %d %s", permissions, dirName);
    if(system(str) < 0) fprintf(stderr, "%s failed\n", str);
  }
#else
  CreateDirectory(dirName, NULL);
#endif
}


inline void ScanDir(const char *dirName, const char *match, const char *outName) {
  char matches[100][400], str[400];
  char tmpMatch[1000];
  char *ptr;
  int numMatches=0;

  strcpy(tmpMatch, match);
  while((ptr=strtok(numMatches ? NULL : tmpMatch,"|"))) {
    if(numMatches) sprintf(str, "ls %s/*%s >> %s", dirName, ptr, outName);
    else sprintf(str, "ls %s/*%s > %s", dirName, ptr, outName);
    if(system(str) < 0) fprintf(stderr, "ScanDir failed\n");
    strcpy(matches[numMatches++], ptr);
  }
}

#endif
