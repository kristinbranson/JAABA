#ifndef __BLOB_H
#define __BLOB_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <assert.h>
#include <string.h>

#ifndef INFINITY
#define INFINITY HUGE_VAL
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifdef WIN32
#include <float.h>
#ifndef isnan
#define isnan(a) _isnan(a)
#endif
#endif
//#define isnan(a) (0)


#define MAX_BEHAVIOR_GROUPS 100
#define MAX_BEHAVIOR_VALUES 100


class CvStatModel;
//class CvFileStorage;

/*
 * A spine is a line passing through the center of the worm,
 * and consists of a sequence of num_skeletal_points control points.
 * A SpinePointLocation specifies where to place a control point
 * of a spine.  This is parameterized by an (x,y) pixel location
 * and the orientation of the line orthogonal to the spine (e.g.
 * a cross-section of the worm) 
 */
typedef struct {
  // A possible location of a spine point
  double x, y;           // x-y coordinates in world coordinates
  int orientation;       // orientation of the worm cross-section.  The actual orientation
                         // is orientation * 2pi / num_orientations

  // Parameters that can be computed directly from x,y,orientation
  // and could be useful in evaluating a cost function
  double width;          // width of the cross-section
  int p1, p2, p3, p4;    // indices of the four points on the contour that intersect with the
                         // ray beginning at x,y in the direction of 0, 180, 90, and 270 degrees
                         // from the angle defined by orientation 
  double d3, d4;         // the length of the rays for p3 and p4

} SpinePointLocation;


/*
 * A blob in a single frame
 */
typedef struct _Blob {
  // The raw data, as returned by Rex's tracker
  double *contour;
  int num_pts;
  double frame_time;

  double *fixed_contour;
  int fixed_num_pts;

  // The fit model
  SpinePointLocation *model_pts;
  int num_model_pts;
  double *features;  // the spine features for each model point 

  int is_manual;  // the model was annotated by hand (as opposed to by a machine)
  int behaviors[MAX_BEHAVIOR_GROUPS];   // the behavior of the blob in this frame
  int meta_behaviors[MAX_BEHAVIOR_GROUPS];
  int is_good;    // did spine fitting succeed on this blob?

  int on_pt;

  void *data;     // extra data for custom storage
} Blob;


/*
 * An array of blobs over some tracked time sequence
 */
typedef struct _BlobSequence {
  int num_frames;
  Blob *frames;

  char from[400];
  char fname[400];
} BlobSequence;


typedef struct _MultiBlobSequence {
  int num_blobs;
  BlobSequence **blobs;

  char fname[400];
} MultiBlobSequence;



// Single behavior class. 
typedef struct _BehaviorValue {
  char name[400];
  char abbreviation[400];
  char regExp[400];
  unsigned int color;

  CvStatModel *classifier;
  int classifier_method;
} BehaviorValue;

// Group of mutually-exclusive behaviors that should be detected dependently. 
typedef struct _BehaviorGroup {
  char name[400];
  BehaviorValue *values;
  int num_values;

  CvStatModel *classifier;
  int classifier_method;
  int is_multiclass;
} BehaviorGroup;

// All groups of behaviors -- this represents all behavior labels.
typedef struct _BehaviorGroups {
  BehaviorGroup *behaviors;
  int num;
  CvStatModel *classifier;
  int classifier_method;
  int is_multiclass;
} BehaviorGroups;



#define  SQR(x) ((x)*(x))
#define my_max(x,y) ((x)>=(y) ? (x) : (y))
#define my_max3(x,y,z) (my_max(my_max(x,y),z))
#define my_min(x,y) ((x)<=(y) ? (x) : (y))
#define my_abs(x) ((x) < 0 ? (-(x)) : (x))
#define ANGLE_DIST(dtheta) \
  (dtheta > M_PI ? SQR(dtheta-2*M_PI) : (dtheta < -M_PI ? SQR(dtheta+2*M_PI) : SQR(dtheta)))
#define ANGLE_MOD(theta) \
  (theta > M_PI ? (theta-2*M_PI) : (theta < -M_PI ? (theta+2*M_PI) : (theta)))




BlobSequence *import_blob_sequence(const char *fname, int num_spine_points);
void free_blob_sequence(BlobSequence *s);
BlobSequence *load_blob_sequence(const char *fname, BehaviorGroups *behaviors);
void save_blob_sequence(const char *fname, BlobSequence *b, BehaviorGroups *behaviors, 
			int num_model_pts, int num_orientations);

void blob_get_bounding_box(Blob **blobs, int num_blobs, double *min_x, 
			   double *min_y, double *max_x, double *max_y);
void free_blob(Blob *b);

MultiBlobSequence *import_multi_blob_sequence(const char *dir_name, int num_spline_points);
void free_multi_blob_sequence(MultiBlobSequence *m);
void save_multi_blob_sequence(MultiBlobSequence *m, BehaviorGroups *behaviors, int num_model_pts, int num_orientations, int save_blobs);
MultiBlobSequence *load_multi_blob_sequence(const char *fname, BehaviorGroups *behaviors);

void free_behaviors(BehaviorGroups *beh);
BehaviorGroups *load_behaviors(const char *fname, const char *classifier_dir);
int save_behaviors(const char *fname, BehaviorGroups *behaviors);
int find_behavior(BehaviorGroups *behaviors, int *inds, const char *str);
void save_blob_behavior_bouts(const char *fname, BlobSequence *b, BehaviorGroups *behaviors, int is_meta);



char *BehaviorCharString(BlobSequence *b, BehaviorGroups *beheviors, int beh);
void RegularExpressionCharString(const char *regExp, char *regA, BehaviorGroup *behaviors);
void DetectRegularExpression(const char *regExp, char *str, BlobSequence *b, int beh, int beh_val);
void ExtractMetaBehaviors(BlobSequence *b, BehaviorGroups *behaviors, BehaviorGroups *meta_behaviors);

void strip_extension(char *fname);

inline void chomp(char *str) {
  int l = (int)strlen(str)-1;
  while(l >= 0 && isspace(str[l])) {
    str[l] = '\0';
    l--;
  }
}

inline char *StringCopy(const char *str) {
  char *retval = (char*)malloc(strlen(str)+1);
  strcpy(retval, str);
  return retval;
}

inline int FileExists(const char *fname) { 
  FILE *fin = fopen(fname, "r"); 
  if(fin) fclose(fin);
  return fin != NULL;
}

inline double point_dist(double x1, double y1, double x2, double y2) {
  return sqrt(SQR(x2-x1)+SQR(y2-y1));
}


// Find the squared distance from (x,y) to the line (x+vx*t,y+vy*t), where [vx,vy] is a unit vector
inline double point_line_dist(double x, double y, double x0, double y0, double vx, double vy, double *t, double *side) {
  double vx1 = x-x0, vy1 = y-y0;
  double p = vx1*vx + vy1*vy;
  double dx = vx1-p*vx, dy = vy1-p*vy;
  if(t) *t = p;
  if(side) *side = vx1*vy - vy1*vx;
  
  return SQR(dx) + SQR(dy);
}


// Find the nearest point to (x,y) in the point array contour
inline int nearest_point(double x, double y, double *contour, int num_pts) {
  double best = INFINITY, dist;
  int retval = -1, i;
  
  for(i = 0; i < num_pts; i++) {
    dist = SQR(contour[2*i]-x)+SQR(contour[2*i+1]-y);
    if(dist < best) {
      best = dist;
      retval = i;
    }
  }
  return retval;
}
extern int g_debug;

// Find the closest contour point that "intersects" the line (x+vx*t,y+vy*t) in the point array contour, where [vx,vy] 
// is a unit vector.  Works by stepping through each contour line segment and seeing if it intersects the line.  If multiple
// such contour segments are found, it chooses the one at minimum distance.  It then selects the contour point on that segment
// with minimum distance to the line
inline int nearest_to_line(double x, double y, double vx, double vy, int forward, double *contour, int num_pts) {
  double best = INFINITY, dist, t, last_dist, last_t, side, last_side;
  int retval = -1, i, last;
  
  last = num_pts-1;
  last_dist = point_line_dist(contour[2*last], contour[2*last+1], x, y, vx, vy, &last_t, &last_side);
  for(i = 0; i < num_pts; i++) {
    dist = point_line_dist(contour[2*i], contour[2*i+1], x, y, vx, vy, &t, &side);
    if(((side <= 0 && last_side >= 0) || (side >= 0 && last_side <= 0)) && ((forward > 0 && t >= 0) || (forward < 0 && t <= 0) || forward == 0)) {
      if(my_abs(t) < best) {
	best = my_abs(t);
	retval = last_dist < dist ? last : i;
      }
    }
    last_side = side;
    last_t = t;
    last_dist = dist;
    last = i;
  }
  return retval;
}


#endif
