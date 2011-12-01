#ifndef __FIT_MODEL_H
#define __FIT_MODEL_H


#include "blob.h"
#include "spine_features.h"
#include "train.h"

struct _SpinePointTransitions;

typedef enum {
  FIT_FRAME_BY_FRAME = 0,
  FIT_GREEDY_ONE_PASS,
  FIT_GREEDY_TWO_PASS,
  FIT_GLOBALLY_OPTIMAL
} FitMethod;

// Parameters used to fit a spine
typedef struct _FitParams {
  // the number of lines that makeup the spine (num_spine_points-1)
  int num_spine_lines; 
  
  // the number of discrete orientations to search through
  int num_orientations;
  double *orientations, *cos_orientations, *sin_orientations;

  int pixel_granularity;  // granularity in pixels for candidate spine point locations

  double expected_length;  // expected length of the worm, in world coordinates

  FitMethod fit_method;

  // Multiply points in frame_contours by this to convert from world 
  // units to pixel units 
  double world_to_pixel_scale;

  ClassifierMethod classifier_method;
  int is_multiclass;
  int feature_diff_frames;

  // These parameters encode the relative cost of each feature type
  double mu[F_NUM_FEATURES];
  double var[F_NUM_FEATURES];

  // Prune the search space by eliminating things with weird geometry
  double max_position_change;  // max distance in world coordinates between consecutive spine points
  double max_angle_change;     // maximum changein angle in radians between consecutive spine points
  double max_width_change;     // maximum width ratio between two adjacent spine cross sections
  double max_width;            // maximum width of any cross section in world coordinates
  double max_endpoint_length;  // maximum length in world coordinates of the spine segment next to an endpoint
  double max_outlier_dist;     // maximum outlier distance in pixels of contour points
  double max_endpoint_ratio;   // maximum ratio of endpoint segment length to the adjacent segment length
  int blob_min_pixel_width;    // the minimum width in pixels of any worm cross section

  double speed_estimate_window;   // the amount of time 
} FitParams;




#define ORI(ori) (params->orientations[ori])
#define COS_ORI(ori) (params->cos_orientations[ori])
#define SIN_ORI(ori) (params->sin_orientations[ori])
#define COS_N_ORI(ori) (params->cos_orientations[ori])
#define SIN_N_ORI(ori) (-params->sin_orientations[ori])

// World to pixel coordinates
#define W2P(x,i) (((x)+offsets[i])*params->world_to_pixel_scale)

// Pixel to pworld coordinates
#define P2W(x,i) ((x)/params->world_to_pixel_scale-offsets[i])


struct _IplImage;

struct _IplImage *draw_blob(struct _IplImage *img, Blob *b, double *offsets, FitParams *params, double zoom, BehaviorGroups *beh);

double fit_model(BlobSequence *s, FitParams *params, int debug, void (*on_frame)(int t), int *keep_going);
double fit_model_one_frame(Blob *b, Blob *prev,  FitParams *params, int free_mem, Blob *result, int debug);

struct _IplImage *binary_mask(double *contour, int num_pts, double *offsets, FitParams *params);
SpinePointLocation *get_legal_internal_points(double *contour, int num_pts, 
						 double *offsets, struct _IplImage *blob, 
						 int *num_legal_internal_points, 
						 FitParams *params);
void compute_intersect_points(SpinePointLocation *pt, double *contour, int num_pts, 
						 FitParams *params);
double guess_world_to_pixel_scale(double *contour, int num_pts);


void compute_deterministic_spine_attributes(SpinePointLocation *pt, SpinePointLocation *prev, 
					    SpinePointLocation *next, double *contour, int num_pts, 
					    FitParams *params, int find_orientation);
void compute_deterministic_spine_attributes_blob(Blob *blob, FitParams *params, int find_orientation);
void compute_deterministic_spine_attributes_blob_sequence(BlobSequence *s, FitParams *params, int find_orientation);
void flip_spine(Blob *blob, FitParams *params);
void flip_spine_sequence(BlobSequence *blob, FitParams *params);



FitParams default_parameters();

#endif
