#ifndef __SPINE_FEATURE_H
#define __SPINE_FEATURE_H

#include "blob.h"

struct _FitParams;

// Let p_s and p_{s-1} be consecutive spine points and let the vector between
//   the two points be l_s
// Let length_l and length_r be the length of the two line segments going from
//   the left and right contour intersection points to the tip of the worm
// Let geo_l and geo_r be the geodesic length of the curve going from
//   the left and right contour intersection points to the tip of the worm
typedef enum  {
  F_WIDTH = 0,               // the width of this cross-section

  /*********** Head/Tail unary features ***********************/
  F_HEAD_LENGTH,             // Length of the head spine segment
  F_TAIL_LENGTH,             // Length of the tail spine segment
  F_HEAD_SYMMETRY,           // (length_l-length_r)/(length_l+length_r)
  F_TAIL_SYMMETRY,
  F_HEAD_GEODESIC_SYMMETRY, // (geo_l-geo_r)/(geo_l+geo_r)
  F_TAIL_GEODESIC_SYMMETRY, // (geo_l-geo_r)/(geo_l+geo_r)
  F_HEAD_PERCENT_OUTLIERS,  // The fraction of contour points lying outside the head region
  F_TAIL_PERCENT_OUTLIERS,
  F_HEAD_GEOMETRY,          // cos of the angle between the head segment and the cross-section
  F_TAIL_GEOMETRY,
  
  /*********** Time unary features ***********************/
  F_WIDTH_VELOCITY,           // Change in width from the previous time step
  F_FORWARD_VELOCITY,         // Change in the global position along the spine axis
  F_SIDEWAYS_VELOCITY,        // Change in the global position, orthogonal to the spine axis
  F_ANGULAR_VELOCITY,         // Change in the absolute angle from the previous time step

  F_NUM_UNARY_FEATURES,



  /*********** Spine Position Pairwise features ***********************/
  F_PERCENT_OUTLIERS,  // the (normalized) number of contour points outside too far away 
                       // from the contour predicted by the spine
  F_WIDTH_CHANGE,   // p_s->width - p{s-1}->width
  F_ANGLE_CHANGE,   // p_s->angle - p{s-1}->angle
  F_LENGTH,         // dist((p_s->x,pt->y),(p_{s-1}->x,p{s-1}->y))
  F_LEFT_PARALLEL,  // Decompose the line segment on the contour to the left
  F_LEFT_ORTHOGONAL,//   l_s into parallel and orthogonal components to l_s
  F_RIGHT_PARALLEL,
  F_RIGHT_ORTHOGONAL,
  F_SYMMETRY,       // F_LEFT_ORTHOGONAL+F_RIGHT_ORTHOGONAL (0 for symmetric contours)
  F_GEOMETRY1,
  F_GEOMETRY2,
  
  /*********** Spine Position Pairwise features, time matters **********/
  F_SEGMENT_ANGULAR_VELOCITY, // Change in F_ANGLE_CHANGE between consecutive time steps
  F_LENGTH_VELOCITY,          // Change in F_LENGTH between consecutive time steps

  F_NUM_FEATURES
} SpineFeatureTypes;

typedef enum {
  F_G_WIDTH_LENGTH_RATIO = 0,
  F_G_FORWARD_VELOCITY,
  F_G_SIDEWAYS_VELOCITY,
  F_G_MAX_LEFT_CURVATURE,
  F_G_MAX_RIGHT_CURVATURE,
  F_G_CURVE1,
  F_G_CURVE2,
  F_G_CURVE3,
  F_G_CURVE4,
  F_G_CURVE5,
  F_G_CURVE6,
  F_NUM_GOOD_GLOBAL_FEATURES,

  F_G_WIDTH1,
  F_G_WIDTH2,
  F_G_WIDTH3,
  F_G_WIDTH4,
  F_G_WIDTH5,
  F_G_WIDTH6,
  F_G_LENGTH,
  F_G_MAX_WIDTH,
  F_G_AVE_WIDTH,
  F_G_AVE_X,
  F_G_AVE_Y,
  F_G_X1,
  F_G_X2,
  F_G_X3,
  F_G_X4,
  F_G_X5,
  F_G_X6,
  F_G_Y1,
  F_G_Y2,
  F_G_Y3,
  F_G_Y4,
  F_G_Y5,
  F_G_Y6,
  F_NUM_GLOBAL_FEATURES
} GlobalSpineFeatureTypes;

#define F_NUM_GLOBAL_POINTS 6


#define NUM_PAIRWISE_FEATURES (F_NUM_FEATURES-F_NUM_UNARY_FEATURES-1)
#define PAIR_SIZE (sizeof(int)+NUM_PAIRWISE_FEATURES*sizeof(double))

#define SQR_FEATURES(b) ((b)->features + (b)->num_model_pts*F_NUM_FEATURES)
#define GLOBAL_FEATURES(b) ((b)->features + 2*(b)->num_model_pts*F_NUM_FEATURES)
#define NORMALIZED_FEATURES(b) ((b)->features + 2*(b)->num_model_pts*F_NUM_FEATURES+F_NUM_GLOBAL_FEATURES)
#define NORMALIZED_SQR_FEATURES(b) ((b)->features + 3*(b)->num_model_pts*F_NUM_FEATURES+F_NUM_GLOBAL_FEATURES)
#define NORMALIZED_GLOBAL_FEATURES(b) ((b)->features + 4*(b)->num_model_pts*F_NUM_FEATURES+F_NUM_GLOBAL_FEATURES)

typedef struct _SpinePointTransitions {
  SpinePointLocation *pt;
  int ind;
  int num;
  double *unary_features;
  unsigned char *pairwise_features;
} SpinePointTransitions;


SpinePointTransitions **compute_pairwise_transitions(SpinePointLocation *legal_points, 
						     int num_legal, double *contour, int num_contour_points, 
						     struct _FitParams *params);
void free_pairwise_spine_transitions(SpinePointTransitions **t, int num);

double eval_unary_cost(SpinePointLocation *pt, int s, double *feat, SpinePointLocation *prev_pt, double dt_inv, 
		       double *mu, double *var, struct _FitParams *params);
double eval_pairwise_cost(SpinePointLocation *pt, int s, double *feat, SpinePointLocation *pt_prev, 
			  double *prev_features, double dt_inv, double *mu, double *var, struct _FitParams *params);

int compute_pairwise_features(SpinePointLocation *t1, SpinePointLocation *t2, 
			       double *features, double *contour, int num_contour_points, 
			       struct _FitParams *params);
int compute_unary_features(SpinePointLocation *pt, double *features, double *contour, int num_contour_points, 
			   struct _FitParams *params);

void print_features(double *feat, int s, bool has_prev, double *mu, double *var, struct _FitParams *params);
void print_global_features(double *feat);
void compute_all_global_features(BlobSequence *s, struct _FitParams *params, int debug);
void compute_all_features(BlobSequence *s, struct _FitParams *params, int debug);
void compute_global_features(BlobSequence *s, int t, struct _FitParams *params, int debug);
void compute_features(Blob *b, Blob *prev, double dt_inv, struct _FitParams *params, int debug);

extern const char *g_global_feature_names[F_NUM_GOOD_GLOBAL_FEATURES];

#endif
