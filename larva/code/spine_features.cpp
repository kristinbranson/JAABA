#include "spine_features.h"
#include "fit_model.h"

#include <cv.h>
#include <highgui.h>
#include <ml.h>  

int g_debug = 0;

#define MYSIGNUM(x) (-(x<0)+(x>0)) // CSC: fast signum function

#define LINE_SIDE(x, y, x0, y0, vx, vy) ((x-x0)*vy - (y-y0)*vx)
  
#define GET_CONTOUR_DIRECTION(p1, p3, dir, diff)			\
  d1 = (p3+num_contour_points-p1) % num_contour_points;			\
  d2 = (p1+num_contour_points-p3) % num_contour_points;			\
  if(d1 < d2) {								\
    dir = 1;								\
    diff = d1;								\
  } else {								\
    dir = -1;								\
    diff = d2;								\
  }

#define COUNT_ENDPOINT_OUTLIERS(p1, p3, dir, p2, p4)			\
  i = p1;								\
  while(i != p3) {							\
    dd1 = LINE_SIDE(contour[2*i], contour[2*i+1], pt->x, pt->y,		\
		   COS_ORI(pt->orientation), SIN_ORI(pt->orientation));	\
    dd2 = LINE_SIDE(contour[2*i], contour[2*i+1], contour[2*p3], contour[2*p3+1], \
		   COS_ORI(pt->orientation), SIN_ORI(pt->orientation));	\
    dd3 = LINE_SIDE(contour[2*i], contour[2*i+1], contour[2*p1], contour[2*p1+1], \
		   COS_N_ORI(pt->orientation), SIN_N_ORI(pt->orientation));	\
    dd4 = LINE_SIDE(contour[2*i], contour[2*i+1], contour[2*p2], contour[2*p2+1], \
		   COS_N_ORI(pt->orientation), SIN_N_ORI(pt->orientation)); \
    i += dir;								\
    if((dd1 < 0 && dd2 < 0) ||  (dd1 > 0 && dd2 > 0) || (dd3 < 0 && dd4 < 0) ||  (dd3 > 0 && dd4 > 0)) \
      num++;								\
    if(i < 0) i += num_contour_points;					\
    else if(i >= num_contour_points) i -= num_contour_points;		\
  }		
 
#define COUNT_OUTLIERS(p1, p3, vx, vy, dir)				\
  i = p1;								\
  while(i != p3) {							\
    d = point_line_dist(contour[2*i], contour[2*i+1], contour[2*p1], contour[2*p1+1], vx, vy, NULL, NULL); \
    i += dir;								\
    if(d > outlier_dist) 						\
      num++;								\
    if(i < 0) i += num_contour_points;					\
    else if(i >= num_contour_points) i -= num_contour_points;		\
  } 

//if(g_debug) fprintf(stderr, "*%d: %f %f\n", i, d, outlier_dist);	
//} else if(g_debug) fprintf(stderr, "%d: %f\n", i, d);		

int compute_unary_features(SpinePointLocation *pt, double *features, double *contour, int num_contour_points, 
			    FitParams *params) {
  int d1, d2, dir1, diff1, dir2, diff2, i, num;
  double dd1, dd2, dd3, dd4;
  double dx, dy, len;

  features[F_WIDTH] = pt->width;
  if(features[F_WIDTH] > params->max_width) {
    features[F_WIDTH] = INFINITY;
    return 0;
  }

  // Compute head features
  if(pt->p4 >= 0 && pt->p1 >= 0 && pt->p2 >= 0 && pt->p3 >= 0) {
    features[F_HEAD_LENGTH] = pt->d4;
    features[F_HEAD_SYMMETRY] = (point_dist(contour[2*pt->p1], contour[2*pt->p1+1], 
					    contour[2*pt->p4], contour[2*pt->p4+1]) -
				 point_dist(contour[2*pt->p2], contour[2*pt->p2+1], 
					    contour[2*pt->p4], contour[2*pt->p4+1]));
    
    
    GET_CONTOUR_DIRECTION(pt->p1, pt->p4, dir1, diff1);
    GET_CONTOUR_DIRECTION(pt->p2, pt->p4, dir2, diff2);
    features[F_HEAD_GEODESIC_SYMMETRY] = (diff1 - diff2) / (double)num_contour_points;
    num = 0;
    COUNT_ENDPOINT_OUTLIERS(pt->p1, pt->p4, dir1, pt->p2, pt->p3);
    COUNT_ENDPOINT_OUTLIERS(pt->p2, pt->p4, dir2, pt->p1, pt->p3);
    features[F_HEAD_PERCENT_OUTLIERS] = num / (double)num_contour_points;
    dx = contour[2*pt->p4] - pt->x;
    dy = contour[2*pt->p4+1] - pt->y;
    len = sqrt(SQR(dx)+SQR(dy));
    features[F_HEAD_GEOMETRY] = (COS_ORI(pt->orientation)*dx/len + SIN_ORI(pt->orientation)*dy/len);
  } else {
    features[F_HEAD_LENGTH] = features[F_HEAD_SYMMETRY] = features[F_HEAD_GEOMETRY] = features[F_HEAD_GEODESIC_SYMMETRY] = 
      features[F_HEAD_PERCENT_OUTLIERS] = INFINITY;
  }
  
  // Compute tail features
  if(pt->p3 >= 0 && pt->p1 >= 0 && pt->p2 >= 0 && pt->p4 >= 0) {
    features[F_TAIL_LENGTH] = pt->d3;
    features[F_TAIL_SYMMETRY] = (point_dist(contour[2*pt->p1], contour[2*pt->p1+1], 
					    contour[2*pt->p3], contour[2*pt->p3+1]) -
				 point_dist(contour[2*pt->p2], contour[2*pt->p2+1], 
					    contour[2*pt->p3], contour[2*pt->p3+1]));
        
    GET_CONTOUR_DIRECTION(pt->p1, pt->p3, dir1, diff1);
    GET_CONTOUR_DIRECTION(pt->p2, pt->p3, dir2, diff2);
    features[F_TAIL_GEODESIC_SYMMETRY] = (diff1 - diff2) / (double)num_contour_points;
    num = 0;
    COUNT_ENDPOINT_OUTLIERS(pt->p1, pt->p3, dir1, pt->p2, pt->p4);
    COUNT_ENDPOINT_OUTLIERS(pt->p2, pt->p3, dir2, pt->p1, pt->p4);
    features[F_TAIL_PERCENT_OUTLIERS] = num / (double)num_contour_points;
    dx = contour[2*pt->p3] - pt->x;
    dy = contour[2*pt->p3+1] - pt->y;
    len = sqrt(SQR(dx)+SQR(dy));
    features[F_TAIL_GEOMETRY] = (COS_ORI(pt->orientation)*dx/len + SIN_ORI(pt->orientation)*dy/len);
  } else {
    features[F_TAIL_LENGTH] = features[F_TAIL_SYMMETRY] = features[F_TAIL_GEOMETRY] = features[F_TAIL_GEODESIC_SYMMETRY] = 
      features[F_TAIL_PERCENT_OUTLIERS] = INFINITY;
  }

  return 1;
}

// Compute pairwise features, aborting early (return false if the two points
// t1 and t2 cannot be consecutive spine points)
int compute_pairwise_features(SpinePointLocation *t1, SpinePointLocation *t2, 
			       double *features, double *contour, int num_contour_points, 
			       FitParams *params) {
  double dx, dy, vx_l, vy_l, proj_l, len_l, vx_r, vy_r, proj_r, len_r, d;
  int dir1, dir2, d1, d2, diff1, diff2, num, i;
  double outlier_dist = SQR(params->max_outlier_dist/params->world_to_pixel_scale);

  if(t1->p1 < 0 || t2->p1 < 0 || t1->p2 < 0 || t2->p2 < 0)
    return 0;

  // Do a quick check to see if the two points might be close enough to be adjacent spine points
  dx = t2->x - t1->x;
  if(my_abs(dx) > params->max_position_change)
    return 0;
  dy = t2->y - t1->y;
  if(my_abs(dy) > params->max_position_change)
    return 0;

  // Check if the change in angle exceeds the maximum angle change
  features[F_ANGLE_CHANGE] = ANGLE_MOD(ORI(t2->orientation)-ORI(t1->orientation));
    // CSC: split signed value into absolute value and sign
//  features[CSCF_ANGLE_CHANGE_ABS] = abs(features[F_ANGLE_CHANGE]); // CSC
//  features(CSCF_ANGLE_CHANGE_SIGN] = MYSIGNUM(features[F_ANGLE_CHANGE]); // CSC
  if(my_abs(features[F_ANGLE_CHANGE]) > params->max_angle_change)
    return 0;

  // Check if the orientation of the spine line is reversed
  if(COS_ORI(t1->orientation)*dy - SIN_ORI(t1->orientation)*dx < 0)
    return 0;

  // Check if the change in width is too big
  features[F_WIDTH_CHANGE] = t2->width/t1->width;
  if(features[F_WIDTH_CHANGE] > params->max_width_change || 
     features[F_WIDTH_CHANGE] < 1/params->max_width_change)
    return 0;
  features[F_WIDTH_CHANGE] -= 1; // make the mean centered at zero

  // Check if the change in angle is too much, such that the two consecutive spine cross-section
  // segments intersect
  vx_l = contour[2*t2->p1] - contour[2*t1->p1];
  vy_l = contour[2*t2->p1+1] - contour[2*t1->p1+1];
  proj_l = -SIN_ORI(t1->orientation)*vx_l + COS_ORI(t1->orientation)*vy_l;
  //if(proj_l < 0) 
  //  return 0;
  vx_r = contour[2*t2->p2] - contour[2*t1->p2];
  vy_r = contour[2*t2->p2+1] - contour[2*t1->p2+1];
  proj_r = -SIN_ORI(t1->orientation)*vx_r + COS_ORI(t1->orientation)*vy_r;
  //if(proj_r < 0) 
  //  return 0;

  // Do an exact check if the two points are close enough to each other
  features[F_LENGTH] = sqrt(SQR(dx)+SQR(dy));
  if(features[F_LENGTH] > params->max_position_change)
    return 0;

  // Compute the symmetry scores.  features[F_LEFT_ORTHOGONAL] is 0 if the left contour segment
  //  is parallel to the spine segment.  features[F_LEFT_ORTHOGONAL]=-features[F_RIGHT_ORTHOGONAL]
  //  if the left and right contour segments are symmetric about the spine segment
  len_l = sqrt(SQR(vx_l)+SQR(vy_l));
  if(len_l > 0) {
    features[F_LEFT_PARALLEL] = proj_l/len_l-1;
    features[F_LEFT_ORTHOGONAL] = (COS_ORI(t1->orientation)*vx_l + 
				   SIN_ORI(t1->orientation)*vy_l)/len_l;
  } else {
    features[F_LEFT_PARALLEL] = 0;
    features[F_LEFT_ORTHOGONAL] = 0;
  }
  len_r = sqrt(SQR(vx_r)+SQR(vy_r));
  if(len_r > 0) {
    features[F_RIGHT_PARALLEL] = proj_r/len_r-1;
    features[F_RIGHT_ORTHOGONAL] = (COS_ORI(t1->orientation)*vx_r + 
				    SIN_ORI(t1->orientation)*vy_r)/len_r;
  } else {
    features[F_RIGHT_PARALLEL] = 0;
    features[F_RIGHT_ORTHOGONAL] = 0;
  }

  features[F_SYMMETRY] = features[F_LEFT_ORTHOGONAL] + features[F_RIGHT_ORTHOGONAL];

  dx /= features[F_LENGTH];  dy /= features[F_LENGTH];  
  features[F_GEOMETRY1] = (COS_ORI(t1->orientation)*dx + SIN_ORI(t1->orientation)*dy);
  features[F_GEOMETRY2] = (COS_ORI(t2->orientation)*dx + SIN_ORI(t2->orientation)*dy);

  
  // Number of points that are too far away from the predicted contour line segments
  vx_l /= len_l;  vx_r /= len_r;  vy_l /= len_l;  vy_r /= len_r;  
  num = 0;
  GET_CONTOUR_DIRECTION(t1->p1, t2->p1, dir1, diff1);
  GET_CONTOUR_DIRECTION(t1->p2, t2->p2, dir2, diff2);
  COUNT_OUTLIERS(t1->p1, t2->p1, vx_l, vy_l, dir1);
  COUNT_OUTLIERS(t1->p2, t2->p2, vx_r, vy_r, dir2);
  features[F_PERCENT_OUTLIERS] = num / (double)num_contour_points;

  return 1;
}

#define COST(ind) (SQR(feat[ind]-mu[ind]) / var[ind])

double eval_unary_cost(SpinePointLocation *pt, int s, double *feat, SpinePointLocation *prev_pt, double dt_inv,
		       double *mu, double *var, FitParams *params) {
  double cost = 0, dx, dy;

  if(feat[F_WIDTH] > params->max_width)
    return INFINITY;

  if(s == 0) {
    // Head spine point
    if(feat[F_HEAD_LENGTH] >= params->max_endpoint_length)
      return INFINITY;
    cost += COST(F_HEAD_LENGTH) + COST(F_HEAD_SYMMETRY) + COST(F_HEAD_GEODESIC_SYMMETRY) +
      COST(F_HEAD_GEOMETRY) + COST(F_HEAD_PERCENT_OUTLIERS);
  } else if(s == params->num_spine_lines) {
    // Tail spine point
    if(feat[F_TAIL_LENGTH] >= params->max_endpoint_length)
      return INFINITY;
    cost += COST(F_TAIL_LENGTH) + COST(F_TAIL_SYMMETRY) + COST(F_TAIL_GEODESIC_SYMMETRY) +
      COST(F_TAIL_GEOMETRY) + COST(F_TAIL_PERCENT_OUTLIERS);
  } 

  cost += COST(F_WIDTH);
  
  // Compute time dependent features
  if(prev_pt) {
    dx = pt->x-prev_pt->x;
    dy = pt->y-prev_pt->y;
    feat[F_FORWARD_VELOCITY] = (-SIN_ORI(prev_pt->orientation)*dx +
				COS_ORI(prev_pt->orientation)*dy) * dt_inv;
    feat[F_SIDEWAYS_VELOCITY] = (COS_ORI(prev_pt->orientation)*dx +
				 SIN_ORI(prev_pt->orientation)*dy) * dt_inv;
    feat[F_WIDTH_VELOCITY] = (pt->width - prev_pt->width) * dt_inv;
    feat[F_ANGULAR_VELOCITY] = ANGLE_MOD(ORI(pt->orientation) - 
					 ORI(prev_pt->orientation)) * dt_inv;
    if(s == 0)
      cost += COST(F_FORWARD_VELOCITY) + COST(F_SIDEWAYS_VELOCITY) + COST(F_ANGULAR_VELOCITY);
    
    cost += COST(F_WIDTH_VELOCITY);
  }

  if(!(cost < INFINITY)){
    // KB WARNING
    fprintf(stderr,"Warning: eval_unary_cost returning INFINITY");
  }

  return cost;
}

double eval_pairwise_cost(SpinePointLocation *pt, int s, double *feat, 
			  SpinePointLocation *pt_prev, double *features_prev_frame, double dt_inv, 
			  double *mu, double *var, FitParams *params) {
  double cost;

  if(s == 1 && pt_prev->d4 > feat[F_LENGTH]*params->max_endpoint_ratio)
    return INFINITY;
  else if(s == params->num_spine_lines && pt->d3 > feat[F_LENGTH]*params->max_endpoint_ratio)
    return INFINITY;

  cost = COST(F_WIDTH_CHANGE) + COST(F_ANGLE_CHANGE) + COST(F_LENGTH) + COST(F_LEFT_PARALLEL) + 
    COST(F_LEFT_ORTHOGONAL) + COST(F_RIGHT_PARALLEL) + COST(F_RIGHT_ORTHOGONAL) + COST(F_SYMMETRY) +
    COST(F_PERCENT_OUTLIERS) + COST(F_GEOMETRY1) + COST(F_GEOMETRY2);


  // Compute time-dependent features
  if(features_prev_frame) {
    feat[F_SEGMENT_ANGULAR_VELOCITY] = ANGLE_MOD(feat[F_ANGLE_CHANGE] - 
						 features_prev_frame[F_ANGLE_CHANGE])*dt_inv;
    feat[F_LENGTH_VELOCITY] = (feat[F_LENGTH] - features_prev_frame[F_LENGTH])*dt_inv;
    if(pt_prev)
      cost += COST(F_SEGMENT_ANGULAR_VELOCITY) + COST(F_LENGTH_VELOCITY);
  }

  return cost;
}

SpinePointTransitions *compute_plausible_transitions(int ind, SpinePointLocation *legal_points, 
						     int num_legal, double *contour, int num_contour_points, 
						     FitParams *params) {
  SpinePointTransitions *retval = (SpinePointTransitions*)malloc(sizeof(SpinePointTransitions)+
								 F_NUM_UNARY_FEATURES*sizeof(double));
  int num_alloc = 0, i, *ind2;
  SpinePointLocation *pt = &legal_points[ind];
  double *features;

  retval->num = 0;
  retval->pt = pt;
  retval->ind = ind;
  retval->unary_features = (double*)(retval+1);
  compute_unary_features(pt, retval->unary_features, contour, num_contour_points, params);
								 
  for(i = 0; i < num_legal; i++) {
    if(i == ind)
      continue;
    if(retval->num >= num_alloc) {
      num_alloc += 16;
      retval = (SpinePointTransitions*)realloc(retval, sizeof(SpinePointTransitions) + 
					       F_NUM_UNARY_FEATURES*sizeof(double) +
					       num_alloc*PAIR_SIZE);
      retval->unary_features = (double*)(retval+1);
      retval->pairwise_features = (unsigned char*)(retval->unary_features + F_NUM_UNARY_FEATURES);
    }
    ind2 = ((int*)(retval->pairwise_features+retval->num*PAIR_SIZE));
    *ind2 = i;
    features = ((double*)(ind2+1))-F_NUM_UNARY_FEATURES-1;
    if(compute_pairwise_features(pt, &legal_points[i], features, contour, num_contour_points, params))
      retval->num++;
  }

  return retval;
}

SpinePointTransitions **compute_pairwise_transitions(SpinePointLocation *legal_points, 
						     int num_legal, double *contour, int num_contour_points, 
						     FitParams *params) {
  SpinePointTransitions **retval = (SpinePointTransitions**)malloc(sizeof(SpinePointTransitions*)*num_legal);
  int i;

  for(i = 0; i < num_legal; i++) {
    retval[i] = compute_plausible_transitions(i, legal_points, num_legal, contour, num_contour_points, params);
  }
  return retval;
}


void free_pairwise_spine_transitions(SpinePointTransitions **t, int num) {
  int n;

  for(n = 0; n < num; n++) 
    free(t[n]);
  free(t);
}

void print_features(double *feat, int s, bool has_prev, double *mu, double *var, FitParams *params) {
  const char *names[F_NUM_FEATURES] = { "width", "head_length", "tail_length", "head_symmetry", "tail_symmetry", 
					"head_geo_symmetry", "tail_geo_symmetry", "head_percent_outliers", 
					"tail_percent_outliers", "head_geometry", "tail_geometry", "width_velocity", 
					"forward_velocity", "sideways_velocity",
					"angular_velocity", "", "percent_outliers", "width_change", "angle_change", 
					"length", "left_parallel", "left_orthogonal", "right_parallel", "right_orthogonal", 
					"symmetry", "geometry1", "geometry2", "segment_angular_velocity", "length_velocity" };
  bool is_head[F_NUM_FEATURES] =     { 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
  bool is_tail[F_NUM_FEATURES] =     { 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
  bool is_temporal[F_NUM_FEATURES] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1 };
  bool is_invalid[F_NUM_FEATURES]  = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
  bool is_pairwise[F_NUM_FEATURES] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };

  for(int i = 0; i < F_NUM_FEATURES; i++) {
    if(!is_invalid[i] && (has_prev || !is_temporal[i]) && (s == 0 || !is_head[i]) && (s == params->num_spine_lines || !is_tail[i]) && 
       (s > 0 || !is_pairwise[i]))
      fprintf(stderr, "  %d.%s: val=%f, cost=%f, mu=%f var=%f\n", s, names[i], (float)feat[i], (float)COST(i), (float)mu[i], (float)var[i]);
  }
}

void compute_features(Blob *b, Blob *prev, double dt_inv, FitParams *params, int debug) {
  double *contour = b->fixed_contour ? b->fixed_contour : b->contour;
  int num_pts =  b->fixed_contour ? b->fixed_num_pts : b->num_pts;
  double cost = 0, c;
  int i, n;

  if(b->features)
    free(b->features);
  b->features = (double*)malloc(2*sizeof(double)*((2*F_NUM_FEATURES)*my_max(params->num_spine_lines+1,b->num_model_pts)+F_NUM_GLOBAL_FEATURES));

  double *sqr_features = SQR_FEATURES(b);
  double *normalized_features = NORMALIZED_FEATURES(b);
  double *normalized_sqr_features = NORMALIZED_SQR_FEATURES(b);

  g_debug = 1;
  for(n = 0; n < b->num_model_pts; n++) {
    g_debug = n == 2;
    compute_unary_features(&b->model_pts[n], b->features+n*F_NUM_FEATURES, contour, num_pts, params);
    c = eval_unary_cost(&b->model_pts[n], n, b->features+n*F_NUM_FEATURES,  
			prev != NULL ? &prev->model_pts[n] : NULL, dt_inv, params->mu, params->var, params);
    if(n > 0) {
      compute_pairwise_features(&b->model_pts[n], &b->model_pts[n-1], b->features+n*F_NUM_FEATURES, 
				contour, num_pts, params);
      c += eval_pairwise_cost(&b->model_pts[n], n, b->features+n*F_NUM_FEATURES,
			      &b->model_pts[n-1], prev != NULL && prev->features ? prev->features+n*F_NUM_FEATURES : NULL, 
			      dt_inv, params->mu, params->var, params);
    }
    if(debug) {
      print_features(b->features+n*F_NUM_FEATURES, n, prev != NULL, params->mu, params->var, params);
      fprintf(stderr, "  cost = %f\n", (float)c);
    }
    cost += c;

    for(i = 0; i < F_NUM_FEATURES; i++) {
      sqr_features[n*F_NUM_FEATURES+i] = SQR(b->features[n*F_NUM_FEATURES+i]);
      normalized_features[n*F_NUM_FEATURES+i] = b->features[n*F_NUM_FEATURES+i] / sqrt(params->var[i]);
      normalized_sqr_features[n*F_NUM_FEATURES+i] = sqr_features[n*F_NUM_FEATURES+i] / params->var[i];
    }
  }

  g_debug = 0;
  if(debug) fprintf(stderr, "total cost = %f\n", (float)cost);
}

void compute_all_features(BlobSequence *s, FitParams *params, int debug) {
  Blob *prev = NULL;
  int t;

  for(t = 0; t < s->num_frames; t++) {
    if(!s->frames[t].features) 
      s->frames[t].features = (double*)malloc(2*sizeof(double)*((2*F_NUM_FEATURES)*(params->num_spine_lines+1)+F_NUM_GLOBAL_FEATURES));
    compute_features(&s->frames[t],  prev, prev ? (1/(s->frames[t].frame_time-prev->frame_time)) : 0, params, debug);
    prev = s->frames[t].is_good ? & s->frames[t] : NULL;
  }
}

void compute_all_global_features(BlobSequence *s, FitParams *params, int debug) {
  int n, i, j, t, prev_n, ori, p1, p2;
  Blob *b;
  double *feat, *global_feat, *norm_global_feat;
  double vx, vy, len, m, curr_len, prev_curr_len, w, x, y, width, d, target_len, x1, y1, x2, y2;
  double *contour;
  int num_pts;
  double dists[1000];
  double sum_length = 0, max_length = 0;

  for(t = 0; t < s->num_frames; t++) {
    b = &s->frames[t];
    feat = b->features;
    if(!feat || !b->is_good)
      continue;
    global_feat = GLOBAL_FEATURES(b);
    norm_global_feat = NORMALIZED_GLOBAL_FEATURES(b);
    contour = b->fixed_contour ? b->fixed_contour : b->contour;
    num_pts =  b->fixed_contour ? b->fixed_num_pts : b->num_pts;

    // Compute the total length of the spine
    len = b->model_pts[0].d4 + b->model_pts[b->num_model_pts-1].d3;
    for(n = 1; n < b->num_model_pts; n++) {
      dists[n] = point_dist(b->model_pts[n-1].x, b->model_pts[n-1].y, b->model_pts[n].x, b->model_pts[n].y);
      len += dists[n];
    }
    global_feat[F_G_LENGTH] = len;
    if(len > max_length) 
      max_length = len;
    sum_length += len;

    global_feat[F_G_AVE_WIDTH] = 0;
    global_feat[F_G_MAX_WIDTH] = 0;
    global_feat[F_G_AVE_X] = 0;
    global_feat[F_G_AVE_Y] = 0;
    prev_n = n = 0;
    prev_curr_len = curr_len = b->model_pts[0].d4;
    target_len = len/(F_NUM_GLOBAL_POINTS+1);
    for(i = 0; i < F_NUM_GLOBAL_POINTS; i++) {
      // Compute the point along the spine at distance target_len from the head
      while(target_len > curr_len && n < b->num_model_pts) {
	prev_n = n;
	prev_curr_len = curr_len;
	n++;
	assert(n < b->num_model_pts);
	curr_len += dists[n];
	if(n == b->num_model_pts-1)
	  curr_len += b->model_pts[b->num_model_pts-1].d3;
      }
      w = (curr_len-prev_curr_len) > 0 ? (target_len-prev_curr_len) / (curr_len-prev_curr_len) : 1;
      if(w < 0 || w > 1) w = 1;
      x = (1-w)*b->model_pts[prev_n].x + w*b->model_pts[n].x;
      y = (1-w)*b->model_pts[prev_n].y + w*b->model_pts[n].y;
      
      // Compute the width of the cross-section at (x,y)
      width = INFINITY;
      ori = -1;
      for(j = 0; j < params->num_orientations; j++) {
	p1 = nearest_to_line(x, y, COS_ORI(j), SIN_ORI(j), 1, contour, num_pts);
	p2 = nearest_to_line(x, y, COS_ORI(j), SIN_ORI(j), -1, contour, num_pts);
	if(p1 < 0) { x1 = x; y1 = y; }
	else { x1 = contour[2*p1]; y1 = contour[2*p1+1]; }
	if(p2 < 0) { x2 = x; y2 = y; }
	else { x2 = contour[2*p2]; y2 = contour[2*p2+1]; }
	d = point_dist(x1, y1, x2, y2);
	if(d < width) {
	  width = d;
	  ori = j;
	}
      }
      assert(ori >= 0);

      global_feat[F_G_X1+i] = x;
      global_feat[F_G_Y1+i] = y;
      global_feat[F_G_WIDTH1+i] = width;
      global_feat[F_G_AVE_WIDTH] += width/F_NUM_GLOBAL_POINTS;
      global_feat[F_G_AVE_X] += x/F_NUM_GLOBAL_POINTS;
      global_feat[F_G_AVE_Y] += y/F_NUM_GLOBAL_POINTS;
      if(width > global_feat[F_G_MAX_WIDTH])
	global_feat[F_G_MAX_WIDTH] = width;

      target_len += len/(F_NUM_GLOBAL_POINTS+1);
    }
    global_feat[F_G_WIDTH_LENGTH_RATIO] = global_feat[F_G_MAX_WIDTH] / global_feat[F_G_LENGTH];

    
    // Instead of measuring the angle or curvature of the worm directly, we measure the (normalized) distance of each
    // spline point from  the line connecting the head and tail points, with the intuition that internal
    // points of turned larva will be on one side of the line, with distance proportional to the strength of turn.
    // Sigmoid shaped larva will be on one side of the line for the tail and the other side for the head.  
    global_feat[F_G_MAX_LEFT_CURVATURE] = 0;
    global_feat[F_G_MAX_RIGHT_CURVATURE] = 0;
    vx = b->model_pts[0].x - b->model_pts[b->num_model_pts-1].x;
    vy = b->model_pts[0].y - b->model_pts[b->num_model_pts-1].y;
    m = sqrt(SQR(vx)+SQR(vy));
    vx /= SQR(m);  vy /= SQR(m);
    for(i = 0; i < F_NUM_GLOBAL_POINTS; i++) {
      global_feat[F_G_CURVE1+i] = (vy*(global_feat[F_G_X1+i]-b->model_pts[b->num_model_pts-1].x) - 
				   vx*(global_feat[F_G_Y1+i]-b->model_pts[b->num_model_pts-1].y));
      if(-global_feat[F_G_CURVE1+i] > global_feat[F_G_MAX_LEFT_CURVATURE])
	global_feat[F_G_MAX_LEFT_CURVATURE] = -global_feat[F_G_CURVE1+i];
      if(global_feat[F_G_CURVE1+i] > global_feat[F_G_MAX_RIGHT_CURVATURE])
	global_feat[F_G_MAX_RIGHT_CURVATURE] = global_feat[F_G_CURVE1+i];
    }
    
  }

  sum_length /= s->num_frames;

  for(t = 0; t < s->num_frames; t++) 
    compute_global_features(s, t, params, sum_length, debug);
}


void compute_global_features(BlobSequence *s, int t, FitParams *params, double ave_length, int debug) {
  int t1, t2, i;
  double vx, vy, m, dt;
  double *global_feat, *global_feat1, *global_feat2, *norm_global_feat;
  Blob *b = &s->frames[t];

  // Compute the velocity averaged over a time window of params->speed_estimate_window seconds
  if(!s->frames[t].is_good)
    return;
  for(t1 = t; t1 >= 0 && s->frames[t1].is_good && s->frames[t1].features &&
	(s->frames[t].frame_time-s->frames[t1].frame_time) < params->speed_estimate_window/2; t1--);
  for(t2 = t; t2 < s->num_frames && s->frames[t2].is_good && s->frames[t2].features &&
	(s->frames[t2].frame_time-s->frames[t].frame_time) < params->speed_estimate_window/2; t2++);
  t1++; t2--;
  
  global_feat1 = GLOBAL_FEATURES(&s->frames[t1]);
  global_feat2 = GLOBAL_FEATURES(&s->frames[t2]);
  global_feat = GLOBAL_FEATURES(b);
  norm_global_feat = NORMALIZED_GLOBAL_FEATURES(b);

  vx = b->model_pts[0].x - b->model_pts[b->num_model_pts-1].x;
  vy = b->model_pts[0].y - b->model_pts[b->num_model_pts-1].y;
  m = sqrt(SQR(vx)+SQR(vy));
  if(m > 0) { vx /= m;  vy /= m; }
  else {
    b->is_good = 0;
    return; 
  }

  // The component of the velocity projected along the larva's tail to head axis
  dt = (s->frames[t2].frame_time-s->frames[t1].frame_time);
  global_feat[F_G_FORWARD_VELOCITY] = (vx*(global_feat2[F_G_AVE_X]-global_feat1[F_G_AVE_X]) + 
				       vy*(global_feat2[F_G_AVE_Y]-global_feat1[F_G_AVE_Y]))/ave_length / dt;
    
  // The component of the velocity orthogonal to the tail to head axis
  global_feat[F_G_SIDEWAYS_VELOCITY] = (vy*(global_feat2[F_G_AVE_X]-global_feat1[F_G_AVE_X]) - 
					vx*(global_feat2[F_G_AVE_Y]-global_feat1[F_G_AVE_Y]))/ave_length / dt;

  memcpy(norm_global_feat, global_feat, sizeof(double)*F_NUM_GLOBAL_FEATURES);
  norm_global_feat[F_G_FORWARD_VELOCITY] /= params->expected_length;
  norm_global_feat[F_G_SIDEWAYS_VELOCITY] /= params->expected_length;
  norm_global_feat[F_G_MAX_WIDTH] /= params->expected_length;
  norm_global_feat[F_G_AVE_WIDTH] /= params->expected_length;
  for(i = 0; i < F_NUM_GLOBAL_POINTS; i++) 
    norm_global_feat[F_G_WIDTH1+i] /= params->expected_length;

  if(debug) {
    fprintf(stderr, "%d:\n", t);
    print_global_features(global_feat);
  }
}

const char *g_global_feature_names[F_NUM_GOOD_GLOBAL_FEATURES] = { "width_length_ratio", "forward_velocity", "sideways_velocity", "max_left_curvature", 
								   "max_right_curvature", "curve1", "curve2", "curve3", "curve4", "curve5", "curve6" };
void print_global_features(double *feat) {
  int i;

  for(i = 0; i < F_NUM_GOOD_GLOBAL_FEATURES; i++)
    fprintf(stderr, "  %s: %f\n", g_global_feature_names[i], (float)feat[i]);
}
