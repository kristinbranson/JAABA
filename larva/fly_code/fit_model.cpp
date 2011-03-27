
#include "fit_model.h"
#include "spine_features.h"

#include <cv.h>
#include <highgui.h>
#include <ml.h>  

#ifdef WIN32

#pragma comment(lib, "cv210")
#pragma comment(lib, "cxcore210")
#pragma comment(lib, "highgui210")
#pragma comment(lib, "ml210")
#endif

#define SCALE 20

extern int g_debug;


#define SHOW_FIXED_CONTOUR 0




// Cached innformation corresponding to one Blob frame used for dynamic programming
typedef struct {
  Blob *blob;
  FitParams *params;
  double *contour;
  int num_pts;

  IplImage *img;     // a binary silhouette image of the blob 
  double offset[2];

  // The set of points that are legal states (used because many spine points have the
  // same set of legal states)
  SpinePointLocation *legal_points;
  int num_legal_points;

  // a num_legal_points array of legal transitions between states.  This is a sparse version
  // of a num_legal_pointsXnum_legal_points array caching all features for evaluating 
  // the transition costs between states
  struct _SpinePointTransitions **transitions;

  double **table_costs;  // a num_skeleton_points X num_states array that caches the minimum
                         //   achievable cost to reach a particular state
  int **table_states;    // a num_skeleton_points X num_states array containing indices into
                         //   states[n] corresponding to the state achieving minimal cost

} DPCache;

DPCache *new_dp_cache(Blob *b, FitParams *params);
void free_dp_cache(DPCache *d);


// Simple function to fit the entire sequence frame by frame, either using temporal info greedily or not using temporal info at all
double fit_model_simple(BlobSequence *s, FitParams *params, double *costs, bool use_temporal_info, bool free_memory,
			int debug, void (*on_frame)(int t), int *keep_going) {
  int t;
  double cost, total_cost = 0;
  Blob *prev = NULL;

  for(t = 0; t < s->num_frames && (!keep_going || *keep_going); t++) {
    if(debug) fprintf(stderr, "%d: ", t);

    if(!s->frames[t].model_pts || !(s->frames[t].is_manual & 1)) {
      cost = fit_model_one_frame(&s->frames[t], prev, params, 0, NULL, debug);
      if(costs) costs[t] = cost;
      if(s->frames[t].is_good) {
	total_cost += cost;
	compute_features(&s->frames[t], NULL, 0, params, 0);
      }
    } else {
      s->frames[t].is_good = 1;
      if(costs) costs[t] = 0;
    }

    if(on_frame)
      (*on_frame)(t);

    if(free_memory) {
      free_dp_cache((DPCache*)s->frames[t].data);
      s->frames[t].data = NULL;
    }

    if(use_temporal_info && s->frames[t].is_good) 
      prev = &s->frames[t];
    else 
      prev = NULL;
  }

  return total_cost;
}

// Two-pass greedy fitting procedure.  First fit spine on each frame independently.  Then greedily propagate temporal information
// frame-to-frame, at each step propagating temporal information from the frame with the highest confidence
double fit_model_temporal_greedy_two_pass(BlobSequence *s, FitParams *params, int debug, void (*on_frame)(int t), int *keep_going) {
  double *costs = (double*)malloc(sizeof(double)*s->num_frames);
  bool *used = (bool*)malloc(2*sizeof(bool)*s->num_frames);
  bool *fixed = used + s->num_frames;
  double cost, total_cost = 0, best_cost;
  int t, t2, best;

  memset(used, 0, 2*sizeof(bool)*s->num_frames);

  // First fit without using temporal information
  fit_model_simple(s, params, costs, false, false, debug, on_frame, keep_going);

  for(t = 0; t < s->num_frames && (!keep_going || *keep_going); t++) {

    // Pick the next frame with lowest cost
    best = -1;
    best_cost = INFINITY;
    for(t2 = 0; t2 < s->num_frames; t2++) {
      if(!used[t2]) {
	if(s->frames[t2].is_manual&1) {
	  best = t2;
	  break;
	}
	if(costs[t2] < best_cost) {
	  best_cost = costs[t2];
	  best = t2;
	}
      }
    }

    if(!fixed[best]) {
      total_cost += costs[best];
      fixed[best] = true;
    }
    used[best] = true;
    if(on_frame)
      (*on_frame)(best);

    if(s->frames[best].is_good) {
      if(best > 0 && !(s->frames[best-1].is_manual&1) && !fixed[best-1]) {
	// Propagate info to the frame to the left
	if(debug) fprintf(stderr, "2nd pass, frame %d: ", best-1);
	cost = fit_model_one_frame(&s->frames[best-1], &s->frames[best], params, 0, NULL, debug);
	total_cost += cost;
	fixed[best-1] = true;
	costs[best-1] = (cost + costs[best])/2;
	compute_features(&s->frames[best-1], &s->frames[best], 1.0/(s->frames[best].frame_time-s->frames[best-1].frame_time), params, 0);
      }

      if(best < s->num_frames-1 && !(s->frames[best+1].is_manual&1) && !fixed[best+1]) {
	// Propagate info to the frame to the right
	if(debug) fprintf(stderr, "2nd pass, frame %d: ", best+1);
	cost = fit_model_one_frame(&s->frames[best+1], &s->frames[best], params, 0, NULL, debug);
	total_cost += cost;
	fixed[best+1] = true;
	costs[best+1] = (cost + costs[best])/2;
	compute_features(&s->frames[best+1], &s->frames[best], 1.0/(s->frames[best+1].frame_time-s->frames[best].frame_time), params, 0);
      }
    }
  }

  free(costs);
  free(used);

  return total_cost;
}


void fit_model_finish(BlobSequence *s, FitParams *params, int can_flip_spine) {
  double total_forward_motion, sum_dist, sum_flip_dist;
  int t, n;
  Blob *prev = NULL;

  // Free cache tables if necessary
  for(t = 0; t < s->num_frames; t++) {
    if(s->frames[t].data) {
      free_dp_cache((DPCache*)s->frames[t].data);
      s->frames[t].data = NULL;
    }
  }

  if(can_flip_spine > 1) {
    // Maintain head/tail consistency across all pairs of frames.  Iterate through the full sequence
    // and flip the spine in frame t if flipping brings the model points closer to frame t-1
    for(t = 0; t < s->num_frames; t++) {
      if(s->frames[t].is_good) {
	sum_dist = sum_flip_dist = 0;
	if(prev) {
	  for(n = 0; n < s->frames[t].num_model_pts; n++) {
	    sum_dist += point_dist(s->frames[t].model_pts[n].x, s->frames[t].model_pts[n].y, 
				   prev->model_pts[n].x, prev->model_pts[n].y);
	    sum_flip_dist += point_dist(s->frames[t].model_pts[s->frames[t].num_model_pts-1-n].x, 
					s->frames[t].model_pts[s->frames[t].num_model_pts-1-n].y, 
					prev->model_pts[n].x, prev->model_pts[n].y);
	  }
	  if(sum_flip_dist < sum_dist) {
	    flip_spine(&s->frames[t], params);
	    compute_features(&s->frames[t], prev, prev ? (1/(s->frames[t].frame_time-prev->frame_time)) : 0, params, 0);
	  }
	}
	prev = &s->frames[t];
      }
    }
  }

  compute_all_global_features(s, params, 0);

  if(can_flip_spine > 0) {
    // Classify the head/tail of the sequence based on the total distance travelled
    total_forward_motion = 0;
    for(t = 1; t < s->num_frames; t++) {
      if(s->frames[t].features && s->frames[t].is_good)
	total_forward_motion += GLOBAL_FEATURES(&s->frames[t])[F_G_FORWARD_VELOCITY]*(s->frames[t].frame_time-s->frames[t-1].frame_time);
      //assert(!isnan(total_forward_motion));
    }
    if(total_forward_motion < 0)
      flip_spine_sequence(s, params);
  }
}

double fit_model(BlobSequence *s, FitParams *params, int debug, void (*on_frame)(int t), int *keep_going) {
  int can_flip_spine = 2;
  double retval = 0;

  switch(params->fit_method) {
  case FIT_FRAME_BY_FRAME:
  case FIT_GREEDY_ONE_PASS:
    retval = fit_model_simple(s, params, NULL, params->fit_method == FIT_GREEDY_ONE_PASS, true,
			      debug, on_frame, keep_going);
    break;
  case FIT_GREEDY_TWO_PASS:
    retval = fit_model_temporal_greedy_two_pass(s, params, debug, on_frame, keep_going);
    break;
  case FIT_GLOBALLY_OPTIMAL:
    fprintf(stderr, "Globally optimal spine fit method not implemented\n");
    can_flip_spine = 1;
    assert(0);
    break;
  }

  fit_model_finish(s, params, can_flip_spine);

  return retval;
}

  
// Fit a spine model to a single frame, assuming the location of the previous
// time frame is known
double fit_model_one_frame(Blob *b, Blob *prev, FitParams *params, int free_cache, Blob *result, int debug) {
  int n, i, j, k, best;
  double unary_cost, pairwise_cost, best_cost;
  DPCache *d;
  int ds = PAIR_SIZE, df = ((int)sizeof(int)) - ((int)sizeof(double))*(F_NUM_UNARY_FEATURES+1);
  double *features;
  unsigned char *ptr;
  double dt_inv = prev ? 1/(b->frame_time-prev->frame_time) : 0;
  Blob *b_orig = b;
  Blob tmp = *b;
  b = &tmp;
  
  if(!b->data) {
    d = new_dp_cache(b, params);
    d->transitions = compute_pairwise_transitions(d->legal_points, d->num_legal_points, d->contour,
						  d->num_pts, params);
  } else 
    d = (DPCache*)b->data;

  // Use dynamic programming to compute the optimal spine points
  for(n = 0; n < params->num_spine_lines+1; n++) {			       
    // Find the optimal cost for every sub-solution through spine points 1...n for
    // every possible legal hidden state
    for(j = 0; j < d->num_legal_points; j++) {
      unary_cost = eval_unary_cost(d->transitions[j]->pt, n, d->transitions[j]->unary_features, 
				   prev != NULL ? &prev->model_pts[n] : NULL, dt_inv,
				   params->mu, params->var, params);

      // Compute the position of the previous spine part that minimizes
      // the cost needed to get to this state of the current spine part
      if(n > 0 && unary_cost < INFINITY) {
	best = -1;
	best_cost = INFINITY;
	for(i = 0, ptr = d->transitions[j]->pairwise_features; i < d->transitions[j]->num; i++, ptr += ds) {
	  k = *(int*)ptr;
	  if(d->table_costs[n-1][k] != INFINITY) {
	    features = (double*)(ptr + df);
	    pairwise_cost = d->table_costs[n-1][k] + eval_pairwise_cost(d->transitions[j]->pt, n, features,
						&d->legal_points[k], prev != NULL && prev->features != NULL ? 
									prev->features+n*F_NUM_FEATURES : NULL, 
									dt_inv, params->mu, params->var, params);
	    if(pairwise_cost < best_cost) {
	      best_cost = pairwise_cost;
	      best = k;
	    } 
	  }
	}
      } else { 
	best_cost = n > 0 ? INFINITY : 0;
	best = -1;
      }
      // Store the cost and hidden state pointers that yield optimal (note: this is
      // not implemented in a memory efficient way)
      d->table_costs[n][j] = unary_cost + best_cost;
      d->table_states[n][j] = best;
    }
  }

  // Allocate data to save the fit model
  b->data = d;
  b->model_pts = (SpinePointLocation*)malloc(sizeof(SpinePointLocation)*(params->num_spine_lines+1));
  b->num_model_pts = params->num_spine_lines+1;

  // Select the hidden state for the last spine point with lowest cost, then
  // trace back the cache tables to extract the corresponding solution for all
  // other points
  best = -1;
  best_cost = INFINITY;
  for(j = 0; j < d->num_legal_points; j++) {
    if(d->table_costs[params->num_spine_lines][j] < best_cost) {
      best_cost = d->table_costs[params->num_spine_lines][j];
      best = j;
    }
  }
  n = params->num_spine_lines;
  while(n >= 0 && best >= 0) { 
    //fprintf(stderr, "line %d: (%f,%f) %f %f, %d\n", n, (float)d->legal_points[best].x, (float)d->legal_points[best].y, 
    //    (float)ORI(d->legal_points[best].orientation), (float)d->legal_points[best].width, d->table_states[n][best]);
    b->model_pts[n] = d->legal_points[best];

    // Backtrack to read out the rest of the optimal solution
    best = d->table_states[n][best];
    n--;
  }

  // Marks a frame as bad if no solution was found.
  // TODO: threshold on the cost?
  b->is_good = n <= 0 && best_cost < INFINITY;  
  if(!b->is_good) {
    free(b->model_pts);
    b->model_pts = NULL;
    b->num_model_pts = 0;
  }

  if(debug) fprintf(stderr, "  %f\n", (float)best_cost);

  if(free_cache) {
    free_dp_cache(d);
    b->data = NULL;
    d = NULL;
  }

  if(!result) {
    if(d) d->blob = b_orig;
    if(b_orig->model_pts)
      free(b_orig->model_pts);
    *b_orig = tmp;
  } else {
    if(d) d->blob = result;
    *result = tmp;
  }

  return best_cost;
}




// Convert contour points to a binary mask image and return it.
// Also stores the x and y offsets used to convert the contour 
// to a cropped image of minimum size
IplImage *binary_mask(double *contour, int num_pts, double *offsets, FitParams *params) {
  double min_x = INFINITY, max_x = -INFINITY, min_y = INFINITY, max_y = -INFINITY;
  int i, w, h;
  IplImage *blob;
  CvPoint *pts = (CvPoint*)malloc(sizeof(CvPoint)*num_pts);

  // Find the bounding box around the set of contour points
  for(i = 0; i < num_pts; i++) {
    if(contour[2*i] < min_x) min_x = contour[2*i];
    if(contour[2*i] > max_x) max_x = contour[2*i];
    if(contour[2*i+1] < min_y) min_y = contour[2*i+1];
    if(contour[2*i+1] > max_y) max_y = contour[2*i+1];
  }
  offsets[0] = -min_x+1/params->world_to_pixel_scale;
  offsets[1] = -min_y+1/params->world_to_pixel_scale;
  w = ceil((max_x-min_x)*params->world_to_pixel_scale)+3;
  h = ceil((max_y-min_y)*params->world_to_pixel_scale)+3;
  
  // Create a binary image of the blob by drawing a closed polygon
  blob = cvCreateImage(cvSize(w, h), IPL_DEPTH_8U, 1);
  cvZero(blob);
  for(i = 0; i < num_pts; i++) {
    pts[i].x = (int)(.5+(contour[2*i]+offsets[0])*params->world_to_pixel_scale);
    pts[i].y = (int)(.5+(contour[2*i+1]+offsets[1])*params->world_to_pixel_scale);
  }
  cvFillPoly(blob,&pts,&num_pts,1,cvScalar(255,255,255));
  
  free(pts);

  return blob;
}


// Routine to remove topological components from a contour if some component has a width below some
// threshold (params->blob_min_pixel_width) at any point along the contour
int fix_contour(double *contour, int num_pts, double *offsets, FitParams *params, double **new_contour, 
		int *new_num_pts, IplImage **new_mask) {
  IplImage *mask = binary_mask(contour, num_pts, offsets, params);
  IplImage *reduced_mask = cvCloneImage(mask);
  IplImage *gb = cvCloneImage(mask);
  CvMemStorage* storage = cvCreateMemStorage(0), *storage2;
  CvSeq *cont = NULL, *ptr, *biggest = NULL;
  int num_biggest = 0;
  CvPoint *pt;
  int retval = 0, i;

  // Shrink the mask by params->blob_min_pixel_width pixels, and select the largest remaining 
  // connected component
  cvSmooth(mask, gb, CV_GAUSSIAN, params->blob_min_pixel_width*2+1);
  cvThreshold(gb, reduced_mask, 254, 255, CV_THRESH_BINARY);
  cvFindContours(reduced_mask, storage, &cont, sizeof(CvContour), CV_RETR_CCOMP, CV_CHAIN_APPROX_NONE );
  for(ptr = cont; ptr != NULL; ptr = ptr->h_next ) {
    if(ptr->total > num_biggest) {
      biggest = ptr;
      num_biggest = ptr->total;
    }
  }

  if(biggest) {
    // Render the largest component and grow the mask back to its original size
    storage2 = cvCreateMemStorage(0);
    //cvSaveImage("reduced.png", reduced_mask);
    cvZero(reduced_mask);
    cvDrawContours(reduced_mask, biggest, CV_RGB(255,255,255), CV_RGB(255,255,255), -1, CV_FILLED);
    cvSmooth(reduced_mask, gb, CV_GAUSSIAN, params->blob_min_pixel_width*2+1);
    cvThreshold(gb, gb, 1, 1, CV_THRESH_BINARY);
    cvMul(mask, gb, reduced_mask);
    //cvSaveImage("gb.png", gb);
    //cvSaveImage("expanded.png", reduced_mask);

    // Extract the new contour and convert back to world coordinates
    cvFindContours(reduced_mask, storage2, &cont, sizeof(CvContour), CV_RETR_CCOMP, CV_CHAIN_APPROX_NONE);    
    *new_contour = (double*)malloc(2*sizeof(double)*cont->total);
    for(i = 0; i < cont->total; i++) {
      pt = CV_GET_SEQ_ELEM(CvPoint, cont, i); 
      (*new_contour)[2*i] = P2W(pt->x,0);
      (*new_contour)[2*i+1] = P2W(pt->y,1);
    }
    *new_num_pts = cont->total;
    *new_mask = reduced_mask;
    retval = true;
    cvReleaseMemStorage(&storage2);
  } else {
    cvReleaseImage(&reduced_mask);
    *new_mask = NULL;
    *new_num_pts = 0;
    *new_contour = NULL;
  }

  cvReleaseMemStorage(&storage);
  cvReleaseImage(&mask);
  cvReleaseImage(&gb);
  
  return retval;
}

// Compute spine point parameters that should be fully determinable given (x,y,theta)
void compute_intersect_points(SpinePointLocation *pt, double *contour, int num_pts, 
			      FitParams *params) {
  int p1, p2, p3, p4;

  pt->p1 = p1 = nearest_to_line(pt->x, pt->y, COS_ORI(pt->orientation),
				SIN_ORI(pt->orientation), 1, contour, num_pts);
  pt->p2 = p2 = nearest_to_line(pt->x, pt->y, COS_ORI(pt->orientation),
				SIN_ORI(pt->orientation), -1, contour, num_pts);
  pt->p3 = p3 = nearest_to_line(pt->x, pt->y, -SIN_ORI(pt->orientation),
				COS_ORI(pt->orientation), -1, contour, num_pts);
  pt->p4 = p4 = nearest_to_line(pt->x, pt->y, -SIN_ORI(pt->orientation),
				COS_ORI(pt->orientation), 1, contour, num_pts);
  if(p1 >= 0 && p2 >= 0)
    pt->width = point_dist(contour[2*p1], contour[2*p1+1], 
			   contour[2*p2], contour[2*p2+1]);
  else
    pt->width = 0;

  if(p3 >= 0)
    pt->d3 = point_dist(pt->x, pt->y, contour[2*p3], contour[2*p3+1]);
  else
    pt->d3 = 0;

  if(p4 >= 0)
    pt->d4 = point_dist(pt->x, pt->y, contour[2*p4], contour[2*p4+1]);
  else
    pt->d4 = 0;
}


void compute_deterministic_spine_attributes(SpinePointLocation *pt, SpinePointLocation *prev, 
					    SpinePointLocation *next, double *contour, int num_pts, 
					    FitParams *params, int find_orientation) {
  int i;
  double width, dy, dx;
  
  if(find_orientation) {
    if(prev) { dy = pt->y-prev->y;  dx = pt->x-prev->x; }
    else { dy = next->y-pt->y;  dx = next->x-pt->x; }
  
    pt->width = 0;
    pt->orientation = 0;
    for(i = 0; i < params->num_orientations; i++) {
      if(dy*COS_ORI(i) - dx*SIN_ORI(i) < 0) {
	pt->p1 = nearest_to_line(pt->x, pt->y, COS_ORI(i), SIN_ORI(i), 1, contour, num_pts);
	pt->p2 = nearest_to_line(pt->x, pt->y, COS_ORI(i), SIN_ORI(i), -1, contour, num_pts);
	if(pt->p1 >= 0 && pt->p2 >= 0) {
	  width = point_dist(contour[2*pt->p1], contour[2*pt->p1+1], 
			     contour[2*pt->p2], contour[2*pt->p2+1]);
	  if(width > 0 && (pt->width == 0 || width < pt->width)) {
	    pt->width = width;
	    pt->orientation = i;
	  }
	}
      }
    }
  }

  compute_intersect_points(pt, contour, num_pts, params);
  
}

void compute_deterministic_spine_attributes_blob(Blob *blob, FitParams *params, int find_orientations) {
  int i;
  double *contour = blob->fixed_contour ? blob->fixed_contour : blob->contour;
  int num_pts =  blob->fixed_contour ? blob->fixed_num_pts : blob->num_pts;

  if(blob->model_pts) 
    for(i = 0; i < blob->num_model_pts; i++) {
      compute_deterministic_spine_attributes(&blob->model_pts[i], i > 0 ? &blob->model_pts[i-1] : NULL,
					     i < blob->num_model_pts-1 ? &blob->model_pts[i+1] : NULL,
					     contour, num_pts, params, find_orientations);
    }
}

void compute_deterministic_spine_attributes_blob_sequence(BlobSequence *s, FitParams *params, int find_orientations) {
  int i;
  
  for(i = 0; i < s->num_frames; i++)
    compute_deterministic_spine_attributes_blob(&s->frames[i], params, find_orientations);
}

void flip_spine(Blob *blob, FitParams *params) {
  SpinePointLocation tmp;
  int i;

  for(i = 0; i < blob->num_model_pts/2; i++) {
    tmp = blob->model_pts[i];
    blob->model_pts[i] = blob->model_pts[blob->num_model_pts-1-i];
    blob->model_pts[blob->num_model_pts-1-i] = tmp;
  }
  for(i = 0; i < blob->num_model_pts; i++) 
    blob->model_pts[i].orientation = (blob->model_pts[i].orientation + params->num_orientations/2) % params->num_orientations;
  
  compute_deterministic_spine_attributes_blob(blob, params, 1);
}

void flip_spine_sequence(BlobSequence *blob, FitParams *params) {
  int i;

  for(i = 0; i < blob->num_frames; i++) {
    flip_spine(&blob->frames[i], params);

    if(blob->frames[i].features)
      compute_features(&blob->frames[i],  i > 0 ? &blob->frames[i] : NULL, 
		       i ? (1/(blob->frames[i].frame_time-blob->frames[i-1].frame_time)) : 0, params, 0);
  }
}

// Helper function for get_legal_internal_points(). Go from a pixel in a rotated blob image
// to the corresponding point in world coordinates
void rotated_image_coords_to_world_coords(double x_img, double y_img, int ori,
                                          double mean_x, double mean_y, double *rot_offsets,
                                          FitParams *params, double *x, double *y) {
  // Convert from pixel coordinates to (rotated) world coordinates
  double x_rot = x_img / params->world_to_pixel_scale - rot_offsets[0];
  double y_rot = y_img / params->world_to_pixel_scale - rot_offsets[1];

  // Convert (x_rot,y_rot) into the coordinate system of the non-rotated blob image
  *x = mean_x + (x_rot-mean_x)*COS_ORI(ori) - (y_rot-mean_y)*SIN_ORI(ori);
  *y = mean_y + (y_rot-mean_y)*COS_ORI(ori) + (x_rot-mean_x)*SIN_ORI(ori);
}




// Return all points that are both
//   1) Inside the blob
//   2) Equidistant to two contour points (for a given angle and position)
// This function works by computing rotated silhouette images for a discrete
// set of orientations, then extracting the midpoint of blobs for each horizontal
// scanline.  The runtime is linear in the number of pixel locations times the
// number of orientations.  The number of legal points is the number of 
// orientations times the length of the curve through the center of the worm.
SpinePointLocation *get_legal_internal_points(double *contour, int num_pts, 
						 double *offsets, IplImage *blob, 
						 int *num_legal_internal_points, 
						 FitParams *params) {
  int i, j, k, xx, yy, l, add, num;
  int in_blob, start, legal_pt_start;
  unsigned char *ptr;
  int numAlloc = 0;
  double mean_x = 0, mean_y = 0;
  double rot_offsets[2];
  SpinePointLocation *pts = NULL;
  double x, y;
  double *rot_pts = (double*)malloc(2*sizeof(double)*num_pts);
  IplImage *rot_blob;
  int **img_map, **center_map;

  *num_legal_internal_points = 0;

  // compute the mean of the center points
  for(i = 0; i < num_pts; i++) {
    mean_x += contour[2*i];
    mean_y += contour[2*i+1];
  }
  mean_x /= num_pts;
  mean_y /= num_pts;

  img_map = (int**)malloc(blob->height*(sizeof(int*)+blob->width*(sizeof(int)))*2);
  center_map = img_map + blob->height;
  for(i = 0; i < blob->height; i++) {
    img_map[i] = (int*)(((int*)(center_map+blob->height)) + 2*i*blob->width);
    center_map[i] = img_map[i] + blob->width;
    for(j = 0; j < blob->width; j++) 
      center_map[i][j] = img_map[i][j] = -1;
  }

  for(k = 0; k < params->num_orientations; k++) {
    legal_pt_start = *num_legal_internal_points;

    // Compute a rotated version of the contour points such that all spine cross-section
    // lines at this orientation become horizontal, and render the corresponding blob silhouette image
    for(i = 0; i < num_pts; i++) {
      rot_pts[2*i]   = mean_x + (contour[2*i]-mean_x)*COS_N_ORI(k) - (contour[2*i+1]-mean_y)*SIN_N_ORI(k);
      rot_pts[2*i+1] = mean_y + (contour[2*i+1]-mean_y)*COS_N_ORI(k) + (contour[2*i]-mean_x)*SIN_N_ORI(k);
    }
    rot_blob = binary_mask(rot_pts, num_pts, rot_offsets, params);
    ptr = (unsigned char*)rot_blob->imageData;
    
    for(i = 0; i < rot_blob->height; i+=params->pixel_granularity, 
	  ptr += rot_blob->widthStep*params->pixel_granularity) {
      // For each horizontal scanline, solve for the point that is equidistant from 
      // the 2 boundary contour points
      in_blob = 0;
      for(j = 0; j <= rot_blob->width; j++) {
	if(!in_blob && j < rot_blob->width && ptr[j]) {
	  // We are inside a blob
	  in_blob = 1;
	  start = j;
	} else if(in_blob && (j == rot_blob->width || !ptr[j])) {
	  // We are exiting a blob, compute the point that is equidistant from the exit point
	  // and the start point
	  in_blob = 0;
	  rotated_image_coords_to_world_coords(((j-1) + start)/2.0, (double)i, k,
					       mean_x, mean_y, rot_offsets, params, &x, &y);
	  
	  // Add a new legal state
	  if(*num_legal_internal_points+1 > numAlloc) {
	    numAlloc += 128;
	    pts = (SpinePointLocation*)realloc(pts, sizeof(SpinePointLocation)*numAlloc);
	  }
	  pts[*num_legal_internal_points].x = x;
	  pts[*num_legal_internal_points].y = y;
	  pts[*num_legal_internal_points].orientation = k;
	  compute_intersect_points(&pts[*num_legal_internal_points], contour, num_pts, params);
	  
	  xx = my_min(my_max(0,(int)(W2P(x,0)+.5)), blob->width-1);
	  yy = my_min(my_max(0,(int)(W2P(y,1)+.5)), blob->height-1);
	  if((center_map[yy][xx] < 0 || (pts[*num_legal_internal_points].width &&
					 pts[*num_legal_internal_points].width <= pts[center_map[yy][xx]].width)) &&
	     (img_map[yy][xx] < 0 || (pts[*num_legal_internal_points].width &&
				      pts[*num_legal_internal_points].width <= pts[img_map[yy][xx]].width))) {
	    add = (*num_legal_internal_points)++;
	    img_map[yy][xx] = center_map[yy][xx] = add;

	    for(l = start; l < j; l++) {
	      rotated_image_coords_to_world_coords(l, (double)i, k,
						   mean_x, mean_y, rot_offsets, params, &x, &y);
	      xx = my_min(my_max(0,(int)(W2P(x,0)+.5)), blob->width-1);
	      yy = my_min(my_max(0,(int)(W2P(y,1)+.5)), blob->height-1);
	      if(img_map[yy][xx] != add && (img_map[yy][xx] < 0 || (pts[add].width &&
								    pts[add].width <= pts[img_map[yy][xx]].width))) {
		img_map[yy][xx] = add;
	      }
	      if(center_map[yy][xx] != add && (center_map[yy][xx] >= 0 && (pts[add].width &&
									pts[add].width <= pts[center_map[yy][xx]].width))) {
		pts[center_map[yy][xx]].orientation = -1;
	      }
	    }
	  } 
	}
      }
    }
    cvReleaseImage(&rot_blob);

    //char fname[400]; sprintf(fname, "debug/tmp%02d.png", k);
    //cvSaveImage(fname, draw_spine_blob(pts+legal_pt_start, *num_legal_internal_points-legal_pt_start-1, contour, num_pts, 20, params));

  }

  num = *num_legal_internal_points;
  *num_legal_internal_points = 0;
  for(i = 0; i < num; i++) {
    if(pts[i].width > 0 && pts[i].orientation >= 0)
      pts[(*num_legal_internal_points)++] = pts[i];
  }
  
  free(rot_pts);
  free(img_map);

  return pts;
}

// Guess the conversion from world to pixel coordinates, which works only if the number
// of unique x and y values in the contour is equal to the width and height in pixels
double guess_world_to_pixel_scale(double *contour, int num_pts) {
  int i, j, is_duplicate_x, is_duplicate_y;
  int w = 0, h = 0;
  double min_x = INFINITY, max_x = -INFINITY, min_y = INFINITY, max_y = -INFINITY;
  double scale_x, scale_y;

  for(i = 0; i < num_pts; i++) {
    if(contour[2*i] < min_x) min_x = contour[2*i];
    if(contour[2*i] > max_x) max_x = contour[2*i];
    if(contour[2*i+1] < min_y) min_y = contour[2*i+1];
    if(contour[2*i+1] > max_y) max_y = contour[2*i+1];

    is_duplicate_x = is_duplicate_y = 0;
    for(j = 0; j < i; j++) {
      if(contour[2*i] == contour[2*j]) 
	is_duplicate_x = 1;
      if(contour[2*i+1] == contour[2*j+1]) 
	is_duplicate_y = 1;
    }
    if(!is_duplicate_x) 
      w++;
    if(!is_duplicate_y) 
      h++;
  }
  scale_x = (w-1) / (max_x-min_x);
  scale_y = (h-1) / (max_y-min_y);
  assert(scale_x/scale_y < 1.05 && scale_x/scale_y > .95);
  
  return (scale_x+scale_y)/2;
}






DPCache *new_dp_cache(Blob *b, FitParams *params) {
  DPCache *d = (DPCache*)malloc(sizeof(DPCache));
  int n;

  d->blob = b;
  d->params = params;
  
  if(!fix_contour(b->contour, b->num_pts, d->offset, params, &d->contour, 
		  &d->num_pts, &d->img)) {
    d->img = binary_mask(b->contour, b->num_pts, d->offset, params);
    d->contour = b->contour;
    d->num_pts = b->num_pts;
  } else {
    b->fixed_contour = d->contour;
    b->fixed_num_pts = d->num_pts;
  }
  d->legal_points = get_legal_internal_points(d->contour, d->num_pts, d->offset, d->img, 
					      &d->num_legal_points, params);
  d->table_costs = (double**)malloc((params->num_spine_lines+1)*sizeof(double*));
  d->table_states = (int**)malloc((params->num_spine_lines+1)*sizeof(int*));
  
  // Allocate cache tables
  for(n = 0; n < params->num_spine_lines+1; n++) {  
    d->table_costs[n] = (double*)malloc(sizeof(double)*d->num_legal_points);
    d->table_states[n] = (int*)malloc(sizeof(int*)*d->num_legal_points);
  }

  return d;
}

void free_dp_cache(DPCache *d) {
  int n;

  if(d) {
    if(d->contour != d->blob->contour && d->contour != d->blob->fixed_contour)
      free(d->contour);
    cvReleaseImage(&d->img);
    free(d->legal_points);
    for(n = 0; n < d->params->num_spine_lines+1; n++) {
      free(d->table_costs[n]);
      free(d->table_states[n]);
    }
    if(d->transitions) 
      free_pairwise_spine_transitions(d->transitions, d->num_legal_points);

    free(d->table_costs);
    free(d->table_states);
    free(d);
  }
}

// Draw a visualization of a blob and its skeleton
IplImage *draw_blob(IplImage *img, Blob *b, double *offsets, FitParams *params, double zoom, BehaviorGroups *behaviors) {
  double tmp_offsets[2];
  int num = b->num_model_pts;
  CvPoint *pts = (CvPoint*)malloc(sizeof(CvPoint)*(num > b->num_pts ? num : b->num_pts));
  double min_x = INFINITY, max_x = -INFINITY, min_y = INFINITY, max_y = -INFINITY;
  int i, w, h;
  SpinePointLocation *spine = b->num_model_pts > 0 ? b->model_pts : NULL;
  double *contour = b->fixed_contour ? b->fixed_contour : b->contour;
  double *contour2 = b->contour;
  int num_pts = b->num_pts;
#if SHOW_FIXED_CONTOUR
  num_pts = b->fixed_contour ? b->fixed_num_pts : b->num_pts;
  contour2 = contour;
#endif

  if(!offsets) {
    offsets = tmp_offsets;
    blob_get_bounding_box(&b, 1, &min_x, &min_y, &max_x, &max_y);
    offsets[0] = -min_x;
    offsets[1] = -min_y;
    w = ceil((max_x-min_x)*params->world_to_pixel_scale*zoom)+1;
    h = ceil((max_y-min_y)*params->world_to_pixel_scale*zoom)+1;
  }

  if(!img) {
    img = cvCreateImage(cvSize(w,h), IPL_DEPTH_8U, 3);
    cvZero(img);
  }

  // Draw the blob contour
  for(i = 0; i < num_pts; i++) {
    pts[i].x = (int)(.5+W2P(contour2[2*i],0)*zoom);
    pts[i].y = (int)(.5+W2P(contour2[2*i+1],1)*zoom);
  }
  cvFillPoly(img, &pts, &num_pts, 1, CV_RGB(128,128,128));
  cvPolyLine(img, &pts, &num_pts, 1, 0, CV_RGB(255,255,255), 2);

  // Draw the spine cross-sections
  for(i = 0; i < num; i++) {
    if(spine[i].p1 >= 0 && spine[i].p2 >= 0) {
      CvPoint p1 = cvPoint(W2P(contour[2*spine[i].p1],0)*zoom, W2P(contour[2*spine[i].p1+1],1)*zoom),
	p2 = cvPoint(W2P(contour[2*spine[i].p2],0)*zoom, W2P(contour[2*spine[i].p2+1],1)*zoom);
      cvLine(img, p1, p2, CV_RGB(0,0,255));
      if(spine[i].p1 >= 0)
	cvCircle(img, p1, 3, CV_RGB(0,255,0), 1);
      if(spine[i].p2 >= 0)
	cvCircle(img, p2, 3, CV_RGB(255,255,0), 1);
    }
  }
  

  // Draw the actual spine
  for(i = 0; i < num; i++) {
    pts[i].x = (int)(.5+W2P(spine[i].x,0)*zoom);
    pts[i].y = (int)(.5+W2P(spine[i].y,1)*zoom);
  }
  cvPolyLine(img, &pts, &num, 1, 0, CV_RGB(255,0,0), 2);
  
  if(spine) {
    if(spine[0].p4 >= 0) {
      cvLine(img, cvPoint(pts[0].x, pts[0].y), 
	     cvPoint(W2P(contour[2*spine[0].p4],0)*zoom, 
		     W2P(contour[2*spine[0].p4+1],1)*zoom), CV_RGB(255,0,255), 2);
      cvCircle(img, cvPoint(W2P(contour[2*spine[0].p4],0)*zoom, 
			    W2P(contour[2*spine[0].p4+1],1)*zoom),
	       5, CV_RGB(255,0,255), -1);
    }
    
    if(spine[num-1].p3 >= 0) {
      cvLine(img, cvPoint(pts[num-1].x, pts[num-1].y), 
	     cvPoint(W2P(contour[2*spine[num-1].p3],0)*zoom, 
		     W2P(contour[2*spine[num-1].p3+1],1)*zoom), CV_RGB(255,255,0), 2);
      cvCircle(img, cvPoint(W2P(contour[2*spine[num-1].p3],0)*zoom, 
			    W2P(contour[2*spine[num-1].p3+1],1)*zoom),
	       5, CV_RGB(255,255,0), -1);
    }
  }

  for(i = 0; i < num; i++) {
    cvCircle(img, cvPoint(pts[i].x, pts[i].y), 5, i == b->on_pt-1 ? CV_RGB(0,255,0) : CV_RGB(0,0,255), -1);
  }
  /*if(g_debug && b->fixed_contour) {
    for(i = 0; i < b->fixed_num_pts; i++) {
      if(i == 27 || (i >= 10 && i <= 13)) {
	cvCircle(img, cvPoint(W2P(b->fixed_contour[2*i],0)*zoom+.5, W2P(b->fixed_contour[2*i+1],1)*zoom+.5), 5, 
		 i==12 || i==13 ? CV_RGB(0,255,255) : CV_RGB(255,0,0), 1);
      }
    }
    cvLine(img, cvPoint(W2P(contour[2*spine[1].p1],0)*zoom, W2P(contour[2*spine[1].p1+1],1)*zoom), 
	   cvPoint(W2P(contour[2*spine[2].p1],0)*zoom, W2P(contour[2*spine[2].p1+1],1)*zoom), cvScalar(0,255,0), 3);
    cvLine(img, cvPoint(W2P(contour[2*spine[1].p2],0)*zoom, W2P(contour[2*spine[1].p2+1],1)*zoom), 
	   cvPoint(W2P(contour[2*spine[2].p2],0)*zoom, W2P(contour[2*spine[2].p2+1],1)*zoom), cvScalar(0,255,0), 3);
	   }*/

  /*if(behaviors) {
    char str[1000], tmp[400];
    strcpy(str, "");
    for(i = 0; i < behaviors->num; i++) {
      if(b->behaviors[i] >= 0) {
	if(i) strcat(str, ", ");
	sprintf(tmp, "%s: %s", behaviors->behaviors[i].name, behaviors->behaviors[i].values[b->behaviors[i]].name);
	strcat(str, tmp);
      }
    }
    CvFont font;
    cvInitFont(&font, CV_FONT_VECTOR0, 0.5f, 0.4f, 0, 2);
    cvPutText(img, str, cvPoint(0, img->height-8), &font, CV_RGB(0,255,0));
    }*/

  free(pts);
	     
  return img;
}





FitParams default_parameters() {
  FitParams params;
  int i;

  double length_prior = .65;

  // Default parameters
  params.num_spine_lines = 7;
  params.num_orientations = 30;
  params.pixel_granularity = 2;
  params.world_to_pixel_scale = 1;
  params.expected_length = length_prior / params.num_spine_lines;
  params.fit_method = FIT_GREEDY_TWO_PASS;//FIT_FRAME_BY_FRAME;//FIT_GREEDY_ONE_PASS;//FIT_GREEDY_TWO_PASS;
  params.speed_estimate_window = .5;

  params.is_multiclass = 1;
  params.classifier_method = CLASSIFIER_SVM_STRUCT_BEHAVIOR;//CLASSIFIER_SVM_LIGHT;//CLASSIFIER_NORMAL_BAYES;//CLASSIFIER_DECISION_TREE;;//CLASSIFIER_SVM;
  params.feature_diff_frames = 10;


  params.max_position_change = 12*params.expected_length;//3*params.expected_length;
  params.max_endpoint_length = 6*params.expected_length;//1.5*params.expected_length;
  params.max_angle_change    = M_PI/3;
  params.max_endpoint_ratio   = INFINITY;//.65;
  params.max_width_change   = INFINITY;//1.5;
  params.max_width = 4*length_prior;//length_prior;
  params.max_outlier_dist = 1.5;
  params.blob_min_pixel_width = 1;

  for(i = 0; i < F_NUM_FEATURES; i++) {
    params.mu[i] = 0;
    params.var[i] = 1;
  }

  params.mu[F_WIDTH] = length_prior / 4.33;
  params.var[F_WIDTH] = SQR(3*params.mu[F_WIDTH]/3);//SQR(params.mu[F_WIDTH]/3);
  params.mu[F_HEAD_LENGTH] = params.mu[F_TAIL_LENGTH] = params.expected_length/2;
  params.var[F_HEAD_LENGTH] = params.var[F_TAIL_LENGTH] = SQR(4*params.expected_length/4); //SQR(params.expected_length/4);
  params.var[F_HEAD_SYMMETRY] = params.var[F_TAIL_SYMMETRY] = SQR(.02);
  params.var[F_HEAD_GEODESIC_SYMMETRY] = params.var[F_TAIL_GEODESIC_SYMMETRY] = SQR(.02);
  params.var[F_HEAD_PERCENT_OUTLIERS] = params.var[F_TAIL_PERCENT_OUTLIERS] = SQR(.15);
  
  params.var[F_WIDTH_VELOCITY] = SQR(3*params.mu[F_WIDTH]/5/.08);//SQR(params.mu[F_WIDTH]/5/.08);
  params.var[F_FORWARD_VELOCITY] = SQR(3*length_prior/25/.08); //SQR(length_prior/25/.08);
  params.var[F_SIDEWAYS_VELOCITY] = SQR(6*length_prior/50/.08); //SQR(length_prior/50/.08);
  params.var[F_ANGULAR_VELOCITY] = SQR(.8/.08);
    
  params.var[F_SEGMENT_ANGULAR_VELOCITY] = SQR(.8/.08);
  params.var[F_LENGTH_VELOCITY] = SQR(length_prior/20/.08);

  params.var[F_PERCENT_OUTLIERS] = SQR(.01);
  params.var[F_WIDTH_CHANGE] = SQR(.15);
  params.var[F_ANGLE_CHANGE] = SQR(.8);
  params.mu[F_LENGTH] = params.expected_length;
  params.var[F_LENGTH] = SQR(3*params.expected_length/3);//SQR(params.expected_length/3);
  params.var[F_LEFT_PARALLEL] = params.var[F_RIGHT_PARALLEL] = 1;
  params.var[F_LEFT_ORTHOGONAL] = params.var[F_RIGHT_ORTHOGONAL] = SQR(.3);
  params.var[F_SYMMETRY] = SQR(.8); //SQR(.03);
  params.var[F_HEAD_GEOMETRY] = params.var[F_TAIL_GEOMETRY] = 
    params.var[F_GEOMETRY1] = params.var[F_GEOMETRY2] = SQR(sin(2*M_PI/(15+params.num_orientations)));




  // Cache values of orientations and their sines and cosines
  params.orientations = (double*)malloc(sizeof(double)*params.num_orientations*3);
  params.cos_orientations = params.orientations + params.num_orientations;
  params.sin_orientations = params.cos_orientations + params.num_orientations;
  for(i = 0; i < params.num_orientations; i++) {
    params.orientations[i] = i * 2*M_PI / params.num_orientations;
    params.cos_orientations[i] = cos(params.orientations[i]);
    params.sin_orientations[i] = sin(params.orientations[i]);
  }

  return params;
}


#ifdef STANDALONE
int main(int argc, char **argv) {
  FitParams params;
  BlobSequence *s;
  char blob_file[400];
  int i;

  if(argc < 2) { 
    fprintf(stderr, "USAGE: ./fit_model file.outline <num_spine_lines>\n");
    return -1;
  }
  strcpy(blob_file, argv[1]);

  s = import_blob_sequence(blob_file);
  if(!s || s->num_frames <= 0) {
    fprintf(stderr, "Couldn't read blob file %s. Aborting...\n", blob_file);
    return -1;
  }

  params = default_parameters();
  params.world_to_pixel_scale = guess_world_to_pixel_scale(s->frames[0].contour, s->frames[0].num_pts);

  fit_model(s, &params);

  free_blob_sequence(s);
  free(params.orientations);
}
#endif
