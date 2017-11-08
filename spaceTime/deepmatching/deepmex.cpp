#include "mex.h"
#include "image.h"
#include "deep_matching.h"
#include "io.h"
#include <stdio.h>

#define I1 0
#define I2 1
#define RX 2
#define RY 3
#define NR 4

#define DEBUG 0

image_t* rescale_image1( image_t* im, int width, int height ) 
{
  image_t* res = image_new(width,height);
  image_resize_bilinear_newsize(res, im, width, height);
  image_delete(im);
  return res;
}

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[]){

  // main inputs
  float *img1;
  float *img2;
  int rescalex;
  int rescaley;
  int ngh_rad;
  image_t *im1=NULL, *im2=NULL;

  int imszx,imszy;

  if(DEBUG)
    mexPrintf("Starting mex\n");

  // check arguments
  if( (nrhs != 5) || (nlhs>1) ){
    mexErrMsgTxt("Usage: matches = deepmex(im1,im2,rescalex,rescaley,ngh_rad)");
    return;
  }


  // check format of arguments
  if( !mxIsSingle(prhs[I1]) ){
    mexErrMsgTxt("img1 must be of class single.");
    return;
  }
  if( !mxIsSingle(prhs[I2]) ){
    mexErrMsgTxt("img2 must be of class single.");
    return;
  }
//  if( !mxIsUint64(prhs[RX]) ){
//    mexErrMsgTxt("RX must be of class uint64.");
//    return;
//  }
//  if( !mxIsUint64(prhs[RY]) ){
//    mexErrMsgTxt("idx must be of class uint64.");
//    return;
//  }
//  if( !mxIsUint64(prhs[NR]) ){
//    mexErrMsgTxt("idx must be of class uint64.");
//    return;
//  }

  // get input sizes
  imszy = (int)mxGetM(prhs[I1]);
  imszx = (int)mxGetN(prhs[I1]);

  // make sure input sizes match
  if(imszy != (int)mxGetM(prhs[I2])){
    mexErrMsgTxt("Number of rows in images dont match");
    return;
  }
  if(imszx != (int)mxGetN(prhs[I2])){
    mexErrMsgTxt("Number of columns in images dont match");
    return;
  }
  if(!mxIsNumeric(prhs[RX]) || !mxIsScalar(prhs[RX])){
    mexErrMsgTxt("rescalex should be numeric and scalar");
    return;
  }
  if(!mxIsNumeric(prhs[RY]) || !mxIsScalar(prhs[RY])){
    mexErrMsgTxt("rescaley should be numeric and scalar");
    return;
  }
  if(!mxIsNumeric(prhs[NR]) || !mxIsScalar(prhs[NR])){
    mexErrMsgTxt("ngh_rad should be numeric and scalar");
    return;
  }

  // get inputs
  img1 = (float *)mxGetData(prhs[I1]);
  img2 = (float *)mxGetData(prhs[I2]);

  rescalex = (int)mxGetScalar(prhs[RX]);
  rescaley = (int)mxGetScalar(prhs[RY]);
  ngh_rad = (int)mxGetScalar(prhs[NR]);
  
  if(DEBUG)
    mexPrintf("Input reading done\n");
//Read the images as color because color_image_t makes more sense than image_t
//size allocation of new image of image_t seems off because it uses sizeof(float) twice. So better copy the data manually
  color_image_t* it1 = color_image_new(imszx,imszy);
  color_image_t* it2 = color_image_new(imszx,imszy);
  for(int i=0; i<imszx*imszy; i++) {
    it1->c1[i] = it1->c2[i] = it1->c3[i] = img1[i];
    it2->c1[i] = it2->c2[i] = it2->c3[i] = img2[i];
    if(DEBUG)
     {//mexPrintf("%i %g\n",i,img1[i]);
     }
  }
  

  im1 = image_gray_from_color(it1);
  im2 = image_gray_from_color(it2);
  color_image_delete(it1);
  color_image_delete(it2);

  if(DEBUG)
    mexPrintf("done converting color imgs to gray\n");
  dm_params_t params;
  set_default_dm_params(&params);
 
  //png settings
  params.desc_params.presmooth_sigma = 0; // no image smoothing since the image is uncompressed
  params.desc_params.hog_sigmoid = 0.2;
  params.desc_params.mid_smoothing = 1.5;
  params.desc_params.post_smoothing = 1;
  params.desc_params.ninth_dim = 0.1; // low ninth_dim since image PSNR is high
  params.ngh_rad = ngh_rad;
  float fx=1, fy=1;
  fx *= im1->width / float(rescalex);
  fy *= im1->height / float(rescaley);
  im1 = rescale_image1(im1, rescalex, rescaley);
  im2 = rescale_image1(im2, rescalex, rescaley);

  if(DEBUG)
    mexPrintf("done resizing the image\n");

  float_image* corres = deep_matching( im1, im2, &params, NULL );  // standard call

  if(DEBUG)
   mexPrintf("done matching the images\n");

  plhs[0] = mxCreateNumericMatrix(corres->tx,corres->ty,mxSINGLE_CLASS,mxREAL);

  float * temp = (float *) mxGetData(plhs[0]);
  for(int i=0; i<corres->ty*corres->tx;i++)
   temp[i] = corres->pixels[i];

  image_delete(im1);
  image_delete(im2);
  free(corres->pixels);
  free(corres);
}

