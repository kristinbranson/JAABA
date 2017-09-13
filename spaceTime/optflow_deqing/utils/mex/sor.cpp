/*
  sor.cpp - Linear equation solver based on successive over-relaxation method

  This is a Matlab v7 MEX plugin.  Compile using "mex -largeArrayDims -O sor.cpp".
  
  Input arguments:
   - A:           Coefficient matrix (transposed)
   - b:           Right-hand side vector
   - w:           Over-relaxation parameter
   - niters:      Maximum number of iterations
   - tol:         Residual convergence tolerance
   - x0:          Starting value for x (optional)
   
  Outputs:
   - x:           Approximate solution for A'x = b
   - flag:        Convergence indicator:  true if converged (optional)
   - res:         Residual (optional)
   - n:           Number of iterations (optional)

  
  Author:  Stefan Roth, Department of Computer Science, TU Darmstadt
  Contact: sroth@cs.tu-darmstadt.de
  $Date$
  $Revision$
  
 * Copyright 2004-2007, Brown University, Providence, RI. USA
 * Copyright 2007-2010 TU Darmstadt, Darmstadt, Germany.
 *
 * All Rights Reserved
 *
 * All commercial use of this software, whether direct or indirect, is
 * strictly prohibited including, without limitation, incorporation into in
 * a commercial product, use in a commercial service, or production of other
 * artifacts for commercial purposes.
 *
 * Permission to use, copy, modify, and distribute this software and its
 * documentation for research purposes is hereby granted without fee,
 * provided that the above copyright notice appears in all copies and that
 * both that copyright notice and this permission notice appear in
 * supporting documentation, and that the name of the author and Brown
 * University not be used in advertising or publicity pertaining to
 * distribution of the software without specific, written prior permission.
 *
 * For commercial uses contact the Technology Venture Office of Brown University
 *
 * THE AUTHOR AND BROWN UNIVERSITY DISCLAIM ALL WARRANTIES WITH REGARD TO
 * THIS SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND
 * FITNESS FOR ANY PARTICULAR PURPOSE.  IN NO EVENT SHALL THE AUTHOR OR
 * BROWN UNIVERSITY BE LIABLE FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL
 * DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR
 * PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS
 * ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF
 * THIS SOFTWARE.
 *
 
 */

#include <iostream>
#include <sstream>
#include "mex.h"

#include "residual.h"


extern "C"
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{

  /* Check for proper number of input and output arguments */
  if (nrhs < 5 || nrhs > 6)
    mexErrMsgTxt("5 or 6 input arguments required.");
  
  /* Check data type of input argument */
  for (int i = 0; i < nrhs; i++)
    {
      if (!mxIsDouble(prhs[i]))
	{
	  std::ostringstream s;
	  s << "Input " << i + 1 << " must be of type double.";
	  mexErrMsgTxt(s.str().c_str());
	}
      if (mxIsComplex(prhs[i]))
	{
	  std::ostringstream s;
	  s << "Input " << i + 1 << " may not be complex.";
	  mexErrMsgTxt(s.str().c_str());
	}
      if (mxGetNumberOfDimensions(prhs[i]) > 2)
	{
	  std::ostringstream s;
	  s << "Input " << i + 1 << " may not have more than 2 dimensions.";
	  mexErrMsgTxt(s.str().c_str());
	}
    }

  for (int i = 2; i < 5; i++)
    if (mxGetM(prhs[i]) != 1 || mxGetN(prhs[i]) != 1)
      {
	std::ostringstream s;
	s << "Input " << i + 1 << " must be scalar.";
	mexErrMsgTxt(s.str().c_str());
      }

  
  mwSize ai = mxGetM(prhs[0]);
  mwSize aj = mxGetN(prhs[0]);
  mwSize bi = mxGetM(prhs[1]);
  mwSize bj = mxGetN(prhs[1]);

  if (ai <= 0 || aj <= 0)
    mexErrMsgTxt("First argument may not be empty.");
  if (bi <= 0 || bj <= 0)
    mexErrMsgTxt("Second argument may not be empty.");
  if (ai != aj)
    mexErrMsgTxt("Coefficient matrix must be square.");    
  if (ai != bi)
    mexErrMsgTxt("Argument sizes incompatible.");
  if (bj > 1)
    mexErrMsgTxt("Right-hand side vector may only have one column.");

  if (nrhs > 5)
    {
      mwSize xi = mxGetM(prhs[5]);
      mwSize xj = mxGetN(prhs[5]);

      if (ai != xi)
	mexErrMsgTxt("Argument sizes incompatible.");
      if (xj > 1)
	mexErrMsgTxt("Initialization vector may only have one column.");
    }
  
  if (!mxIsSparse(prhs[0]))
    mexErrMsgTxt("Coefficient matrix must be sparse.");    
  if (mxIsSparse(prhs[1]))
    mexErrMsgTxt("Right hand size may not be sparse.");    
  

  double*  APr = mxGetPr(prhs[0]);
  mwIndex* AJc = mxGetJc(prhs[0]);
  mwIndex* AIr = mxGetIr(prhs[0]);
  double*  bPr = mxGetPr(prhs[1]);

  double w      = mxGetScalar(prhs[2]);
  int    niters = (int) mxGetScalar(prhs[3]);
  double tol    = mxGetScalar(prhs[4]);


  // -----------------------------------------------------------------
  // Construct output
  // -----------------------------------------------------------------
  
  plhs[0] = mxCreateDoubleMatrix(bi, 1, mxREAL);
  double* xPr  = mxGetPr(plhs[0]);

  if (nrhs > 5)
    memcpy(xPr, mxGetPr(prhs[5]), bi * sizeof(double));
  
  // -----------------------------------------------------------------
  // Perform SOR 
  // -----------------------------------------------------------------

  for (int n = 0; n < niters; n++)
    {
      if ((n % 5) == 0)
	{
	  // Terminate if residual treshold reached
	  double res = residual(prhs[0], plhs[0], prhs[1]);
	  if (res < tol)
	    {
	      if (nlhs > 1)
		plhs[1] = mxCreateLogicalScalar(true);
	      if (nlhs > 2)
		plhs[2] = mxCreateDoubleScalar(res);
	      if (nlhs > 3)
		plhs[3] = mxCreateDoubleScalar(n);
	      return;
	    }
	}
      
      double* A = APr;
      
      // Iterate over all the columns of A
      for (mwSize col = 0; col < ai; col++)
	{
	  mwIndex row_start = AJc[col];
	  mwIndex nrows = AJc[col+1] - row_start;
	  mwIndex* row = AIr + row_start;
	  
	  double sigma = 0.0;	    
	  double diag  = 0.0;

	  // Iterate over all non-zero row indices
	  for (mwSize row_idx = nrows; row_idx > 0; row_idx--)
	    {
	      double  tmp = *A;
	      mwIndex r = *row;

	      if (r == col)
		diag = tmp;
	      else
		sigma += tmp * xPr[r];

	      A++;	
	      row++;
	    }
	  
	  sigma = (bPr[col] - sigma) / diag;
	  xPr[col] += w * (sigma - xPr[col]);
	}
    }
  
  if (nlhs > 1)
    plhs[1] = mxCreateLogicalScalar(false);
  if (nlhs > 2)
    plhs[2] = mxCreateDoubleScalar(residual(prhs[0], plhs[0], prhs[1]));
  if (nlhs > 3)
    plhs[3] = mxCreateDoubleScalar(niters);
}
