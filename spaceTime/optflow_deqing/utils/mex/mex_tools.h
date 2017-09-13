#ifndef _MEX_TOOLS_H_
#define _MEX_TOOLS_H_

/*
  mex_tools.h - MEX utility functions

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

#include "mex.h"

#ifdef __cplusplus
extern "C"
{
#endif
  
/*
  Create a double matrix of size m x n with given complexity flag, but
  in contrast to MEX-builtin don't initialize values.
 */    
mxArray*
mxCreateDoubleMatrixNoinit(mwSize m, mwSize n, mxComplexity ComplexFlag);

  
/*
  Create a double array with ndim dimensions and sizes dims with given
  complexity flag, but in contrast to MEX-builtin don't initialize
  values.
 */    
mxArray*
mxCreateDoubleArrayNoinit(mwSize ndim, const mwSize* dims, mxComplexity ComplexFlag);

#ifdef __cplusplus
}
#endif

#endif /* _MEX_TOOLS_H_ */
