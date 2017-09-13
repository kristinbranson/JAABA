Implementation of the Deep Matching algorithm, published at ICCV 2013 in
"DeepFlow: Large displacement optical flow with deep matching" by Philippe 
Weinzaepfel, Jerome Revaud, Zaid Harchaoui and Cordelia Schmid.
Code and idea by Jerome Revaud, INRIA. The code is only for scientific 
or personnal use. Please contact me/INRIA for commercial use.
Email: jerome.revaud@inria.fr

Copyright (C) 2015 Jerome Revaud

Version 1.2

License:

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>


Installation:
  
  make clean all
  
  This program has been built on a fedora18 x64 machine and tested on Mac OS X. 
  *No assistance* will be given to compile the code on other OS. However, if 
  you are able to sucessfully adapt the code for other platforms (Windows, 
  Matlab wrapper etc.), please notify me so that I can release these versions on 
  the webpage:
  
    http://lear.inrialpes.fr/src/deepmatching/



Example usages and explanations:
  
  To get detailed information on parameters:
    ./deepmatching -h
    ./deepmatching --help
  
  
  * Build verification:
      ./deepmatching liberty1.png liberty2.png -downscale 2 -v
      
    should produce the following output:
      layer 0, patch_size = 16x16
      remaining 16 big cells (actually, 16 are unique)
      layer 1, patch_size = 32x32
      remaining 25 big cells (actually, 25 are unique)
      layer 2, patch_size = 64x64
      remaining 25 big cells (actually, 25 are unique)
      found 625 local matches
      gathering correspondences 96%...
      8 8 0 12 2.56704 9
      8 40 4 48 2.52681 13
      8 24 8 32 2.39697 13
      40 40 40 32 2.58871 0
      40 56 40 48 2.51677 0
      40 24 40 12 2.59455 0
      56 40 56 28 2.56667 0
      56 24 56 12 2.62056 0
      24 40 24 32 2.54901 2
      24 56 24 48 2.49202 4
      24 24 24 16 2.4571 2
  
  * To visualize the output correspondences:
    Use the "viz.py" python script provided.
      ./deepmatching climb1.png climb2.png -nt 0 | python viz.py climb1.png climb2.png
  
  * To restrict matching to local neighborhood:
    The "-ngh_rad <D>" option restricts the matching to a radius of <D> pixels.
    It uses less memory and is faster. For instance, This should produce about 
    the same output as before but costs about 2 times less memory and computations:
    
      ./deepmatching climb1.png climb2.png -nt 0 -ngh_rad 192 | python viz.py climb1.png climb2.png
  
 * To rescore matches prior to calling deepflow / epicflow:
    simply pipe the output correspondences in 'rescore.py'
      ./deepmatching img1 img2 [args] | python rescore.py img1 img2
  
  
 * Scale and invariant version: (see the --help)
      ./deepmatching dino1.jpg dino2.jpg -nt 0 -downscale 1 -max_scale 2 -rot_range -45 +45 -v | python viz.py dino1.jpg dino2.jpg
    
    param -max_scale: maximum scale factor (here x2, default = x5)
    param -rot_range: rotation range in degrees (default = from 0 to 360)


For details about the options, please refer to the help, the papers or the code.


Important tip:
  If the program stops with "segmentation fault", then it means that your machine 
  does not have enough memory. In this case, you should consider increasing the 
  "-downscale" parameter.


Version history:

  version 1.0.2:
    Many thanks to Bowen Zhang from Tongji University for reporting an issue with the makefile

  version 1.1:
  - New mode added for "fully scale & rotation invariant DeepMatching".
  - Improved visualisation (viz.py) 
  - Removed useless/suboptimal options (-iccv_settings)
  - Fixed a bug related to memory allocation for large images

  version 1.2:
  - Added a new option "-ngh_rad" to restrict the matching to a local neighborhood, which allows
    much reduced memory usage and computations.
  - static-compiled version is now fully multhi-threaded with BLAS
  - few minor bugfix, code cleaning and updates.

































