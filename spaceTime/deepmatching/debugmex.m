i1 = imread('im1.png');
i2 = imread('im2.png');
%setenv('BLAS_VERSION','/usr/lib/atlas-base/libsatlas.so');
k = deepmex(single(i1),single(i2),320,320,122);
%setenv('BLAS_VERSION','/usr/local/MATLAB/R2015a/bin/glnxa64/mkl.so');

exit
