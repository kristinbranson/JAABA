#mex -v -c conv.cpp *.o  -ljpeg -lpng CXXFLAGS=' -fopenmp -fPIC' LDFLAGS='/usr/lib/atlas-base/libsatlas.so' -lgomp
#mex -v -c deep_matching.cpp *.o  -ljpeg -lpng CXXFLAGS=' -fopenmp -fPIC' LDFLAGS='/usr/lib/atlas-base/libsatlas.so' -lgomp
#bash mexCmd
#mex -v -c conv.cpp *.o  -ljpeg -lpng -lmwblas -lmwlapack CXXFLAGS=' -fopenmp -fPIC' -lgomp 
#mex -v -c deep_matching.cpp *.o  -ljpeg -lpng -lmwblas -lmwlapack CXXFLAGS=' -fopenmp -fPIC' -lgomp 
#mex -v deepmex.cpp *.o  -ljpeg -lpng -lmwblas -lmwlapack CXXFLAGS=' -fopenmp -fPIC' -lgomp
#without openmp
export FLAGS=" -lmwblas -lmwlapack "

rm -rf *.o deepmex
mex -v -c hog.cpp $FLAGS CXXFLAGS=' -fPIC'  LDFLAGS=' -fPIC'
mex -v -c maxfilter.cpp $FLAGS CXXFLAGS=' -fPIC'  LDFLAGS=' -fPIC'
mex -v -c io.cpp $FLAGS CXXFLAGS=' -fPIC'  LDFLAGS=' -fPIC'
mex -v -c image.cpp $FLAGS CXXFLAGS=' -fPIC'  LDFLAGS=' -fPIC'
mex -v -c pixel_desc.cpp $FLAGS CXXFLAGS=' -fPIC'  LDFLAGS=' -fPIC'
mex -v -c conv.cpp *.o  $FLAGS CXXFLAGS=' -fPIC'  LDFLAGS=' -fPIC'
mex -v -c conv.cpp *.o  $FLAGS CXXFLAGS=' -fPIC' LDFLAGS=' -fPIC'
mex -v -c deep_matching.cpp *.o $FLAGS CXXFLAGS=' -fPIC'  LDFLAGS=' -fPIC' 
mex -v deepmex.cpp *.o  $FLAGS CXXFLAGS=' -fPIC' LDFLAGS=' -fPIC' -output deepmex_cluster
/usr/local/matlab-2015a/bin/matlab -nodesktop -nosplash -r debugmex_cluster 
