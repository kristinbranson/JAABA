INCDIR = -I../include -I../../online_structured_svm/json -I../../online_structured_svm/include -I../../blob/include

INCDIR_EXTRA = -I/opt/matlab/extern/include -I/usr/local/MATLAB/R2012a/extern/include
LIBS_EXTRA = -L/opt/matlab/extern/lib -L/opt/matlab/bin/glnxa64 -lmat -lmx -Wl,-rpath,/opt/matlab/bin/glnxa64 -L/usr/local/MATLAB/R2012a/bin/glnxa64 -Wl,-rpath,/usr/local/MATLAB/R2012a/bin/glnxa64

#LIBS_CV = -lcv -lcxcore -lhighgui # opencv-1.0
LIBS_CV = -lopencv_highgui -lopencv_imgproc -lopencv_core # opencv-2.0


CC=g++
CFLAGS= $(INCDIR) $(INCDIR_EXTRA) -fopenmp -DUSE_OPENMP 
CFLAGS_DEBUG= -g -Wall
CFLAGS_RELEASE= -O2 -fomit-frame-pointer -ffast-math
CFLAGS_SHARED= -fPIC
LIBS= -L../lib -L../../online_structured_svm/lib -L../../online_structured_svm/json/lib  -L../../blob/lib -lm $(LIBS_EXTRA)
LIBS_RELEASE= -lonline_structured_svm -ljson -lbehavior -lblob $(LIBS_CV)
LIBS_DEBUG= -ldonline_structured_svm -ldjson -ldbehavior -ldblob $(LIBS_CV)
LD=g++
LDFLAGS= -shared
AR=ar
ARFLAGS=rcs



CPP_FILES := $(wildcard *.cpp)
HEADER_FILES = *.h ../include/*.h ../../online_structured_svm/include/*.h ../../online_structured_svm/json/json/*.h ../../blob/include/*.h
OUT_FILES := $(addprefix bin/$(CONFIGURATION)/,$(notdir $(CPP_FILES:.cpp=.out)))


all: debug_static 


debug_shared: OUT_FILES=$(addprefix bin/debug_shared/,$(notdir $(CPP_FILES:.cpp=.out)))
debug_shared: LIBS += $(LIBS_DEBUG)
debug_shared: $(addprefix bin/debug_shared/,$(notdir $(CPP_FILES:.cpp=.out)))

release_shared: OUT_FILES=$(addprefix bin/release_shared/,$(notdir $(CPP_FILES:.cpp=.out)))
release_shared: LIBS += $(LIBS_RELEASE)
release_shared: $(addprefix bin/release_shared/,$(notdir $(CPP_FILES:.cpp=.out)))

debug_static: OUT_FILES=$(addprefix bin/debug_static/,$(notdir $(CPP_FILES:.cpp=.out)))
debug_static: LIBS += $(LIBS_DEBUG)
debug_static: $(addprefix bin/debug_static/,$(notdir $(CPP_FILES:.cpp=.out)))

release_static: OUT_FILES=$(addprefix bin/release_static/,$(notdir $(CPP_FILES:.cpp=.out)))
release_static: LIBS += $(LIBS_RELEASE)
release_static: $(addprefix bin/release_static/,$(notdir $(CPP_FILES:.cpp=.out)))


clean: 
	rm -f bin/*/*.out


bin/debug_static/%.out: %.cpp $(HEADER_FILES) ../lib/libdbehavior.a ../../blob/lib/libdblob.a ../../online_structured_svm/lib/libdonline_structured_svm.a ../../online_structured_svm/json/lib/libdjson.a 
	$(CC) $(CFLAGS) $(CFLAGS_DEBUG) $< -o $@ $(LIBS)

bin/debug_shared/%.out: %.cpp $(HEADER_FILES) 
	$(CC) $(CFLAGS) $(CFLAGS_DEBUG) $(CFLAGS_RELEASE) $< -o $@ $(LIBS)

bin/release_static/%.out: %.cpp $(HEADER_FILES) ../lib/libbehavior.a ../../blob/lib/libblob.a ../../online_structured_svm/lib/libonline_structured_svm.a ../../online_structured_svm/json/lib/libjson.a 
	$(CC) $(CFLAGS) $< -o $@ $(LIBS)

bin/release_shared/%.out: %.cpp $(HEADER_FILES) 
	$(CC) $(CFLAGS) $(CFLAGS_RELEASE) $< -o $@ $(LIBS)


