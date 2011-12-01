INCDIR = -I../include -I../json

CC=g++
CFLAGS= $(INCDIR) -fopenmp -DUSE_OPENMP
CFLAGS_DEBUG= -g -Wall
CFLAGS_RELEASE= -O2 -fomit-frame-pointer -ffast-math
CFLAGS_SHARED= -fPIC
LIBS= -L../lib -L../json/lib -lm
LIBS_RELEASE= -lonline_structured_svm -ljson
LIBS_DEBUG= -ldonline_structured_svm -ldjson
LD=g++
LDFLAGS= -shared
AR=ar
ARFLAGS=rcs



CPP_FILES := $(wildcard *.cpp)
HEADER_FILES = *.h ../include/*.h ../json/json/*.h 
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

bin/debug_static/test_client.out: test_client.cpp structured_svm_multiclass.cpp $(HEADER_FILES) ../lib/libdonline_structured_svm.a ../json/lib/libdjson.a
	$(CC) $(CFLAGS) $(CFLAGS_DEBUG) -DNO_SERVER $< structured_svm_multiclass.cpp -o $@ $(LIBS)

bin/debug_shared/test_client.out: test_client.cpp structured_svm_multiclass.cpp $(HEADER_FILES) 
	$(CC) $(CFLAGS) $(CFLAGS_DEBUG) $(CFLAGS_RELEASE) -DNO_SERVER $< structured_svm_multiclass.cpp -o $@ $(LIBS)

bin/release_static/test_client.out: test_client.cpp structured_svm_multiclass.cpp $(HEADER_FILES) ../lib/libonline_structured_svm.a ../json/lib/libjson.a
	$(CC) $(CFLAGS) -DNO_SERVER $< structured_svm_multiclass.cpp -o $@ $(LIBS)

bin/release_shared/test_client.out: test_client.cpp structured_svm_multiclass.cpp $(HEADER_FILES) 
	$(CC) $(CFLAGS) $(CFLAGS_RELEASE) -DNO_SERVER $< structured_svm_multiclass.cpp -o $@ $(LIBS)



bin/debug_static/%.out: %.cpp $(HEADER_FILES) ../lib/libdonline_structured_svm.a ../json/lib/libdjson.a 
	$(CC) $(CFLAGS) $(CFLAGS_DEBUG) $< -o $@ $(LIBS)

bin/debug_shared/%.out: %.cpp $(HEADER_FILES) 
	$(CC) $(CFLAGS) $(CFLAGS_DEBUG) $(CFLAGS_RELEASE) $< -o $@ $(LIBS)

bin/release_static/%.out: %.cpp $(HEADER_FILES) ../lib/libonline_structured_svm.a ../json/lib/libjson.a
	$(CC) $(CFLAGS) $< -o $@ $(LIBS)

bin/release_shared/%.out: %.cpp $(HEADER_FILES) 
	$(CC) $(CFLAGS) $(CFLAGS_RELEASE) $< -o $@ $(LIBS)


