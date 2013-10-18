CC = gcc

INCDIR_EXTRA = -I/opt/matlab/extern/include -I/usr/local/MATLAB/R2012a/extern/include
LIBS_EXTRA = -L/opt/matlab/extern/lib -L/opt/matlab/bin/glnxa64 -lmat -lmx -Wl,-rpath,/opt/matlab/bin/glnxa64 -L/usr/local/MATLAB/R2012a/bin/glnxa64 -Wl,-rpath,/usr/local/MATLAB/R2012a/bin/glnxa64

INCDIR = -Iinclude -Ionline_structured_svm/json -Ionline_structured_svm/include 
CFLAGS=$(INCDIR) $(INCDIR_EXTRA) -fopenmp -I/usr/include/opencv
CFLAGS_DEBUG= -g -Wall
CFLAGS_RELEASE= -O2 -fomit-frame-pointer -ffast-math
CFLAGS_SHARED= -fPIC
LIBS= -Llib -Lonline_structured_svm/lib -Lonline_structured_svm/json/lib  -lm $(LIBS_EXTRA)
LIBS_RELEASE= -lonline_structured_svm -ljson
LIBS_DEBUG= -ldonline_structured_svm -ldjson
LD= gcc
LDFLAGS= -shared
AR=ar
ARFLAGS=rcs


CONFIGURATION=debug_static

CPP_FILES := $(wildcard src/*.cpp)
HEADER_FILES := $(wildcard include/*.h)  $(wildcard online_structured_svm/include/*.h) $(wildcard online_structured_svm/json/json/*.h)
OBJ_FILES=$(addprefix obj/$(CONFIGURATION)/,$(notdir $(CPP_FILES:.cpp=.o)))


all: debug_static 

libs:
	cd online_structured_svm; make $(CONFIGURATION)

clean:
	cd online_structured_svm; make clean
	rm -f obj/*/*.out
	rm -f bin/*/*.out

debug_shared: CONFIGURATION=debug_shared
debug_shared: OBJ_FILES=$(addprefix obj/debug_shared/,$(notdir $(CPP_FILES:.cpp=.o)))
debug_shared: LIBS += $(LIBS_DEBUG)
debug_shared: $(addprefix obj/debug_shared/,$(notdir $(CPP_FILES:.cpp=.o))) libs bin/debug_shared/svm_jaaba_behavior_sequence.out

release_shared: CONFIGURATION=release_shared
release_shared: OBJ_FILES=$(addprefix obj/release_shared/,$(notdir $(CPP_FILES:.cpp=.o)))
release_shared: LIBS += $(LIBS_RELEASE)
release_shared: $(addprefix obj/release_shared/,$(notdir $(CPP_FILES:.cpp=.o))) libs bin/release_shared/svm_jaaba_behavior_sequence.out

debug_static: CONFIGURATION=debug_static
debug_static: OBJ_FILES=$(addprefix obj/debug_static/,$(notdir $(CPP_FILES:.cpp=.o)))
debug_static: LIBS += $(LIBS_DEBUG)
debug_static: $(addprefix obj/debug_static/,$(notdir $(CPP_FILES:.cpp=.o))) libs bin/debug_static/svm_jaaba_behavior_sequence.out

release_static: CONFIGURATION=release_static
release_static: OBJ_FILES=$(addprefix obj/release_static/,$(notdir $(CPP_FILES:.cpp=.o)))
release_static: LIBS += $(LIBS_RELEASE)
release_static: $(addprefix obj/release_static/,$(notdir $(CPP_FILES:.cpp=.o))) libs bin/release_static/svm_jaaba_behavior_sequence.out

obj/debug_static/%.o: src/%.cpp $(HEADER_FILES) 
	$(CC) $(CFLAGS) $(CFLAGS_DEBUG) -c $< -o $@

obj/debug_shared/%.o: src/%.cpp $(HEADER_FILES) 
	$(CC) $(CFLAGS) $(CFLAGS_DEBUG) $(CFLAGS_SHARED) -c $< -o $@

obj/release_static/%.o: src/%.cpp $(HEADER_FILES) 
	$(CC) $(CFLAGS) -c $< -o $@

obj/release_shared/%.o: src/%.cpp $(HEADER_FILES) 
	$(CC) $(CFLAGS) $(CFLAGS_SHARED) -c $< -o $@


bin/debug_static/svm_jaaba_behavior_sequence.out: $(HEADER_FILES) online_structured_svm/lib/libdonline_structured_svm.a online_structured_svm/json/lib/libdjson.a $(addprefix obj/debug_static/,$(notdir $(CPP_FILES:.cpp=.o)))
	$(CC) $(CFLAGS) $(CFLAGS_DEBUG) $(OBJ_FILES) -o bin/debug_static/svm_jaaba_behavior_sequence.out $(LIBS)

bin/debug_shared/svm_jaaba_behavior_sequence.out: $(HEADER_FILES) $(addprefix obj/debug_shared/,$(notdir $(CPP_FILES:.cpp=.o)))
	$(CC) $(CFLAGS) $(CFLAGS_DEBUG) $(OBJ_FILES) -o bin/debug_shared/svm_jaaba_behavior_sequence.out $(LIBS)

bin/release_static/svm_jaaba_behavior_sequence.out: $(HEADER_FILES) online_structured_svm/lib/libonline_structured_svm.a online_structured_svm/json/lib/libjson.a $(addprefix obj/release_static/,$(notdir $(CPP_FILES:.cpp=.o)))
	$(CC) $(CFLAGS) $(CFLAGS_RELEASE) $(OBJ_FILES) -o bin/release_static/svm_jaaba_behavior_sequence.out $(LIBS)

bin/release_shared/svm_jaaba_behavior_sequence.out: $(HEADER_FILES) $(addprefix obj/release_shared/,$(notdir $(CPP_FILES:.cpp=.o)))
	$(CC) $(CFLAGS) $(CFLAGS_RELEASE) $(OBJ_FILES) -o bin/release_shared/svm_jaaba_behavior_sequence.out $(LIBS)


