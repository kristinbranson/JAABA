CC = g++
INCDIR = -Iinclude -I../online_structured_svm/json -I../online_structured_svm/include -I../blob/include
CFLAGS=$(INCDIR) -fopenmp -I/usr/include/opencv
CFLAGS_DEBUG= -g -Wall
CFLAGS_RELEASE= -O2 -fomit-frame-pointer -ffast-math
CFLAGS_SHARED= -fPIC
LD= g++
LDFLAGS= -shared
AR=ar
ARFLAGS=rcs


CONFIGURATION=debug_static

CPP_FILES := $(wildcard src/*.cpp)
HEADER_FILES := $(wildcard include/*.h)  
OBJ_FILES=$(addprefix obj/$(CONFIGURATION)/,$(notdir $(CPP_FILES:.cpp=.o)))

all: $(CONFIGURATION)

example:
	cd examples; make $(CONFIGURATION)

debug_shared: CONFIGURATION=debug_shared
debug_shared: OBJ_FILES=$(addprefix obj/debug_shared/,$(notdir $(CPP_FILES:.cpp=.o)))
debug_shared: $(addprefix obj/debug_shared/,$(notdir $(CPP_FILES:.cpp=.o))) lib/libdbehavior.so example

release_shared: CONFIGURATION=release_shared
release_shared: OBJ_FILES=$(addprefix obj/release_shared/,$(notdir $(CPP_FILES:.cpp=.o)))
release_shared: $(addprefix obj/release_shared/,$(notdir $(CPP_FILES:.cpp=.o))) lib/libbehavior.so example

debug_static: CONFIGURATION=debug_static
debug_static: OBJ_FILES=$(addprefix obj/debug_static/,$(notdir $(CPP_FILES:.cpp=.o)))
debug_static: $(addprefix obj/debug_static/,$(notdir $(CPP_FILES:.cpp=.o))) lib/libdbehavior.a example

release_static: CONFIGURATION=release_static
release_static: OBJ_FILES=$(addprefix obj/release_static/,$(notdir $(CPP_FILES:.cpp=.o)))
release_static: $(addprefix obj/release_static/,$(notdir $(CPP_FILES:.cpp=.o))) lib/libbehavior.a example

clean: 
	rm -f *.o obj/*/*.o lib/*
	cd examples; make clean

obj/debug_static/%.o: src/%.cpp $(HEADER_FILES) 
	$(CC) $(CFLAGS) $(CFLAGS_DEBUG) -c $< -o $@

obj/debug_shared/%.o: src/%.cpp $(HEADER_FILES) 
	$(CC) $(CFLAGS) $(CFLAGS_DEBUG) $(CFLAGS_SHARED) -c $< -o $@

obj/release_static/%.o: src/%.cpp $(HEADER_FILES) 
	$(CC) $(CFLAGS) -c $< -o $@

obj/release_shared/%.o: src/%.cpp $(HEADER_FILES) 
	$(CC) $(CFLAGS) $(CFLAGS_SHARED) -c $< -o $@


lib/libbehavior.so: $(OBJ_FILES) 
	$(LD) $(LDFLAGS) $(OBJ_FILES) -o lib/libbehavior.so

lib/libbehavior.a: $(OBJ_FILES) 
	$(AR) $(ARFLAGS) lib/libbehavior.a $(OBJ_FILES) 

lib/libdbehavior.so: $(OBJ_FILES) 
	$(LD) $(LDFLAGS) $(OBJ_FILES) -o lib/libdbehavior.so

lib/libdbehavior.a: $(OBJ_FILES) 
	$(AR) $(ARFLAGS) lib/libdbehavior.a $(OBJ_FILES) 


