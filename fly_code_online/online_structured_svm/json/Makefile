CC = g++
INCDIR = -I.
CFLAGS=$(INCDIR) -fopenmp
CFLAGS_DEBUG= -g -Wall
CFLAGS_RELEASE= -O2 -fomit-frame-pointer -ffast-math
CFLAGS_SHARED= -fPIC
LD= g++
LDFLAGS= -shared
AR=ar
ARFLAGS=rcs


CONFIGURATION=debug_static

CPP_FILES := $(wildcard src/*.cpp)
HEADER_FILES := $(wildcard json/json/*.h)
OBJ_FILES := $(addprefix obj/$(CONFIGURATION)/,$(notdir $(CPP_FILES:.cpp=.o)))

all: $(CONFIGURATION)

doc: *.h src/*.cpp json/*.h doxygen.config
	doxygen doxygen.config

debug_shared: OBJ_FILES=$(addprefix obj/debug_shared/,$(notdir $(CPP_FILES:.cpp=.o)))
debug_shared: $(addprefix obj/debug_shared/,$(notdir $(CPP_FILES:.cpp=.o))) lib/libdjson.so

release_shared: OBJ_FILES=$(addprefix obj/release_shared/,$(notdir $(CPP_FILES:.cpp=.o)))
release_shared: $(addprefix obj/release_shared/,$(notdir $(CPP_FILES:.cpp=.o))) lib/libjson.so

debug_static: OBJ_FILES=$(addprefix obj/debug_static/,$(notdir $(CPP_FILES:.cpp=.o)))
debug_static: $(addprefix obj/debug_static/,$(notdir $(CPP_FILES:.cpp=.o))) lib/libdjson.a

release_static: OBJ_FILES=$(addprefix obj/release_static/,$(notdir $(CPP_FILES:.cpp=.o)))
release_static: $(addprefix obj/release_static/,$(notdir $(CPP_FILES:.cpp=.o))) lib/libjson.a

clean: 
	rm -f obj/*/*.o lib/*


obj/debug_static/%.o: src/%.cpp $(HEADER_FILES) 
	$(CC) $(CFLAGS) $(CFLAGS_DEBUG) -c $< -o $@

obj/debug_shared/%.o: src/%.cpp $(HEADER_FILES) 
	$(CC) $(CFLAGS) $(CFLAGS_DEBUG) $(CFLAGS_SHARED) -c $< -o $@

obj/release_static/%.o: src/%.cpp $(HEADER_FILES) 
	$(CC) $(CFLAGS) -c $< -o $@

obj/release_shared/%.o: src/%.cpp $(HEADER_FILES) 
	$(CC) $(CFLAGS) $(CFLAGS_SHARED) -c $< -o $@


lib/libjson.so:  $(OBJ_FILES)
	$(LD) $(LDFLAGS) $(OBJ_FILES) -o lib/libjson.so

lib/libjson.a:  $(OBJ_FILES)
	$(AR) $(ARFLAGS) lib/libjson.a $(OBJ_FILES) 

lib/libdjson.so:  $(OBJ_FILES) 
	$(LD) $(LDFLAGS) $(OBJ_FILES) -o lib/libdjson.so

lib/libdjson.a:  $(OBJ_FILES) 
	$(AR) $(ARFLAGS) lib/libdjson.a $(OBJ_FILES) 


