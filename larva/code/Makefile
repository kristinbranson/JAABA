CC = g++
LD = g++
UNAME := $(shell uname -s)
WXCONFIG = wx-config

ARCH = 
COPTS = -g -Wall -fomit-frame-pointer -ffast-math 
LIBS = -lcv -lhighgui -lml -lcxcore
LOPTS = -g -Wall -fopenmp
CFLAGS = $(COPTS) -I/usr/include/opencv -DUSE_OPENMP -fopenmp
LFLAGS = $(LOPTS) $(LIBS) -fopenmp
EXTRA_PROGS =

ifeq ($(UNAME),Darwin)
  export ARCH = -arch i386 -fopenmp
  CFLAGS = $(COPTS) $(ARCH) -I../mac_local/include/opencv 
  WXCONFIG = ../mac_local/bin/wx-config
  LFLAGS = $(ARCH) $(LOPTS) $(LIBS)  -L../mac_local/lib #../mac_local/lib/libiconv.a 
  EXTRA_PROGS = mac
endif

WX_HEADERS = $(shell $(WXCONFIG) --cxxflags)
WX_LIBS = $(shell $(WXCONFIG) --libs)
O_FILES = blob.o fit_model.o spine_features.o train.o

all: gui svm_struct_classify svm_struct_learn $(EXTRA_PROGS)
#all: svm_struct_classify svm_struct_learn $(EXTRA_PROGS)

libsvmstructmulti.a: 
	cd svm/svm_struct; make

clean:
	rm -f gui svm_struct_classify svm_struct_learn *.o; cd svm/svm_struct; make clean

mac:
	cp -f gui ../mac/

#test.o: test.cpp
#	$(CC) $(CFLAGS) -c test.cpp
#
#test: test.o
#	$(CC) $(LFLAGS) $(WX_LIBS) -o test test.o

blob.o: blob.cpp blob.h
	$(CC) $(CFLAGS) -c blob.cpp 

fit_model.o: fit_model.cpp fit_model.h blob.h spine_features.h
	$(CC) $(CFLAGS) -c fit_model.cpp 

spine_features.o: spine_features.cpp fit_model.h blob.h spine_features.h
	$(CC) $(CFLAGS) -c spine_features.cpp 

train.o: train.cpp fit_model.h blob.h spine_features.h train.h
	$(CC) $(CFLAGS) -c train.cpp 

libbehaviors.a: $(O_FILES)
	ar rcs libbehaviors.a $(O_FILES)

gui.o: gui.cpp fit_model.h gui.h blob.h spine_features.h train.h
	$(CC) $(CFLAGS) $(WX_HEADERS) -c gui.cpp 

gui: $(O_FILES) gui.o libsvmstructmulti.a libbehaviors.a
	$(CC) $(LFLAGS) $(WX_LIBS) -o gui $(O_FILES) gui.o svm/svm_struct/libsvmstructmulti.a

svm_struct_classify: libsvmstructmulti.a libbehaviors.a svm/svm_struct/svm_struct/svm_struct_classify.o 
	$(LD) $(LFLAGS) -o svm_struct_classify svm/svm_struct/svm_struct/svm_struct_classify.o svm/svm_struct/libsvmstructmulti.a libbehaviors.a

svm_struct_learn: libsvmstructmulti.a libbehaviors.a svm/svm_struct/svm_struct/svm_struct_main.o 
	$(LD) $(LFLAGS) -o svm_struct_learn svm/svm_struct/svm_struct/svm_struct_main.o svm/svm_struct/libsvmstructmulti.a libbehaviors.a
