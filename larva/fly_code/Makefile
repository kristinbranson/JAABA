CC = g++
LD = g++
UNAME := $(shell uname -s)
WXCONFIG = wx-config

ARCH = 
COPTS = -g -Wall -fomit-frame-pointer -ffast-math 
LIBS = -lcv -lhighgui -lml -lcxcore
LOPTS = -g -Wall
CFLAGS = $(COPTS) -I/usr/include/opencv
LFLAGS = $(LOPTS) $(LIBS) -fopenmp
EXTRA_PROGS = server client

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

#all: gui svm_struct_classify svm_struct_learn $(EXTRA_PROGS)
all: svm_struct_classify svm_struct_learn $(EXTRA_PROGS)

libsvmstructmulti.a: 
	cd svm/svm_struct; make

libjsonrpc.a:
	cd behaviors/json; make

clean:
	rm -f gui svm_struct_classify svm_struct_learn *.o; cd svm/svm_struct; make clean

mac:
	cp -f gui ../mac/

blob.o: blob.cpp blob.h
	$(CC) $(CFLAGS) -c blob.cpp 

fit_model.o: fit_model.cpp fit_model.h blob.h spine_features.h
	$(CC) $(CFLAGS) -c fit_model.cpp 

train.o: train.cpp fit_model.h blob.h spine_features.h train.h
	$(CC) $(CFLAGS) -c train.cpp 

spine_features.o: spine_features.cpp fit_model.h blob.h spine_features.h
	$(CC) $(CFLAGS) -c spine_features.cpp 

svgPlotter.o: behavior/svgPlotter.cpp blob.h 
	$(CC) $(CFLAGS) -c behavior/svgPlotter.cpp 

server.o: behavior/server.cpp fit_model.h blob.h spine_features.h behavior/svgPlotter.h libsvmstructmulti.a libjsonrpc.a
	$(CC) $(CFLAGS) -c behavior/server.cpp 

client.o: behavior/client.cpp blob.h libjsonrpc.a
	$(CC) $(CFLAGS) -c behavior/client.cpp 

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

server: libsvmstructmulti.a libbehaviors.a libjsonrpc.a server.o 
	$(LD) $(LFLAGS) -o server  svgPlotter.o server.o svm/svm_struct/libsvmstructmulti.a libbehaviors.a libjsonrpc.a 

client: libjsonrpc.a client.o 
	$(LD) $(LFLAGS) -o client client.o libjsonrpc.a 

libjsonrpc.a: behavior/json/*.cpp behavior/json/*.h
	cd behavior/json/; $(CC) $(CFLAGS) -c *.cpp -I.; ar rcs ../../libjsonrpc.a *.o