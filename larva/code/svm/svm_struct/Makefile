# Makefile for SVM-multiclass, 03.07.04

#Use the following to compile under unix or cygwin
CC = g++
LD = g++

#Call 'make' using the following line to make CYGWIN produce stand-alone Windows executables
#		make 'SFLAGS=-mno-cygwin'

CFLAGS =   $(SFLAGS) -g -fomit-frame-pointer -ffast-math -Wall $(ARCH)
LDFLAGS =  $(SFLAGS) -g -lm -Wall $(ARCH)
#CFLAGS =  $(SFLAGS) -pg -Wall
#LDFLAGS = $(SFLAGS) -pg -lm -Wall 

all: libsvmstructmulti.a 


.PHONY: clean
clean: svm_light_clean svm_struct_clean
	rm -f *.o *.a *.tcov *.d core gmon.out *.stackdump 

#-----------------------#
#----   SVM-light   ----#
#-----------------------#
svm_light_hideo_noexe: 
	cd svm_light; make svm_learn_hideo_noexe

svm_light_clean: 
	cd svm_light; make clean

#----------------------#
#----  STRUCT SVM  ----#
#----------------------#

svm_struct_noexe: 
	cd svm_struct; make svm_struct_noexe

svm_struct_clean: 
	cd svm_struct; make clean

#-------------------------#
#----  SVM MULTICLASS ----#
#-------------------------#

libsvmstructmulti.a: svm_light_hideo_noexe svm_struct_noexe svm_struct_api_multiclass.o svm_struct_api_behavior_sequence.o svm_struct_api_blob_behavior_sequence.o svm_struct_learn_custom.o svm_struct/svm_struct_learn.o svm_struct/svm_struct_common.o svm_light/svm_learn.o svm_light/svm_common.o svm_light/svm_hideo.o svm_struct/svm_struct_common.o 
	ar rcs libsvmstructmulti.a svm_struct_api_multiclass.o svm_struct_api_behavior_sequence.o svm_struct_api_blob_behavior_sequence.o svm_struct_learn_custom.o svm_struct/svm_struct_learn.o svm_struct/svm_struct_common.o svm_light/svm_learn.o svm_light/svm_common.o svm_light/svm_hideo.o svm_struct/svm_struct_common.o 

svm_struct_api_multiclass.o: svm_struct_api_multiclass.cpp svm_struct_api.h svm_struct_api_multiclass.h svm_struct_api_types.h svm_struct/svm_struct_common.h
	$(CC) -c $(CFLAGS) svm_struct_api_multiclass.cpp -o svm_struct_api_multiclass.o

svm_struct_api_behavior_sequence.o: svm_struct_api_behavior_sequence.cpp svm_struct_api.h svm_struct_api_behavior_sequence.h svm_struct_api_types.h svm_struct/svm_struct_common.h
	$(CC) -c $(CFLAGS) svm_struct_api_behavior_sequence.cpp -o svm_struct_api_behavior_sequence.o

svm_struct_api_blob_behavior_sequence.o: svm_struct_api_blob_behavior_sequence.cpp svm_struct_api.h svm_struct_api_blob_behavior_sequence.h svm_struct_api_behavior_sequence.h svm_struct_api_types.h svm_struct/svm_struct_common.h
	$(CC) -c $(CFLAGS) svm_struct_api_blob_behavior_sequence.cpp -o svm_struct_api_blob_behavior_sequence.o

svm_struct_learn_custom.o: svm_struct_learn_custom.cpp svm_struct_api.h svm_light/svm_common.h svm_struct_api_types.h svm_struct/svm_struct_common.h
	$(CC) -c $(CFLAGS) svm_struct_learn_custom.cpp -o svm_struct_learn_custom.o
