# Makefile for SVM-struct, 03.10.06

#Use the following to compile under unix or cygwin
CC = g++
LD = g++

#Call 'make' using the following line to make CYGWIN produce stand-alone Windows executables
#		make 'SFLAGS=-mno-cygwin'


CFLAGS =    $(SFLAGS) -g -fomit-frame-pointer -ffast-math -Wall $(ARCH) -fopenmp -DUSE_OPENMP 
LDFLAGS =   $(SFLAGS) -g -lm -Wall -fopenmp

#CFLAGS =   $(SFLAGS) -pg -Wall
#LDFLAGS =  $(SFLAGS) -pg -Wall

all: svm_struct_noexe

svm_struct_noexe: svm_struct_learn.o svm_struct_classify.o svm_struct_common.o svm_struct_main.o 

.PHONY: clean
clean:
	rm -f *.o *.tcov *.d core gmon.out *.stackdump


#----------------------#
#----  STRUCT SVM  ----#
#----------------------#

svm_struct_common.o: svm_struct_common.cpp svm_struct_common.h ../svm_struct_api_types.h ../svm_struct_api.h ../svm_struct_api_multiclass.h  ../svm_struct_api_behavior_sequence.h ../svm_struct_api_blob_behavior_sequence.h 
	$(CC) -c $(CFLAGS) svm_struct_common.cpp -o svm_struct_common.o

svm_struct_learn.o: svm_struct_learn.cpp ../svm_light/svm_learn.h svm_struct_common.h ../svm_struct_api.h ../svm_struct_api_types.h ../svm_struct_api_multiclass.h  ../svm_struct_api_behavior_sequence.h  ../svm_struct_api_blob_behavior_sequence.h 
	$(CC) -c $(CFLAGS) svm_struct_learn.cpp -o svm_struct_learn.o

svm_struct_main.o: svm_struct_main.cpp ../svm_light/svm_common.h ../svm_light/svm_learn.h svm_struct_learn.h svm_struct_common.h ../svm_struct_api.h ../svm_struct_api_types.h ../svm_struct_api_multiclass.h  ../svm_struct_api_behavior_sequence.h  ../svm_struct_api_blob_behavior_sequence.h 
	$(CC) -c $(CFLAGS) svm_struct_main.cpp -o svm_struct_main.o

svm_struct_classify.o: svm_struct_classify.cpp svm_struct_common.h ../svm_struct_api_types.h ../svm_struct_api.h ../svm_light/svm_common.h  ../svm_struct_api_multiclass.h  ../svm_struct_api_behavior_sequence.h  ../svm_struct_api_blob_behavior_sequence.h 
	$(CC) -c $(CFLAGS) svm_struct_classify.cpp -o svm_struct_classify.o
