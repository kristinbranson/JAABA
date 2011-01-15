# Makefile for Latent Structural SVM

CC=gcc
#CFLAGS= -g -Wall
CFLAGS= -O3 -fomit-frame-pointer -ffast-math
#CFLAGS = -O3 -pg
LD=gcc
#LDFLAGS= -g
LDFLAGS= -O3
#LDFLAGS = -O3 -pg
LIBS= -lm

all: svm_latent_learn svm_latent_classify

clean: 
	rm -f *.o
	rm -f svm_latent_learn svm_latent_classify

svm_light_hideo_noexe: 
	cd svm_light; make svm_learn_hideo_noexe

svm_latent_learn: svm_light_hideo_noexe svm_struct_latent_cccp.o svm_struct_latent_api.o 
	$(LD) $(LDFLAGS) svm_struct_latent_cccp.o svm_light/svm_hideo.o svm_light/svm_learn.o svm_light/svm_common.o svm_struct_latent_api.o -o svm_latent_learn -lm

svm_latent_classify: svm_light_hideo_noexe svm_struct_latent_classify.o svm_struct_latent_api.o 
	$(LD) $(LDFLAGS) svm_struct_latent_classify.o svm_light/svm_common.o svm_struct_latent_api.o -o svm_latent_classify $(LIBS)

svm_struct_latent_cccp.o: svm_struct_latent_cccp.c
	$(CC) -c $(CFLAGS) svm_struct_latent_cccp.c -o svm_struct_latent_cccp.o

svm_struct_latent_api.o: svm_struct_latent_api.c svm_struct_latent_api_types.h
	$(CC) -c $(CFLAGS) svm_struct_latent_api.c -o svm_struct_latent_api.o

svm_struct_latent_classify.o: svm_struct_latent_classify.c
	$(CC) -c $(CFLAGS) svm_struct_latent_classify.c -o svm_struct_latent_classify.o
