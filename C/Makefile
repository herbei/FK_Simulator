# CONFIG

#compiler
CC = gcc
GPUCC = /usr/local/cuda/9.2.88/bin/nvcc



FILE	= fk_simulator


OBJECTS  = $(FILE).o fk_functions.o
#PROGS    = $(FILE).e



INCLUDE_FLAGS = -I/usr/local/include -I/usr/local/cuda/9.2.88/include

LINK_FLAGS = -L/usr/local/lib -lgsl -lgslcblas -lm -L/usr/local/cuda/9.2.88/lib64 -lcuda -lcudart -lcurand

all:            	fk_run

clean:
			rm -f *.o
clear:			
			rm *.*~


o:		
			open $(FILE).c

fk_run:			$(FILE).o fk_functions.o
			$(CC) $(LINK_FLAGS)  -o fk_run $(FILE).o fk_functions.o

fk_functions.o:		fk_functions.cu
			$(GPUCC) --gpu-architecture=sm_60 -c -o fk_functions.o fk_functions.cu





$(FILE).o:		$(FILE).c
			$(CC) -Wall $(INCLUDE_FLAGS)  -c -o $(FILE).o $(FILE).c

