
CC = h5cc -Isrc -Imodel
CFLAGS = -fopenmp  -std=c99  -lm -lgsl -Wall
LDFLAGS = -lm -lgsl

VPATH=src
CPATH=$(RAPTOR)/src

SRC=main.c core.c io.c GRmath.c gr_integrator.c rte_integrator.c pol_rte_integrator.c metric.c emission.c pol_emission.c tetrad.c raptor_model.c utilities.c camera.c
OBJ=main.o core.o io.o GRmath.o gr_integrator.o rte_integrator.o pol_rte_integrator.o metric.o emission.o pol_emission.o tetrad.o raptor_model.o utilities.o camera.o

grmhd: $(OBJ) makefile
	$(CC) $(CFLAGS) -o RAPTOR $(OBJ) $(LDFLAGS)


$(OBJ): makefile functions.h parameters.h raptor_model.h

clean:
	rm *.o
	rm RAPTOR
