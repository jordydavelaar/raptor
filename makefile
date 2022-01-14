
CC = h5cc
CFLAGS = -openmp -lomp -std=c99 -I/usr/include -Wall
LDFLAGS = -lm -lgsl -lcblas -g

SRC2=main.c core.c io.c GRmath.c gr_integrator.c rte_integrator.c pol_rte_integrator.c metric.c emission.c pol_emission.c tetrad.c raptor_harm3d_model.c utilities.c camera.c
OBJ2=main.o core.o io.o GRmath.o gr_integrator.o rte_integrator.o pol_rte_integrator.o metric.o emission.o pol_emission.o tetrad.o raptor_harm3d_model.o utilities.o camera.o

img: $(OBJ2) makefile
	$(CC) $(CFLAGS) -o RAPTOR $(OBJ2) $(LDFLAGS)

$(OBJ): makefile functions.h parameters.h raptor_harm3d_model.h

clean:
	rm *.o
	make img
