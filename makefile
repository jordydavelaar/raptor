
CC = h5cc
CFLAGS = -fopenmp  -std=c99  -lm -lgsl -Wall
LDFLAGS = -lm -lgsl

SRC2=main.c core.c io.c GRmath.c gr_integrator.c rte_integrator.c pol_rte_integrator.c metric.c emission.c pol_emission.c tetrad.c raptor_harm3d_model.c utilities.c camera.c
OBJ2=main.o core.o io.o GRmath.o gr_integrator.o rte_integrator.o pol_rte_integrator.o metric.o emission.o pol_emission.o tetrad.o raptor_harm3d_model.o utilities.o camera.o

SRC3=main.c core.c io.c GRmath.c gr_integrator.c rte_integrator.c pol_rte_integrator.c metric.c emission.c pol_emission.c tetrad.c raptor_bhac3d_model.c utilities.c camera.c
OBJ3=main.o core.o io.o GRmath.o gr_integrator.o rte_integrator.o pol_rte_integrator.o metric.o emission.o pol_emission.o tetrad.o raptor_bhac3d_model.o utilities.o camera.o

harm: $(OBJ2) makefile
	$(CC) $(CFLAGS) -o RAPTOR $(OBJ2) $(LDFLAGS)

bhac: $(OBJ3) makefile
	$(CC) $(CFLAGS) -o RAPTOR $(OBJ3) $(LDFLAGS)

$(OBJ): makefile functions.h parameters.h raptor_bhac3d_model.h

clean:
	rm *.o
