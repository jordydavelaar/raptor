
CC = h5cc
CFLAGS = -openmp -lomp -std=c99 -I/usr/include -Wall
LDFLAGS = -lm -lgsl -lcblas -g

SRC=orbit_tracer.c GRmath.c integrator.c metric.c plasma.c radiative_transfer.c
OBJ=orbit_tracer.o GRmath.o integrator.o metric.o plasma.o radiative_transfer.o

SRC1=orbit_tracer.c GRmath.c integrator.c metric.c plasma.c radiative_transfer.c
OBJ1=orbit_tracer.o GRmath.o integrator.o metric.o plasma.o radiative_transfer.o

SRC2=main.c core.c io.c GRmath.c integrator.c metric.c emission.c pol_emission.c tetrad.c raptor_harm3d_model.c utilities.c
OBJ2=main.o core.o io.o GRmath.o integrator.o metric.o emission.o pol_emission.o tetrad.o raptor_harm3d_model.o utilities.o

run: $(OBJ) makefile
	$(CC) $(CFLAGS) -o RAPTOR $(OBJ) $(LDFLAGS)

orbit: $(OBJ1) makefile
	$(CC) $(CFLAGS) -o RAPTOR $(OBJ1) $(LDFLAGS)

img: $(OBJ2) makefile
	$(CC) $(CFLAGS) -o RAPTOR $(OBJ2) $(LDFLAGS)

$(OBJ): makefile functions.h constants.h parameters.h raptor_harm3d_model.h

clean:
	rm *.o
	make img
