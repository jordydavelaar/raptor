
CC = gcc
CFLAGS = -g -openmp -std=c99 -I/usr/include
LDFLAGS = -lm -lgsl -lcblas

SRC=orbit_tracer.c GRmath.c integrator.c metric.c plasma.c radiative_transfer.c RCARRY.c
OBJ=orbit_tracer.o GRmath.o integrator.o metric.o plasma.o radiative_transfer.o RCARRY.o

SRC1=orbit_tracer.c GRmath.c integrator.c metric.c plasma.c radiative_transfer.c 
OBJ1=orbit_tracer.o GRmath.o integrator.o metric.o plasma.o radiative_transfer.o

SRC2=img_renderer.c GRmath.c integrator.c metric.c plasma.c radiative_transfer.c raptor_harm_model.c utilities.c RCARRY.c
OBJ2=img_renderer.o GRmath.o integrator.o metric.o plasma.o radiative_transfer.o raptor_harm_model.o utilities.o RCARRY.o

run: $(OBJ) makefile
	$(CC) $(CFLAGS) -o RAPTOR $(OBJ) $(LDFLAGS)

orbit: $(OBJ1) makefile
	$(CC) $(CFLAGS) -o RAPTOR $(OBJ1) $(LDFLAGS)

img: $(OBJ2) makefile
	$(CC) $(CFLAGS) -o RAPTOR $(OBJ2) $(LDFLAGS)

$(OBJ): makefile functions.h constants.h parameters.h raptor_harm_model.h

clean:
	rm *.o
