CC = h5cc -I$(RAPTOR)/src -I$(PWD)
<<<<<<< HEAD
CFLAGS = -fopenmp  -std=c99  -lm -lgsl -Wall -Wno-unused-but-set-variable
=======
CFLAGS = -fopenmp  -std=c99  -lm -lgsl -Wall
>>>>>>> c427bcca7e35ae8c365b572e0bce03f433112ccb
LDFLAGS = -fopenmp -lm -lgsl

VPATH=$(RAPTOR)/src:
CPATH=$(RAPTOR)/src:

OBJDIR=build

TARGET=RAPTOR

<<<<<<< HEAD
SOURCES=main.c core.c io.c GRmath.c gr_integrator.c rte_integrator.c pol_rte_integrator.c metric.c pol_emission.c tetrad.c raptor_model.c utilities.c camera.c
=======
SOURCES=main.c core.c io.c GRmath.c gr_integrator.c rte_integrator.c pol_rte_integrator.c metric.c emission.c pol_emission.c tetrad.c raptor_model.c utilities.c camera.c
>>>>>>> c427bcca7e35ae8c365b572e0bce03f433112ccb
OBJECTS := $(patsubst %.c,$(OBJDIR)/%.o,$(SOURCES))

all: create_directories $(SOURCES) $(TARGET)

$(TARGET): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

create_directories:
	@test -d $(OBJDIR) || mkdir -v $(OBJDIR)


$(OBJECTS): $(OBJDIR)/%.o: %.c
	$(CC) $(CFLAGS)  -c $^ -o $@


clean:
	rm -rf $(OBJECTS) $(TARGET)
