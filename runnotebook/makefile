CC = h5cc -I$(RAPTOR)/src -I$(PWD)
CFLAGS = -fopenmp  -std=c99  -O2 -lm -lgsl -Wall -Wno-unused-but-set-variable
LDFLAGS = -fopenmp -lm -lgsl 

VPATH=$(RAPTOR)/src:
CPATH=$(RAPTOR)/src:

OBJDIR=build

TARGET=RAPTOR

SOURCES=main.c core.c io.c GRmath.c gr_integrator.c rte_integrator.c pol_rte_integrator.c metric.c pol_emission.c tetrad.c model.c constants.c camera.c
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
