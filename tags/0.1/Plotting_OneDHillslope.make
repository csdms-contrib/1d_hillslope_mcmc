# make with make -f Plotting_OneDHillslope.make

CC=g++
CFLAGS=-c -Wall -O3 -pg
OFLAGS = -Wall -O3 -pg
LDFLAGS= -Wall
SOURCES=Es_Rs_Curve_Files.cpp OneDImplicitHillslope.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=EsRs_plotting.exe
# EXECUTABLE=EsRs_plotting.out

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
