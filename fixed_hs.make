# make with make -f OneDHillslope.make

CC=g++
CFLAGS=-c -Wall -O3 -pg
OFLAGS = -Wall -O3 -pg
LDFLAGS= -Wall
SOURCES=DBPR_MCMC.cpp OneDImplicitHillslope.cpp MCMCHillslopeDriver.cpp normal.cpp
OBJECTS=$(SOURCES:.cpp=.o)
# EXECUTABLE=DBPR_MCMC.exe
EXECUTABLE=gauss_DBPR.out

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
