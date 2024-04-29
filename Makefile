CC = g++
CXXFLAGS ?= -std=c++20
# PACS_ROOT= ../../Examples/pacs-examples/Examples
CPPFLAGS ?= -O3 -Wall -I. -I$(PACS_ROOT)/src/Utilities
LDFLAGS ?= 
LIBS ?= 

SOURCES = main.cpp 
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = main
HEADERS = algebra.hpp

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(CXXFLAGS) $(LDFLAGS) $(LIBS) -o $@ $(OBJECTS)

%.o: %.cpp
	$(CC) $(CPPFLAGS) -c -o $@ $<

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)
