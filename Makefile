CC = g++
CXXFLAGS ?= -std=c++20
CPPFLAGS ?= -O3 -I. #-Wall 
PACS_ROOT= ../../Examples/pacs-examples/Examples/lib
LDFLAGS ?= -L$(PACS_ROOT)
LIBS ?= 

SOURCES = main.cpp # algebra.cpp 
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
