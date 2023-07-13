# Build executable with:
# % make
# Delete object files and executable with:
# % make clean
# Rebuild all objects and executable with:
# % make -B

SRC_DIR := src
OBJ_DIR := build
BIN_DIR := bin
ARMA_DIR := /home/wangd17/armadillo-10.2.2

EXECUTABLE := run_swdca
SOURCES := $(wildcard $(SRC_DIR)/*.cpp)
HEADERS := $(wildcard $(SRC_DIR)/*.h)
CXX := g++ #g++-9

SHELL = /bin/sh

# Flags to pass to the compiler; per the reccomendations of the GNU Scientific Library
CXXFLAGS:= -std=c++17 -Wextra -pedantic -Wall -W -Wmissing-declarations -Wuninitialized -Wshadow -Wpointer-arith -Wcast-align -Wwrite-strings -fshort-enums -fno-common -m64 -fopenmp -I$(ARMA_DIR)/include #-I$(HOME)/include

# Compiler flags controling optimization levels. Use -O3 for full optimization,
# but make sure your results are consistent
# -g includes debugging information. You can also add -pg here for profiling 
PROFILE=-pg
OPTFLAGS:=$(PROFILE) -O2

# Flags to pass to the linker; -lm links in the standard c math library
LDFLAGS:= -fopenmp -lm -lgsl -lgslcblas -llapack -lblas $(PROFILE) -L$(ARMA_DIR) #-L$(HOME)/lib

# Variable to compose names of object files from the names of sources
OBJECTS := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SOURCES))

# Default target depends on sources and headers to detect changes
all: $(SOURCES) $(HEADERS)  $(BIN_DIR)/$(EXECUTABLE)

# Rule to compile a source file to object code
$(OBJECTS): $(OBJ_DIR)/%.o : $(SRC_DIR)/%.cpp
	$(CXX) -c $(CXXFLAGS) $(OPTFLAGS) $< -o $@

# Build the executable by linking all objects
$(BIN_DIR)/$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(OBJECTS) $(LDFLAGS) -o $@

# clean up so we can start over (removes executable!)
clean:
	rm -f $(OBJ_DIR)/*.o $(BIN_DIR)/$(EXECUTABLE)
