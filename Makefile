# ==============================================================================
# FEATURE TOGGLES (Comment out or change to 0/1 to toggle features)
# ==============================================================================
# Set to 1 for debugging mode (-O0 -g), or 0 for optimized release (-O3):
DEBUG := 0

# Comment out this line if you do not want OpenVDB support:
USE_VDB := 1

# Comment out this line if your compiler doesn't support OpenMP parallelization:
USE_OMP := 1
# ==============================================================================

# Compiler definition
CXX := g++

# Base compiler flags (always used)
CXXFLAGS := -Wall -Wextra
LDFLAGS  := 
LIBS     := 

# ------------------------------------------------------------------------------
# Debug vs Release Flag Management
# ------------------------------------------------------------------------------
ifeq ($(DEBUG),1)
    CXXFLAGS += -O0 -g -DDEBUG
else
    CXXFLAGS += -O3
endif

# ------------------------------------------------------------------------------
# Conditional OpenMP Setup
# ------------------------------------------------------------------------------
ifdef USE_OMP
    CXXFLAGS += -fopenmp -DUSE_OMP
    LDFLAGS  += -fopenmp
endif

# ------------------------------------------------------------------------------
# Conditional OpenVDB Setup
# ------------------------------------------------------------------------------
ifdef USE_VDB
    LOCAL_DIR := $(HOME)/local
    
    # OpenVDB requires at least C++17
    CXXFLAGS += -std=c++17 -I$(LOCAL_DIR)/include -DUSE_VDB
    LDFLAGS  += -L$(LOCAL_DIR)/lib -L$(LOCAL_DIR)/lib64
    LIBS     += -lopenvdb -ltbb -lblosc -lz
endif

# ------------------------------------------------------------------------------
# Build Targets and Logic
# ------------------------------------------------------------------------------
SRCS := main.cpp EquationOfState.cpp Euler.cpp FluxSolver.cpp FileHandler.cpp \
        Mesh.cpp Reconstruction.cpp Solver.cpp STLReader.cpp TestProblems.cpp
OBJS := $(SRCS:.cpp=.o)
TARGET := ilmatar-cfd

all: $(TARGET)

# Link step
$(TARGET): $(OBJS)
	$(CXX) $(LDFLAGS) $^ -o $@ $(LIBS)

# Compile step
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)

# ------------------------------------------------------------------------------
# C++ unit tests (Catch2 single-header)
# ------------------------------------------------------------------------------
UNIT_DEPS := EquationOfState.o Euler.o FluxSolver.o FileHandler.o \
             Mesh.o Reconstruction.o STLReader.o
UNIT_SRCS := tests/unit/test_main.cpp \
             tests/unit/test_eos.cpp \
             tests/unit/test_flux.cpp \
             tests/unit/test_reconstruction.cpp \
             tests/unit/test_mesh.cpp
UNIT_OBJS := $(UNIT_SRCS:.cpp=.o)
UNIT_BIN  := tests/unit/run_tests

unit_tests: $(UNIT_BIN)

$(UNIT_BIN): $(UNIT_OBJS) $(UNIT_DEPS)
	$(CXX) $(LDFLAGS) $^ -o $@ $(LIBS)

tests/unit/%.o: tests/unit/%.cpp
	$(CXX) $(CXXFLAGS) -I. -c $< -o $@

run_unit_tests: unit_tests
	./$(UNIT_BIN)

clean_unit:
	rm -f $(UNIT_OBJS) $(UNIT_BIN)

.PHONY: all clean unit_tests run_unit_tests clean_unit
