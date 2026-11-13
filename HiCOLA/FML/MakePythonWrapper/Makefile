# Hans A. Winther (hans.a.winther@gmail.com)

SHELL := /bin/bash

#===================================================
# The name of the library and the interface file
#===================================================

PYTHONLIBNAME = _MyLibrary.so
SWIGFILE      = InterfaceFile.i
# These are generated
SWIGWRAPPER   = InterfaceFile_wrap.cxx
SWIGWRAPPERO  = InterfaceFile_wrap.o

#===================================================
# Set c++11 compliant compiler. If USE_MPI we use MPICC 
#===================================================

CC      = g++ -std=c++1z -Wall 

#===================================================
# Options
#===================================================

# Required
USE_PYTHON       = true

#===================================================
# Include and library paths
#===================================================

# Main library include (path to folder containin FML/).
# Default: the parent of this directory, so the wrapper builds in-tree without
# editing. Override on the command line for a separately-installed FML:
#   make FML_INCLUDE=/path/to/FML
FML_INCLUDE    ?= $(abspath $(dir $(lastword $(MAKEFILE_LIST)))/..)

# Python and numpy includes/libs.
# Detected from the active interpreter at build time (no hardcoded framework
# paths). Override any of these on the command line if your setup differs:
#   make PYTHON_INCLUDE=/opt/python/include PYTHON_LIB=/opt/python/lib \
#        NUMPY_INCLUDE=$(python -c "import numpy; print(numpy.get_include())") \
#        PYTHON_LINK="-lpython3.10 -ldl"
PYTHON_BIN     ?= python3
PYTHON_VERSION ?= $(shell $(PYTHON_BIN) -c "import sys; print('python%d.%d' % sys.version_info[:2])")
PYTHON_INCLUDE ?= $(shell $(PYTHON_BIN) -c "import sysconfig; print(sysconfig.get_paths()['include'])")
PYTHON_LIB     ?= $(shell $(PYTHON_BIN) -c "import sysconfig; print(sysconfig.get_config_var('LIBDIR') or '')")
NUMPY_INCLUDE  ?= $(shell $(PYTHON_BIN) -c "import numpy; print(numpy.get_include())")
PYTHON_LINK    ?= -l$(PYTHON_VERSION) -ldl

#===================================================
# Compile up all library defines from options above
#===================================================

INC     = -I$(FML_INCLUDE) 
LIB     =
LINK    = 
OPTIONS = 

ifeq ($(USE_PYTHON),true)
INC  += -I$(PYTHON_INCLUDE) -I$(NUMPY_INCLUDE)
LIB  += -L$(PYTHON_LIB) 
LINK += $(PYTHON_LINK)
endif

#===================================================
# Object files to be compiled
#===================================================

VPATH := 
FILE1 = Main.o
OBJS = $(FILE1)

TARGETS := test
all: $(TARGETS)
.PHONY: all clean

# Compile all the object files as a shared library
$(FILE1) : Main.cpp
	$(CC) $(OPTIONS) -fPIC -c $< -o $@ $(INC)

# Run swig and generate a library that can be called from python
test: $(OBJS)
	swig -c++ -python $(OPTIONS) $(SWIGFILE)
	$(CC) -fPIC -c $(SWIGWRAPPER) $(INC)
	$(CC) -fPIC -shared $(OPTIONS) $(OBJS) $(SWIGWRAPPERO) -o $(PYTHONLIBNAME) $(INC) $(LIB) $(LINK)

clean:
	rm -rf $(TARGETS) *.o $(SWIGWRAPPER) $(PYTHONLIBNAME) *.pyc

%.o: %.cpp 
	${CC} -c -o $@ $< $(OPTIONS) $(INC) 

