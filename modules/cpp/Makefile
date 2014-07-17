# Makefile for C++ module.  The source C++ (.cpp) files are compiled
# to object files (.o).  These are linked with the Boost.Python
# library, which creates a shared object library (.so) in this
# directory that can be imported into Python.  This process requires
# Boost.Python wrapper code for the functions in the shared library,
# which is contained in the source file ops/$(LIBNAME).cpp. The name
# of the shared library that is built is $(LIBNAME).so, where LIBNAME
# is defined below.

PYTHON_VERSION := $(shell python -c \
                   'import sys; print str(sys.version_info[0]) + "." + str(sys.version_info[1])')
PYTHON_INCLUDE = /usr/include/python$(PYTHON_VERSION)

# location of the Boost Python include files and library
BOOST_INC = /usr/include
BOOST_LIB = /usr/lib

# compiler options
CXX = g++
CXXFLAGS = -O3 -fPIC -I$(PYTHON_INCLUDE) -I$(BOOST_INC)

SRCDIR = ops
OBJDIR = ops

# the name of the shared library that Python can import
LIBNAME = op_ext

# the Boost.Python wrapper code
TARGET = $(addprefix $(OBJDIR)/, $(LIBNAME))

# objects to build
OBJ = $(addprefix $(OBJDIR)/, conncomponents.o opfunctions.o \
		  pyfunctions.o qlmfunctions.o pyutil.o neighbours.o)

$(TARGET).so: $(TARGET).o $(OBJ)
	g++ -shared -Wl,-soname,"$(LIBNAME).so" -L$(BOOST_LIB) $(OBJ) $(TARGET).o -lboost_python -fPIC -o $(LIBNAME).so

conncomponents.o: conncomponents.h typedefs.h particle.h box.h

opfunctions.o: opfunctions.cpp constants.h

pyfunctions.o: pyfunctions.cpp qlmfunctions.h box.h constants.h particle.h\
               pyutil.h

qlmfunctions.o: qlmfunctions.cpp particle.h box.h constants.h opfunctions.h\
                conncomponents.h utility.h	

neighbours.o: neighbours.cpp particle.h box.h constants.h

pyutil.o: pyutil.cpp particle.h

clean:
	rm $(LIBNAME).so $(OBJDIR)/*.o
