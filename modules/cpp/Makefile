# location of the Python header files

PYTHON_VERSION = 2.7
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
TARGET = $(addprefix $(OBJDIR)/, op_ext)
OBJ = $(addprefix $(OBJDIR)/, conncomponents.o opfunctions.o pyfunctions.o qlmfunctions.o pyutil.o)

$(TARGET).so: $(TARGET).o $(OBJ)
	g++ -shared -Wl,--export-dynamic $(TARGET).o $(OBJ) -L$(BOOST_LIB) -l:libboost_python.so -L/usr/lib/python$(PYTHON_VERSION)/config -lpython$(PYTHON_VERSION) -o $(LIBNAME).so

#$(TARGET).o: $(TARGET).c
#	g++ -I$(PYTHON_INCLUDE) -I$(BOOST_INC) -fPIC -c $(TARGET).c

$(TARGET): $(TARGET).cpp pyfunctions.h

conncomponents.o: conncomponents.h typedefs.h particle.h box.h

opfunctions.o: opfunctions.cpp constants.h

pyfunctions.o: pyfunctions.cpp qlmfunctions.h box.h constants.h particle.h\
               pyutil.h

qlmfunctions.o: qlmfunctions.cpp particle.h box.h constants.h opfunctions.h\
                conncomponents.h utility.h	

pyutil.o: pyutil.cpp particle.h

#$(OBJ) : $(SRC)
#	g++ -fPIC -c $(SRC)

clean:
	rm $(TARGET).so $(OBJDIR)/*.o
