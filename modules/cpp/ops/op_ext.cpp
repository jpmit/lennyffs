// op_ext.cpp
// James Mithen
// j.mithen@surrey.ac.uk

// Code for generating order parameter extension module using
// Boost.Python.

#include <boost/python.hpp>
#include "boost/python/numeric.hpp"
#include "pyfunctions.h"

using namespace boost::python;
 
BOOST_PYTHON_MODULE(op_ext)
{
   numeric::array::set_module_and_type("numpy", "ndarray");

   // Registration of converters.  These are defined in pyfunctions.h
   to_python_converter<TFCLASS, enum_to_python_int<TFCLASS> >();
   to_python_converter<LDCLASS, enum_to_python_int<LDCLASS> >();
   to_python_converter<std::vector<TFCLASS>, vector_to_python_list<TFCLASS> >();
   to_python_converter<std::vector<LDCLASS>, vector_to_python_list<LDCLASS> >();
   to_python_converter<std::vector<int>, vector_to_python_list<int> >();
   to_python_converter<std::vector<double>, vector_to_python_list<double> >();

   // The string is the function name called from python
   def("q6global", py_q6global);   
   def("nclustf", py_nclustf);
   def("tfclass", py_tfclass);
   def("nclusld", py_nclusld);
   def("ncluspolyld", py_ncluspolyld);   
   def("fracsolidld", py_fracsolidld);
   def("fracsolidtf", py_fracsolidtf);
   def("ldclass", py_ldclass);
   def("largestcluster", py_largestcluster);
   def("q4w4q6w6", py_q4w4q6w6);
   def("numneighcut", py_numneighcut);
}
