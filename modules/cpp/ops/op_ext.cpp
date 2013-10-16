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
	  //enum_<LDCLASS>("LDCLASS")
	  //	 .value("FCC", FCC)
	  //		 .value("HCP", HCP)
	  //	 .value("BCC", BCC)
	  //	 .value("LIQUID", LIQUID)
	  //	 .value("ICOS", ICOS)
	  //	 .value("SURFACE", SURFACE);
	  to_python_converter<LDCLASS, enum_to_python_int<LDCLASS> >();
	  to_python_converter<std::vector<LDCLASS>,
								 vector_to_python_list<LDCLASS> >();
	  to_python_converter<std::vector<int>,
								 vector_to_python_list<int> >();

	  def("nclustf", py_nclustf);
	  def("nclusld", py_nclusld);	  
	  def("fracsolidld", py_fracsolidld);
	  def("fracsolidtf", py_fracsolidtf);
	  def("ldclass", py_ldclass);
	  def("largestcluster", py_largestcluster);
}
