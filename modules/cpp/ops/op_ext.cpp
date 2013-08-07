#include <boost/python.hpp>
#include "boost/python/numeric.hpp"
#include "pyfunctions.h"
 
BOOST_PYTHON_MODULE(op_ext)
{
	  boost::python::numeric::array::set_module_and_type("numpy", "ndarray");

	  def("nclustf", py_nclustf);
	  def("nclusld", py_nclusld);	  
	  def("fracsolidld", py_fracsolidld);
	  def("fracsolidtf", py_fracsolidtf);	  

}
