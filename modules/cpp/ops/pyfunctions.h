// pyfunctions.h
// James Mithen
// j.mithen@surrey.ac.uk

// Declarations for functions to be called by Python via Boost.Python

#include <boost/python.hpp>
#include "boost/python/extract.hpp"
#include "boost/python/numeric.hpp"

double py_fracsolidld(boost::python::numeric::array,
							 boost::python::numeric::array,
							 boost::python::numeric::array, 
						    const int, const int, const double, const double,
							 const double, const bool, const double);
double py_fracsolidtf(boost::python::numeric::array,
							 boost::python::numeric::array,
							 boost::python::numeric::array, 
						    const int, const int, const double, const double,
							 const double, const bool, const double, const int,
							 const double);
double py_nclustf(boost::python::numeric::array,
						boost::python::numeric::array,
						boost::python::numeric::array, 
						const int, const int,
						const double, const double,
						const double, const bool, const double,
						const int, const double);
double py_nclusld(boost::python::numeric::array,
						boost::python::numeric::array,
						boost::python::numeric::array, 
						const int, const int, const double, const double,
						const double, const bool, const double);
