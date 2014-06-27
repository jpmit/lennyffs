// typedefs.h
// James Mithen
// j.mithen@surrey.ac.uk

// typedefs for use in the rest of the code.  array2d is a 2d array of
// complex doubles, used for computing the spherical harmonics, and
// qlm values; tensor is a 2d array of real numbers (double precision),
// used for the gyration tensor.  graph is a boost adjacency list, used
// for finding the number of connected components (and hence the size of
// the largest crystalline cluster).

#ifndef TYPEDEFS_H
#define TYPEDEFS_H

#include <boost/graph/adjacency_list.hpp>
#include <boost/multi_array.hpp>

typedef boost::multi_array<std::complex<double>,2> array2d;
typedef boost::multi_array<double,2> tensor;
typedef boost:: adjacency_list <boost::vecS, boost::vecS, boost::undirectedS> graph;

#endif
