// typedefs.h
// James Mithen
// j.mithen@surrey.ac.uk

#ifndef TYPEDEFS_H
#define TYPEDEFS_H

#include <boost/graph/adjacency_list.hpp>
#include <boost/multi_array.hpp>

typedef boost::multi_array<std::complex<double>,2> array2d;
typedef boost::multi_array<double,2> tensor;
typedef boost:: adjacency_list <boost::vecS, boost::vecS,
										  boost::undirectedS> graph;

#endif
