// pyutil.h
// James Mithen
// j.mithen@surrey.ac.uk

#ifndef PYUTIL_H
#define PYUTIL_H

#include "particle.h"
#include <vector>
#include "boost/python/numeric.hpp"

std::vector<Particle> getparticles(boost::python::numeric::array&,
                                   boost::python::numeric::array&,
                                   boost::python::numeric::array&,
                                   const int);

#endif
