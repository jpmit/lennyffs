// pyutil.cpp
// James Mithen
// j.mithen@surrey.ac.uk

// Utility functions needed to make the Boost.Python interface work
// properly.

#include <vector>
#include <boost/python.hpp>
#include "boost/python/extract.hpp"
#include "boost/python/numeric.hpp"
#include "particle.h"

using std::vector;

// get vector of "Particle" type from numeric arrays of particle
// positions.

vector<Particle> getparticles(boost::python::numeric::array& xpos,
                              boost::python::numeric::array& ypos,
                              boost::python::numeric::array& zpos,
                              const int npartot)
{

   vector<Particle> allpars(npartot);
          
   // we will probably want to remove the "Particle" class from the
   // code at some point to prevent the conversion below.  Also, I
   // have no idea how expensive the 'extract' functionality
   // provided by boost is...
   for (vector<Particle>::size_type i = 0; i != allpars.size(); ++i) {
      allpars[i].pos[0] = boost::python::extract<double>(xpos[i]);
      allpars[i].pos[1] = boost::python::extract<double>(ypos[i]);
      allpars[i].pos[2] = boost::python::extract<double>(zpos[i]);
   }

   return allpars;
}
