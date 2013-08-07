// conncomponents.h
// James Mithen
// j.mithen@surrey.ac.uk

#ifndef CONNCOMPONENTS_H
#define CONNCOMPONENTS_H

#include "typedefs.h"
#include "particle.h"
#include "box.h"

graph getxgraph(const std::vector<Particle>&, const std::vector<int>&,
					 const Box&);
std::vector<int> largestcomponent(const graph&);

#endif
