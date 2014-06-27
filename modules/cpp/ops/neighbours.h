// neighbours.h
// James Mithen
// j.mithen@surrey.ac.uk

#ifndef NEIGHBOURS_H
#define NEIGHBOURS_H

#include <vector>
#include "particle.h"
#include "box.h"

void neighnearest(const std::vector<Particle>&, const Box&,
                  std::vector<int>&, std::vector<std::vector<int> >&,
                  const int);
void neighcut(const std::vector<Particle>&, const Box&,
              std::vector<int>&, std::vector<std::vector<int> >&);

#endif

