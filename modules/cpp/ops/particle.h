// particle.h
// James Mithen
// j.mithen@surrey.ac.uk

// A simple particle class for particle simulations e.g. MC/MD

#ifndef PARTICLE_H
#define PARTICLE_H

struct Particle
{
	  double pos[3];
	  char symbol; // for outputting e.g. jmol
};

#endif
