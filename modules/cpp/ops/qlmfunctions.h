// qlmfunctions.h
// James Mithen
// j.mithen@surrey.ac.uk

#ifndef QLMFUNCTIONS_H
#define QLMFUNCTIONS_H

#include <vector>
#include <complex>
#include "particle.h"
#include "box.h"
#include "typedefs.h"
#include "constants.h"

// Compute the spherical harmonics etc.
array2d qlms(const std::vector<Particle>&, const Box&,
				 const std::vector<int>&,
				 const std::vector<std::vector<int> >&, const int);

// Get normalised qlm matrix 
array2d qlmtildes(const array2d&, const std::vector<int>&,
						const int);

// Computes the LD 'averaged' qlms from the qlm matrix
array2d qlmbars(const array2d&, const std::vector<std::vector<int> >&,
					 const int);

// Q value of some particles
double Qpars(const array2d&, const std::vector<int>&, const int);

// vector of ql(i) values for every particle
std::vector<double> qls(const array2d&);

// vector of wl(i) values for every particle
std::vector<double> wls(const array2d&);

// Number of crystalline links that each particle has (tenWolde-Frenkel)
std::vector<int> getnlinks(const array2d&, const std::vector<int>&,
									const std::vector<std::vector<int> >&, const int,
									const int, const double, const int);

// Classify particles as BCC, HCP, FCC, etc. according to LD method
std::vector<LDCLASS> classifyparticlesld(const int,
													  const std::vector<double>&,
													  const std::vector<double>&,
													  const std::vector<double>&,
													  const std::vector<double>&);

// Classify particles as Liquid or Solid according to TF method
std::vector<TFCLASS> classifyparticlestf(const std::vector<int>&, const int,
													  const int);

// Largest cluster of XTAL particles, and fraction of XTAL particles in system
std::vector<int> largestclusterld(const std::vector<Particle>&, const Box&,
											 const std::vector<LDCLASS>&);
std::vector<int> largestclustertf(const std::vector<Particle>&, const Box&,
											 const std::vector<TFCLASS>&);
double fracsolidtf(const std::vector<TFCLASS>&, const int);
double fracsolidld(const std::vector<LDCLASS>&, const int);

#endif
