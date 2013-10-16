// pyfunctions.cpp
// James Mithen
// j.mithen@surrey.ac.uk

// Functions that can be called by Python via the Boost.Python
// interface should be defined here.  Functions to be called
// by Python are prefixed by py_ .
// Currently implemented:
// py_nclustf     - size of largest crystalline cluster, according to
//                  ten-Wolde Frenkel method.
// py_fracsolidtf - fraction of crystalline particles (excluding surface
//                  particles) according to ten-Wolde Frenkel method.
// py_nclusld     - size of largest crystalline cluster, according to
//                  Lechner Dellago method.
// py_fracsolidld - fraction of crystalline particles (excluding surface
//                  particles) according to Lecher Dellago method.

#include <iostream>
#include <vector>
#include <boost/python.hpp>
#include "boost/python/numeric.hpp"
#include "qlmfunctions.h"
#include "box.h"
#include "constants.h"
#include "particle.h"
#include "pyutil.h"
#include "neighbours.h"
#include "typedefs.h"
#include "conncomponents.h"
#include "utility.h"

using std::vector;
using std::cout;

// size of largest crystalline cluster, according to TF method.

double py_nclustf(boost::python::numeric::array xpos,
						boost::python::numeric::array ypos,
						boost::python::numeric::array zpos,
						const int npartot, const int nparsurf,
						const double lboxx, const double lboxy,
						const double lboxz, const bool zperiodic,
						const double nsep, const int nlinks,
						const double linkval)
{
	  // create vector of type "Particle"
	  vector<Particle> allpars = getparticles(xpos, ypos, zpos,
															npartot);

	  // create "Box"
	  Box simbox(lboxx, lboxy, lboxz, nsep, zperiodic);

	  // store number of neighbours and neighbour list
  	  vector<int> numneigh(npartot, 0); // num neighbours for each particle
	  vector<vector<int> > lneigh(npartot); // vector of neighbour particle nums for
	                                        // each par

	  // fill up numneigh and lneigh, we can use either neighcut of
	  // neighnearest for this
	  neighnearest(allpars, simbox, numneigh, lneigh, 12);

	  // matrix of qlm values, for l = 6 only
	  array2d q6lm(boost::extents[npartot][13]);
	  q6lm = qlms(allpars, simbox, numneigh, lneigh, 6);

	  // compute number of crystalline 'links'
	  // first get normalised vectors qlm (-l <= m <= l) for computing
	  // dot product Sij
	  array2d qlmt = qlmtildes(q6lm, numneigh, 6);

	  // do dot products Sij to get number of links
	  vector<int> numlinks = getnlinks(qlmt, numneigh, lneigh, nparsurf,
												  nlinks, linkval, 6);

	  // classify particles using q4lbar etc.
	  vector<TFCLASS> tfclass = classifyparticlestf(numlinks, nlinks, nparsurf);

	  // indices of particles in the largest cluster
	  vector<int> tfcnums = largestclustertf(allpars, simbox, tfclass);
	  return tfcnums.size();


}

// fraction of solid particles (excluding surface particles) according
// to TF method.

double py_fracsolidtf(boost::python::numeric::array xpos,
							 boost::python::numeric::array ypos,
							 boost::python::numeric::array zpos,
						    const int npartot,  const int nparsurf,
							 const double lboxx, const double lboxy,
							 const double lboxz, const bool zperiodic,
							 const double nsep, const int nlinks,
	                   const double linval)
{
	  // create vector of type "Particle"
	  vector<Particle> allpars = getparticles(xpos, ypos, zpos,
															npartot);


	  // create "Box"
	  Box simbox(lboxx, lboxy, lboxz, nsep, zperiodic);

	  // store number of neighbours and neighbour list
  	  vector<int> numneigh(npartot, 0); // num neighbours for each particle
	  vector<vector<int> > lneigh(npartot); // vector of neighbour particle nums for
	                                        // each par

	  // fill up numneigh and lneigh, we can use either neighcut of
	  // neighnearest for this
	  neighnearest(allpars, simbox, numneigh, lneigh, 12);

	  // matrix of qlm values, for l = 6 only
	  array2d q6lm(boost::extents[npartot][13]);
	  q6lm = qlms(allpars, simbox, numneigh, lneigh, 6);

	  // compute number of crystalline 'links'
	  // first get normalised vectors qlm (-l <= m <= l) for computing
	  // dot product Sij
	  array2d qlmt = qlmtildes(q6lm, numneigh, 6);

	  // do dot products Sij to get number of links
	  vector<int> numlinks = getnlinks(qlmt, numneigh, lneigh, nparsurf,
												  nlinks, linval, 6);

	  // classify particles using q4lbar etc.
	  vector<TFCLASS> tfclass = classifyparticlestf(numlinks, nlinks,
																	nparsurf);

	  return fracsolidtf(tfclass, nparsurf);
}

// size of largest crystalline cluster, according to LD method.

double py_nclusld(boost::python::numeric::array xpos,
						boost::python::numeric::array ypos,
						boost::python::numeric::array zpos,
						const int npartot, const int nparsurf,
						const double lboxx, const double lboxy,
						const double lboxz, const bool zperiodic,
						const double nsep)
{
	  // create vector of type "Particle"
	  vector<Particle> allpars = getparticles(xpos, ypos, zpos,
															npartot);


	  // create "Box"
	  Box simbox(lboxx, lboxy, lboxz, nsep, zperiodic);

	  // store number of neighbours and neighbour list
  	  vector<int> numneigh(npartot, 0); // num neighbours for each particle
	  vector<vector<int> > lneigh(npartot); // vector of neighbour particle nums for
	                                         // each par

	  // fill up numneigh and lneigh, we can use either neighcut of
	  // neighnearest for this
	  neighnearest(allpars, simbox, numneigh, lneigh, 12);

	  // matrix of qlm values, for l = 4 and l = 6
	  array2d q4lm(boost::extents[npartot][9]);
	  array2d q6lm(boost::extents[npartot][13]);
	  q4lm = qlms(allpars, simbox, numneigh, lneigh, 4);	  
	  q6lm = qlms(allpars, simbox, numneigh, lneigh, 6);
	  
	  // Lechner dellago eq 6, for l = 4 and l = 6
	  array2d q4lmb = qlmbars(q4lm, lneigh, 4);
	  array2d q6lmb = qlmbars(q6lm, lneigh, 6);

	  // lechner dellago eq 5, for l = 4 and l = 6
	  vector<double> q4lbar = qls(q4lmb);
	  vector<double> w4lbar = wls(q4lmb);
	  vector<double> q6lbar = qls(q6lmb);
	  vector<double> w6lbar = wls(q6lmb);

	  // classify particles using q4lbar etc.
	  vector<LDCLASS> ldclass = classifyparticlesld(nparsurf, q4lbar,
																	q6lbar, w4lbar, w6lbar);

	  // indices of particles in the largest cluster
	  vector<int> ldcnums = largestclusterld(allpars, simbox, ldclass);
	  return ldcnums.size();
}

// fraction of solid particles (excluding surface particles) according
// to LD method.

double py_fracsolidld(boost::python::numeric::array xpos,
							 boost::python::numeric::array ypos,
							 boost::python::numeric::array zpos,
						    const int npartot,  const int nparsurf,
							 const double lboxx, const double lboxy,
							 const double lboxz, const bool zperiodic,
							 const double nsep)
{
	  // create vector of type "Particle"
	  vector<Particle> allpars = getparticles(xpos, ypos, zpos,
															npartot);


	  // create "Box"
	  Box simbox(lboxx, lboxy, lboxz, nsep, zperiodic);

	  // store number of neighbours and neighbour list
  	  vector<int> numneigh(npartot, 0); // num neighbours for each particle
	  vector<vector<int> > lneigh(npartot); // vector of neighbour particle nums for
	                                         // each par

	  // fill up numneigh and lneigh, we can use either neighcut of
	  // neighnearest for this
	  neighnearest(allpars, simbox, numneigh, lneigh, 12);

	  // matrix of qlm values, for l = 4 and l = 6
	  array2d q4lm(boost::extents[npartot][9]);
	  array2d q6lm(boost::extents[npartot][13]);
	  q4lm = qlms(allpars, simbox, numneigh, lneigh, 4);	  
	  q6lm = qlms(allpars, simbox, numneigh, lneigh, 6);
	  
	  // Lechner dellago eq 6, for l = 4 and l = 6
	  array2d q4lmb = qlmbars(q4lm, lneigh, 4);
	  array2d q6lmb = qlmbars(q6lm, lneigh, 6);

	  // lechner dellago eq 5, for l = 4 and l = 6
	  vector<double> q4lbar = qls(q4lmb);
	  vector<double> w4lbar = wls(q4lmb);
	  vector<double> q6lbar = qls(q6lmb);
	  vector<double> w6lbar = wls(q6lmb);

	  // classify particles using q4lbar etc.
	  vector<LDCLASS> ldclass = classifyparticlesld(nparsurf, q4lbar,
																	q6lbar, w4lbar, w6lbar);
	  
	  return fracsolidld(ldclass, nparsurf);
}

// classification particles using LD

vector<LDCLASS> py_ldclass(boost::python::numeric::array xpos,
									boost::python::numeric::array ypos,
									boost::python::numeric::array zpos,
									const int npartot, const int nparsurf,
									const double lboxx, const double lboxy,
									const double lboxz, const bool zperiodic,
									const double nsep)
{
	  // create vector of type "Particle"
	  vector<Particle> allpars = getparticles(xpos, ypos, zpos,
															npartot);


	  // create "Box"
	  Box simbox(lboxx, lboxy, lboxz, nsep, zperiodic);

	  // store number of neighbours and neighbour list
  	  vector<int> numneigh(npartot, 0); // num neighbours for each particle
	  vector<vector<int> > lneigh(npartot); // vector of neighbour
														 // particle nums for each
														 // par

	  // fill up numneigh and lneigh, we can use either neighcut of
	  // neighnearest for this
	  //neighcut(allpars, simbox, numneigh, lneigh);	  
	  neighnearest(allpars, simbox, numneigh, lneigh, 12);

	  // matrix of qlm values, for l = 4 and l = 6
	  array2d q4lm(boost::extents[npartot][9]);
	  array2d q6lm(boost::extents[npartot][13]);
	  q4lm = qlms(allpars, simbox, numneigh, lneigh, 4);	  
	  q6lm = qlms(allpars, simbox, numneigh, lneigh, 6);
	  
	  // Lechner dellago eq 6, for l = 4 and l = 6
	  array2d q4lmb = qlmbars(q4lm, lneigh, 4);
	  array2d q6lmb = qlmbars(q6lm, lneigh, 6);

	  // lechner dellago eq 5, for l = 4 and l = 6
	  vector<double> q4lbar = qls(q4lmb);
	  vector<double> w4lbar = wls(q4lmb);
	  vector<double> q6lbar = qls(q6lmb);
	  vector<double> w6lbar = wls(q6lmb);

	  // classify particles using q4lbar etc.
	  vector<LDCLASS> ldclass = classifyparticlesld(nparsurf, q4lbar,
																	q6lbar, w4lbar, w6lbar);

	  return ldclass;
}

// indices of particles in largest cluster
vector<int> py_largestcluster(boost::python::numeric::array cposx,
										boost::python::numeric::array cposy,
										boost::python::numeric::array cposz,
										const int npar, const double lboxx,
										const double lboxy, const double lboxz,
										const bool zperiodic, const double nsep)
{
	  // create vector of type "Particle"
	  vector<Particle> cpars = getparticles(cposx, cposy, cposz,
														 npar);

	  // create "Box"
	  Box simbox(lboxx, lboxy, lboxz, nsep, zperiodic);

	  // graph of xtal particles, with each particle a vertex and each
	  // link an edge
	  graph xgraph = getxgraph(cpars, range(0, npar), simbox);

	  // largest cluster is the largest connected component of graph
	  vector<int> cnums = largestcomponent(xgraph);

	  return cnums;
}
