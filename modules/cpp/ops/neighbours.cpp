// neighbours.cpp
// James Mithen
// j.mithen@surrey.ac.uk

// Functions to compute neighbour separations.

#include <iostream>
#include <vector>
#include "particle.h"
#include "box.h"
#include "float.h"

using std::vector;

// get neighbour list and num neighbours where neighbours are defined
// as being the closest 'n' particles to a given particle (usually we
// will use n = 12).

void neighnearest(const vector<Particle>& allpars, const Box& simbox,
                  vector<int>& numneigh, vector<vector<int> >& lneigh,
                  const int n)
{
   vector<Particle>::size_type i, j;
   vector<int>::size_type k, m, indx;
   double sepsq;
   vector<vector<int> > bestsep = lneigh;

   const vector<Particle>::size_type npar = allpars.size();

   // fill up numneigh, since by construction here every particle
   // has 'n' neighbours, this is easy.
   for (k = 0; k != numneigh.size(); ++k) {
      numneigh[k] = n;
   }

   // fill up each neighbour list with token large value and with
   // neighbour set to 0 (this could be any number)
   for (k = 0; k != lneigh.size(); ++k) {
      for (m = 0; m != n; ++m) {
         lneigh[k].push_back(0);
         bestsep[k].push_back(DBL_MAX - 1);
      }
   }

   // go through each particle in turn, and find the 'n' nearest
   // neighbours.
   for (i = 0; i != npar; ++i) {
      for (j = 0; j != npar; ++j) {
         if (i != j) {
            sepsq = simbox.sepsq(allpars[i], allpars[j]);

            // if sep is less than largest current neighbour
            // separation, then this particle is a neighbour
            if (sepsq < bestsep[i][n - 1]) {

               // find the correct index
               indx = n - 1;
               while (sepsq < bestsep[i][indx] and indx != -1) {
                  --indx;
               }
               // we want to insert the new neighbour at
               // position indx + 1, and shift everything to
               // the right of it along one place, knocking
               // off the end one..

               lneigh[i].insert(lneigh[i].begin() + indx + 1, j);
               bestsep[i].insert(bestsep[i].begin() + indx + 1, sepsq);

               // delete last element
               lneigh[i].pop_back();
               bestsep[i].pop_back();
            }
         }
      }
   }
}

// get neighbour list and num neighbours where neighbours are defined
// as being all particles within some cutoff radius 'rcut' of a given
// particle.

void neighcut(const vector<Particle>& allpars, const Box& simbox,
              vector<int>& numneigh, vector<vector<int> >& lneigh)
{
   const vector<Particle>::size_type npar = allpars.size();

   double r2, sep[3];
   vector<Particle>::size_type i,j;
          
   for (i = 0; i != npar; ++i) {
      for (j = 0; j != npar; ++j) {
         if (i != j) {
            if (simbox.isneigh(allpars[i], allpars[j], r2)) {
               // particles i and j are neighbours
               lneigh[i].push_back(j);
            }
         }
      }
      // set number of neighbours of particle i
      numneigh[i] = lneigh[i].size();
   }
}
