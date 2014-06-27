// box.h
// James Mithen
// j.mithen@surrey.ac.uk

// A simple class for a (cuboid shaped) simulation box. This makes
// applying the periodic bounary conditions a little more elegant.

#ifndef BOX_H
#define BOX_H

#include <cmath>
#include <vector>
#include "particle.h"
#include "box.h"

class Box
{
public:
   Box(double lx, double ly, double lz, double ns = 1.5, bool pz = false) :
   lboxx(lx), lboxy(ly), lboxz(lz), nsep(ns), nsepsq(ns * ns), periodicz(pz){}
   inline void sep(const Particle& p1, const Particle& p2, double* s) const;
   inline double sepsq(const Particle& p1, const Particle& p2) const;
   inline bool isneigh(const Particle& p1, const Particle& p2, double& r2) const;
   inline bool isneigh(double* s, double&r2) const;

private:
   double lboxx;
   double lboxy;
   double lboxz;       
   double nsep;
   double nsepsq;
   bool periodicz;
};

// Square of separation between two particles.

inline double Box::sepsq(const Particle& p1, const Particle& p2) const
{
   double s[3];
   sep(p1 ,p2, s);
   return (s[0] * s[0] + s[1] * s[1] + s[2] * s[2]);
}

// Separation between two particles modulo periodic bcs.

inline void Box::sep(const Particle& p1, const Particle& p2, double* s) const
{
   double sepx,sepy,sepz;
   sepx = p1.pos[0] - p2.pos[0];
   sepy = p1.pos[1] - p2.pos[1];
   sepz = p1.pos[2] - p2.pos[2];
     
   if (sepx > 0.5 * lboxx) {
      sepx = sepx - lboxx;
   }
   else if (sepx < -0.5 * lboxx) {
      sepx = sepx + lboxx;
   }
   if (sepy > 0.5 * lboxy) {
      sepy = sepy - lboxy;
   }
   else if (sepy < -0.5 * lboxy) {
      sepy = sepy + lboxy;
   }
   if (periodicz) {              
      if (sepz > 0.5 * lboxz) {
         sepz = sepz - lboxz;
      }
      else if (sepz < -0.5 * lboxz) {
         sepz = sepz + lboxz;
      }
   }

   s[0] = sepx;
   s[1] = sepy;
   s[2] = sepz;
   return;
}

// Are p1 and p2 neighbours?  Also return square separation modulo periodic bcs.

inline bool Box::isneigh(const Particle& p1, const Particle& p2, double& rsq) const
{
   double s[3];

   // compute separation between particles (store in s)
   sep(p1, p2, s);

   if ((std::abs(s[0]) < nsep) && (std::abs(s[1]) < nsep)
       && (std::abs(s[2]) < nsep)) {
      rsq = s[0] * s[0] + s[1] * s[1] + s[2] * s[2];
      if (rsq < nsepsq) {
         return true;
      }
   }
   return false;
}

// Same as above but first argument points to separation array.

inline bool Box::isneigh(double *s, double& rsq) const
{
   if ((std::abs(s[0]) < nsep) && (std::abs(s[1]) < nsep)
       && (std::abs(s[2]) < nsep)) {
      rsq = s[0] * s[0] + s[1] * s[1] + s[2] * s[2];
      if (rsq < nsepsq) {
         return true;
      }
   }
   return false;
}

#endif
