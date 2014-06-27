// opfunctions.cpp
// James Mithen
// j.mithen@surrey.ac.uk
//
// Mathematical functions needed to compute the 'Steindhardt' order
// parameters.  See e.g. Steinhardt, Nelson, Ronchetti, Phys. Rev. B
// 28, 784 (1983).

#include <cmath>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <complex>
#include "constants.h"

using std::abs;
using std::pow;
using std::vector;
using std::complex;

// Factorial of an integer.

inline int fact(const int n)
{
   int ifact = 1;
   if (n < 2) {
      return ifact;
   }
   for (int j = 2; j != n + 1; ++j) {
      ifact *= j;
   }
   return ifact;
}

// Computes the associated Legendre polynomial P_l^m(x), for positive
// m only.

double plm(const int lval, const int mval, const double x)
{
   double somx2,pm1m,pll;

   // first compute P_m^m
   double pmm = 1.0;
   double fact = 1.0;

   if (mval > 0) {
      somx2 = sqrt(1.0 - x*x);
      for (int i = 1; i <= mval; ++i) {
         pmm = -pmm*fact*somx2;
         fact = fact + 2.0;
      }
   }

   if (lval == mval) {
      return pmm;
   }
   else {
      // compute P_m+1^m
      pm1m = (2*mval + 1)*x*pmm;
      if (lval == mval + 1) {
         return pm1m;
      }
      else {
         // compute P_l^m if not already returned
         for (int i = mval + 2; i <= lval; ++i) {
            pll = (x*(2*i - 1)*pm1m - (i+mval-1)*pmm)/(i-mval);
            pmm = pm1m;
            pm1m = pll;
         }
         return pll;
      }
   }
}

// Spherical harmonic Y(l,m,theta,phi).

complex<double> ylm(const int lval, const int mval,
                    const double costheta, const double phi)
{
   int absm = abs(mval);
   double coeff = sqrt((2.0 * lval + 1.0)*fact(lval - absm) /
                       (4.0 * PI * fact(lval + absm)));
   if (mval < 0) {
      coeff = coeff * pow(-1, mval);
   }

   double plmval = plm(lval, absm, costheta);

   complex<double> res(cos(mval * phi), sin(mval * phi));
   res = res * coeff * plmval;

   return res;
}
