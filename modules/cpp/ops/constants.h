// constants.h
// James Mithen
// j.mithen@surrey.ac.uk

// Enums for classifying particles, numbers for computing Wigner
// symbols, and a few other handy constants.

#ifndef CONSTANTS_H
#define CONSTANTS_H

const double PI = 3.14159265358979;

// Classification constants for Lechner Dellago (LD)
// and Ten-Wolde Frenkel (TF) criterion.  For LD, we classify
// each particle as FCC (0), HCP (1), BCC (2), LIQUID (3), ICOS (4),
// OTHER (5).  For TF, we classify each particle as LIQUID (0) or
// XTAL (1). Unfortunately we need different names for LD and TF to
// avoid collision here (this is a language constraint for enums).

enum LDCLASS {FCC, HCP, BCC, LIQUID, ICOS, SURFACE};

enum TFCLASS {LIQ, XTAL, SURF};

// Wigner symbols (l  l  l )
//                (m1 m2 m3) 
// for l = 4 and 6, and integers m1 m2 m3
//
// We have |m1| <= l, |m2| <= l, |m3| <= l,
// and symbols are only non-zero when m1 + m2 + m3 = 0. 
// Note that from properties of Wigner symbols
// (l  l  l )  == (l    l   l )
// (m1 m2 m3)     (-m1 -m2 -m3)  (when l = 4 or 6)
// The values are hard coded below, along with the number of
//  distinct permutations (including m1 m2 m3 and -m1 -m2 -m3).

// for l = 6 there are 16 different values
const double WIGNER6[] = { 0.0511827,  //6 -6  0 (6)
								  -0.0957541,  //6 -5 -1 (12)
								   0.129115,   //6 -4 -2 (12)
								  -0.141438,   //6 -3 -3 (6)
								   0.127957,   //5 -5  0 (6)
								  -0.106079,   //5 -4 -1 (12)
 								   0.0408297,  //5 -3 -2 (12)
								   0.0186119,  //4 -4  0 (6)
								   0.0688184,  //4 -3 -1 (12)
								  -0.104459,   //4 -2 -2 (6)
								  -0.100039,   //3 -3  0 (6)
								   0.0452321,  //3 -2 -1 (12)
								   0.0511827,  //2 -2  0 (6)
								  -0.0953576,  //2 -1 -1 (6)
								   0.0465298,  //1 -1  0 (6)
								  -0.0930595   //0  0  0 (1)
};
const int WIGNERINDX6[16][3] = { {6, -6, 0},
											{6, -5, -1},
											{6, -4, -2},
											{6, -3, -3},
											{5, -5, 0},
											{5, -4, -1},
											{5, -3, -2},
											{4, -4, 0},
											{4, -3, -1},
											{4, -2, -2},
											{3, -3, 0},
											{3, -2, -1},
											{2, -2, 0},
											{2, -1, -1},
											{1, -1, 0},
											{0, 0, 0}
};
const int WIGNERPERM6[] = {6,12,12,6,6,12,12,6,12,6,6,12,6,6,6,1};
	  
// for l = 4 there are 9 different values
const double WIGNER4[] = { 0.104298,   // 4 -4  0 (6)
								  -0.164909,   // 4 -3 -1 (12)
								   0.186989,   // 4 -2 -2 (6)
								   0.156447,   // 3 -3  0 (6)
								  -0.0623298,  // 3 -2 -1 (12)
								  -0.0819482,  // 2 -2  0 (6)
								   0.141351,   // 2 -1 -1 (6)
								  -0.0670485,  // 1 -1  0 (6)
								   0.134097,   // 0  0  0 (1)
};
const int WIGNERINDX4[9][3] = { {4, -4,  0},
										  {4, -3, -1},
										  {4, -2, -2},
										  {3, -3,  0},
										  {3, -2, -1},
										  {2, -2,  0},
										  {2, -1, -1},
										  {1, -1,  0},
										  {0,  0,  0},
};
const int WIGNERPERM4[] = {6,12,6,6,12,6,6,6,1};

#endif
