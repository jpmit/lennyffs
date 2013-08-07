// utility.h
// James Mithen
// j.mithen@surrey.ac.uk

// Function definitions for some handy 'utility' functions.

#ifndef UTILITY_H
#define UTILITY_H

#include<vector>

inline std::vector<int> range(int start, int end)
{
	  std::vector<int> ret;
	  ret.reserve(end - start);
	  for (int i = start; i != end; ++i)
			 ret.push_back(i);
	  return ret;
}

inline void reindex(std::vector<int>& indx1, const std::vector<int>& indx2)
{
	  for (std::vector<int>::size_type i = 0; i != indx1.size(); ++i) {
			 indx1[i] = indx2[indx1[i]];
	  }
}

#endif
