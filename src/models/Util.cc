#include "Util.h"
#include <algorithm>

std::vector<real> findNMaximum(std::vector<real> vals, int count) {
	assert(count<vals.size());
	std::sort(vals.begin(),vals.end());
	std::vector<real> newVec = std::vector<real>(vals.end()-count,vals.end());
	std::reverse(newVec.begin(),newVec.end());
	return newVec;
}