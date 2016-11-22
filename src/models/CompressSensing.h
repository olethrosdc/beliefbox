#ifndef COMPRESS_SENSING_H
#define COMPRESS_SENSING_H

#include <vector>
#include <map>

#include "real.h"

struct compressor {
	std::vector<std::vector<real> > compressionMatrix;
	std::map<int,int> indexMap;
};
compressor createCompressionMatrix(std::vector<std::vector<int> > &indices);
std::vector<real> compressSparseArray(std::vector<int> &values, std::vector<int> &indices, compressor c);

#endif