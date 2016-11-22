#include <vector>
#include <map>
#include <random>
#include "CompressSensing.h"
#include "real.h"

compressor createCompressionMatrix(std::vector<std::vector<int> > &indices) {
	compressor c;
	std::map<int,int> indexMap;

	int distinct = 0;
	for(int i = 0;i<indices.size();i++) {
		for(int j=0;j<indices[i].size();j++) {
			if(indexMap.find(indices[i][j])==indexMap.end()) {
				indexMap.insert(std::map<int,int>::value_type(indices[i][j], distinct));
				distinct++;
			}
		}
	}
	int tmpSize = 0;
	// find argmax_i S-sparse 
	for(int i = 0;i<indices.size();i++) {
		if(indices[i].size() > tmpSize) {
			tmpSize = indices[i].size();
		}
	}
	
	// to have the RIP-property K = O(S * log(N/S) for S-sparse
	double ratio = (double) distinct/ (double) tmpSize;
	int resultSize = (int)(tmpSize * std::log(ratio)+0.5);

	std::vector<std::vector<real> > compressionMatrix(resultSize);
	std::normal_distribution<double> distribution(0.0,1.0);
	std::default_random_engine gen;
	for(int i = 0;i<resultSize;i++) {
		compressionMatrix[i].resize(distinct);
		for(int j = 0;j<distinct;j++) {
			compressionMatrix[i][j] = distribution(gen);
		}
	}
	c.compressionMatrix = compressionMatrix;
	c.indexMap = indexMap;
	return c;
}
std::vector<real> compressSparseArray(std::vector<int> &values, std::vector<int> &indices, compressor c) {
	int resultingSize = c.compressionMatrix.size();
	int distinct = c.compressionMatrix[0].size();

	std::vector<int> tmpValueArray(distinct);
	for(int i=0;i<values.size();i++) {
		int realIndex = c.indexMap.at(indices[i]);
		tmpValueArray[realIndex] = values[i];
	}

	std::vector<real> compressedArray(resultingSize);
	
	for(int i=0;i<resultingSize;i++) {
		real tmp = 0.0;
		for(int j=0;j<distinct;j++) {
			tmp+=c.compressionMatrix[i][j] * tmpValueArray[j];
		}
		compressedArray[i] = tmp;
	}
	return compressedArray;
}