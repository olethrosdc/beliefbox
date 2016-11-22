#ifndef READ_DATA_H
#define READ_DATA_H

#include <vector>

#include "real.h"
#include "GPBandit.h"
#include "CompressSensing.h"


int readIntoTupleVector(std::vector<std::vector<int> > &descriptorCount, std::vector<std::vector<int> > &descriptorIndex, 
			std::vector<real> &rewardData, const char* fname);
void generateData(Matrix oldPoints_, Vector oldRewards_, Matrix &newPoints_, Vector &newRewards, int num);

#endif
