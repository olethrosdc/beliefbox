/* -*- Mode: C++; -*- */
// copyright (c) 2011 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include <vector>
class Matrix;

/// A simple K-Fold class for easy access
class KFold
{
protected:
	Matrix& data; ///< the data, to which we need access
	int K; ///< number of folds
	int T; ///< total amount of data
	int N; ///< number of columns in data
	std::vector<int> assignment; ///< assignment of data to folds
public:
	KFold(Matrix& data_, int K_);
	Matrix getTrainFold(int n);
	Matrix getTestFold(int n);
};
