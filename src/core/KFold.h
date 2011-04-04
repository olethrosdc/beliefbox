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
	int n_folds; ///< number of folds
	int n_records; ///< total amount of data
	int n_columns; ///< number of columns in data
	std::vector<int> assignment; ///< assignment of data to folds
	std::vector<int> totals; ///< the total number of data per assignment
public:
	KFold(Matrix& data_, int K_);
	Matrix getTrainFold(int n, int T=0);
	Matrix getTestFold(int n, int T=0);
};
