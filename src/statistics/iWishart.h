/* -*- Mode: C++; -*- */
// copyright (c) 2012 by Nikolaos Tziortziotis <ntziorzi@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef iWishart_H
#define iWishart_H

#include "Distribution.h"
#include "Matrix.h"


/// Inverse Wishart probability distribution
class iWishart : public AbstractDistribution<Matrix>
{
protected:
    Matrix Precision; ///< precision matrix
    Matrix Covariance; ///< covariance matrix
public:
	int k;		///< dimensionality
	real n;		///< degrees of freedom
	iWishart(); 
	iWishart(real n_, const Matrix& V, bool is_covariance = false);
	virtual ~iWishart();
	virtual void generate(Matrix& X) const;
	virtual Matrix generate() const;
	virtual real pdf(const Matrix& X) const
	{
		return exp(log_pdf(X));
	}
	virtual real log_pdf(const Matrix& X) const;
	void setCovariance(const Matrix& V)
    {
        Covariance = V;
        Precision = V.Inverse_LU();
    }
    void setPrecision(const Matrix& V)
    {
        Covariance = V.Inverse_LU();
        Precision = V;
    }
    void Show()
    {
        logmsg("Precision:");
        Precision.print(stdout);
		
        logmsg("Covariance:");
        Covariance.print(stdout);
    }
	
};

#endif
