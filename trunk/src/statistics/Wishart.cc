/* -*- Mode: C++; -*- */
// copyright (c) 2012 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "Wishart.h"
#include "ranlib.h"
#include "SpecialFunctions.h"
#include "NormalDistribution.h"

Wishart::Wishart()
    : k(1),
      n(1),
      V(Matrix::Unity(1,1))
{
    
}

Wishart::Wishart(real n_, const Matrix& V, bool is_covariance)
    : k(V_.Rows()),
      n(n_),
{
    assert(V.Rows() == V.Columns());
    if (is_coveriance) {
        setCovariance(V);
    } else {
        setPrecision(V);
    }
}

Wishart::~Wishart()
{
    
}

void Wishart::generate(Matrix& X) const
{
    Serror("Not implemented\n");
}


///Smith & Hocking, "Wishart Variate Generator"
Matrix Wishart::generate() const
{
	//time_t ltime;
//	struct tm *today;
//	time(&ltime);
//	today = localtime(&ltime);
//	setall(today->tm_sec, 0.1*today->tm_sec);
//	
	NormalDistribution norm;
	Matrix T = Covariance.Cholesky();
	Matrix B(k,k);
	
	for(int i = 0; i < k; ++i){
	    real r = (real)genchi((real)(n - i));
		B(i,i) = sqrt(r);
	}
	
	for(int i = 0; i < k; ++i){
		for(int j = (i + 1); j < k; ++j){
			B(i,j) = norm.generate();
		}
	}
    
	Matrix X = B*T;
	return (Transpose(X) * X) / n;
}

/** the log pdf */
real Wishart::log_pdf(const Matrix& X) const
{
    assert(X.isSymmetric());
    static real log_2 = log(2.0);
    static real log_pi = log(M_PI);

    real rk = (real) k;

    real log_c = - (0.5 * rk * n * log_2
                    + 0.25 * rk * (rk - 1.0) * log_pi);

    for (int j=0; j<k; ++j) {
        log_c -= logGamma(0.5 * (n - j));
    }

    real trace_VX = 0.0;
    for (int i=0; i<k; ++i) {
        for (int j=0; j<k; ++j) {
            trace_VX += Precision(i,j) * X(i,j);
        }
    }

    real det_V = Precision.det();
    real det_X = X.det();

    real log_p = log_c 
        + 0.5 * n * log(det_V)
        + 0.5 * (n - rk - 1.0) * log(det_X)
        - 0.5 * trace_VX;

    return log_p;

}


