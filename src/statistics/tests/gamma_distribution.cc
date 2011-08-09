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

#ifdef MAKE_MAIN
#include "Random.h"
#include <vector>
#include "BetaDistribution.h"
#include "GammaDistribution.h"
#include "NormalDistribution.h"
#include "ExponentialDistribution.h"
#include "ReadFile.h"

Vector Test(std::vector<real>& x, int T, int K);

int main (int argc, char** argv)
{

    int number_of_samples = 1000;

#if 0
    int T = 1000;

	std::vector<real> x(T);    
	real p = 0.0;
	//GammaDistribution gamma_source(1.0, 1.0);
	NormalDistribution normal_source(-3.0, 1.0);
	BetaDistribution gamma_source(1.0, 1.0);
	for (int t=0; t<T; ++t) {
		if (urandom() < p) {
			x[t] = gamma_source.generate();
		} else {
			x[t] = exp(normal_source.generate());
		}
		//printf ("%f #data\n", x[t]);
	}
	Vector posterior = Test(x, T, number_of_samples);
    printf("# posterior (lognorm, exp, gamma)\n");
	exp(posterior).print(stdout);
    printf("# log posterior (lognorm, exp, gamma)\n");
    posterior.print(stdout);

#else
    Matrix data;
    ReadFloatDataASCII(data, "./weights_raw.dat");
    int T = data.Rows();
    std::vector<real> x(T);    
    for (int t=0; t<T; ++t) {
        x[t] = (0.001 +data(t, 0))*0.999; // + 0.1*urandom(); // maybe add a bit of noise
    }
    real S = Max(x);
    for (int t=0; t<T; ++t) {
        x[t] /= S;
    }
    for (int k=1; k<T; ++k) {
        int t = k * k;
        if (t > T) {
            t = T;
        }
        Vector posterior = Test(x, t, number_of_samples);
        if (t >= T) {
            break;
        }
    }

#endif


    return 0;
}

Vector Test(std::vector<real>& x_orig, int T, int K)
{
    std::vector<real> x(T);
    for (int t=0; t<T; ++t) {
        x[t] = x_orig[t];
    }
    //printf ("Running test for T = %d\n", T);
    //printf ("=======================\n\n");
    NormalUnknownMeanPrecision normal_prior;
    //real log_norm_pdf = 0;
    std::vector<real> z(T);
    for (int t=0; t<T; ++t) {
        z[t] = log(x[t]);
    }
    std::vector<real> y(T);
    real min_x = Min(x);
    real scale_x = 1.0 / (Max(x) - min_x);
    for (int t=0; t<T; ++t) {
        y[t] = scale_x * (x[t] - min_x);
        //printf ("%f # Y\n", y[t]);
    }
    real log_norm_pdf = normal_prior.LogLikelihood(z, K);

    real a = 1;
    GammaDistributionUnknownShapeScale gamma_prior(a, a, a);
    real log_gamma_pdf =  gamma_prior.LogLikelihood(x, K);

    GammaExponentialPrior exponential_prior(a, a);
    real log_exp_pdf = exponential_prior.LogLikelihood(x, K);
    
    Vector log_P (3);
    log_P(0) = log_norm_pdf;
    log_P(1) = log_exp_pdf;
    log_P(2) = log_gamma_pdf;
    Vector log_prior (3);
    for (int i=0; i<3; ++i) {
        log_prior(i) = log(1.0/3.0);
    }
    Vector log_posterior = log_prior + log_P;
    log_posterior -= log_posterior.logSum();
    
    Vector posterior = exp(log_posterior);
    //log_posterior.print(stdout);

    NormalDistribution ML_normal;
    GammaDistribution ML_gamma;
    ExponentialDistribution ML_exp;
    BetaDistribution ML_beta;
    
    real log_ML_norm = ML_normal.setMaximumLikelihoodParameters(z);
    real log_ML_gamma = ML_gamma.setMaximumLikelihoodParameters(x, K);
    real log_ML_beta = ML_beta.setMaximumLikelihoodParameters(y, K);
    real log_ML_exp = ML_exp.setMaximumLikelihoodParameters(x);
    printf("ML Beta: %f %f\n", ML_beta.alpha, ML_beta.beta);
    printf("%f %f %f %f # ML log-norm, exp, gamma, beta\n", log_ML_norm, log_ML_exp, log_ML_gamma, log_ML_beta);
    return log_posterior;
};

#endif
