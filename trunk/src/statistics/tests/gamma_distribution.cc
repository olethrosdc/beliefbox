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

    real scaling = atof(argv[1]);
    int number_of_samples = atoi(argv[2]);
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
    ReadFloatDataASCII(data, "./weights_scaled.dat");
    int T = data.Rows();
    std::vector<real> x(T);    
    for (int t=0; t<T; ++t) {
        x[t] = data(t, 0) * scaling; // + 0.1*urandom(); // maybe add a bit of noise
    }
    //real S = Max(x);
    //for (int t=0; t<T; ++t) {
    //x[t] /= S;
    //}
    for (int k=T; k<=T; ++k) {
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
    std::vector<real> x(T); // x is the set of truncated data
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
    for (int t=0; t<T; ++t) {
        y[t] = x[t];
    }
    real log_norm_pdf = normal_prior.LogLikelihoodLogNormal(z, K);

    real a = 1;
    GammaDistributionUnknownShapeScale gamma_prior(a, a, a);
    real log_gamma_pdf =  gamma_prior.LogLikelihood(x, K);

    GammaExponentialPrior exponential_prior(a, a);
    real log_exp_pdf = exponential_prior.LogLikelihood(x, K);

    BetaDistributionMCPrior beta_prior(a, a);
    real log_beta_pdf = beta_prior.LogLikelihood(x, K);
    
    int n_c = 4;
    Vector log_P (n_c);
    log_P(0) = log_norm_pdf;
    log_P(1) = log_exp_pdf;
    log_P(2) = log_gamma_pdf;
    log_P(3) = log_beta_pdf;
    Vector log_prior (n_c);
    for (int i=0; i<n_c; ++i) {
        log_prior(i) = log(1.0/(real) n_c);
    }
    Vector log_posterior = log_prior + log_P;
    log_posterior -= log_posterior.logSum();
    
    Vector posterior = exp(log_posterior);
    log_posterior.print(stdout); printf("# Posterior (norm, exp, gamma, beta)\n");

    NormalDistribution ML_normal;
    GammaDistribution ML_gamma;
    ExponentialDistribution ML_exp;
    BetaDistribution ML_beta;
    
    real log_ML_norm = ML_normal.setMaximumLikelihoodParameters(z);
    real log_ML_gamma = ML_gamma.setMaximumLikelihoodParameters(x, K);
    real log_ML_beta = ML_beta.setMaximumLikelihoodParameters(y, K);
    real log_ML_exp = ML_exp.setMaximumLikelihoodParameters(x);
    printf("norm_ml_params = [%f %f];\n", ML_normal.m, ML_normal.s);
    printf("gamma_ml_params = [%f %f];\n", ML_gamma.alpha, ML_gamma.beta);
    printf("beta_ml_params=[%f %f];\n", ML_beta.alpha, ML_beta.beta);
    printf("exp_ml_params=[%f];\n", ML_exp.l);
    printf("%f %f %f %f # ML log-norm, exp, gamma, beta\n", log_ML_norm, log_ML_exp, log_ML_gamma, log_ML_beta);
    return log_posterior;
};

#endif
