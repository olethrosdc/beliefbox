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

    int number_of_samples = 10;

#if 0
    int T = 100;
    int N = 101;
    int max_iter = 100;
    Vector mix_probabilities(N);
    Vector P_p(N);
    for (int n=0; n<N; ++n) {
        mix_probabilities[n] = (real) n / (real) (N - 1);
        P_p[n] = 0.0;
    }

    for (int iter=0; iter < max_iter; ++iter) {
        std::vector<real> x(T);    
        for (real i = 0; i<N; ++i) {
            real p = mix_probabilities[i];
            GammaDistribution gamma_source(1.0, 1.0);
            NormalDistribution normal_source(0.0, 1.0);
            for (int t=0; t<T; ++t) {
                if (urandom() < p) {
                    x[t] = gamma_source.generate();
                } else {
                    x[t] = exp(normal_source.generate());
                }
            }
            real posterior = Test(x, T, number_of_samples);
            P_p[i] += posterior;
        }
    }
    P_p /= max_iter;
    mix_probabilities.print(stdout);
    P_p.print(stdout);
#else
    Matrix data;
    ReadFloatDataASCII(data, "./weights_scaled.dat");
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
        printf ("%d ", t);
        Vector posterior = Test(x, t, number_of_samples);
        if (t >= T) {
            break;
        }
    }

#endif


    return 0;
}

Vector Test(std::vector<real>& x, int T, int K)
{

    NormalUnknownMeanPrecision normal_prior;
    //real log_norm_pdf = 0;
    std::vector<real> z(T);
    for (int t=0; t<T; ++t) {
        z[t] = log(x[t]);
    }
    real log_norm_pdf = normal_prior.LogLikelihood(z, K);
    //printf("# %f # normal likelihood\n", log_norm_pdf);
    
    real a = 1;
    GammaDistributionUnknownShapeScale gamma_prior(a, a, a);
    real log_gamma_pdf =  gamma_prior.LogLikelihood(x, K);
    //printf ("%f %d %f\n", a, K, log_gamma_pdf);

    real prior_gamma = 0.5;
    real log_prior_gamma = log(prior_gamma);
    real log_prior_norm = log(1.0 - prior_gamma);

    real log_posterior_gamma = log_prior_gamma + log_gamma_pdf
        - logAdd(log_prior_gamma + log_gamma_pdf, log_prior_norm + log_norm_pdf);
    real posterior_gamma = exp(log_posterior_gamma);

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
    log_posterior.print(stdout);

    NormalDistribution ML_normal;
    GammaDistribution ML_gamma;
    ExponentialDistribution ML_exp;
    BetaDistribution ML_beta;
    
    real log_ML_norm = ML_normal.setMaximumLikelihoodParameters(z);
    real log_ML_gamma = ML_gamma.setMaximumLikelihoodParameters(x, K);
    real log_ML_beta = ML_beta.setMaximumLikelihoodParameters(x, K);
    real log_ML_exp = ML_exp.setMaximumLikelihoodParameters(x);
    
    printf("%f %f %f %f # ML\n", log_ML_norm, log_ML_exp, log_ML_gamma, log_ML_beta);
    return log_posterior;
};

#endif
