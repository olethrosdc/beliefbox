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
#include "GammaDistribution.h"
#include "NormalDistribution.h"
#include "ReadFile.h"

real Test(std::vector<real>& x, int T, int K);

int main (int argc, char** argv)
{

    int number_of_samples = 10e2;

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
        x[t] = data(t, 0);
    }
    real posterior = Test(x, T, number_of_samples);
    printf ("%f\n", posterior);
#endif


    return 0;
}

real Test(std::vector<real>& x, int T, int K)
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
#if 1
    printf ("%d %f %f %f %f\n",
            T, log_gamma_pdf, log_norm_pdf,
            log_posterior_gamma, posterior_gamma);
#endif
    NormalDistribution ML_normal;
    GammaDistribution ML_gamma;
    
    real log_ML_norm = ML_normal.setMaximumLikelihoodParameters(z);
    real log_ML_gamma = ML_gamma.setMaximumLikelihoodParameters(x, K);
    
    printf("%f %f\n", log_ML_gamma, log_ML_norm);
    return posterior_gamma;
};

#endif
