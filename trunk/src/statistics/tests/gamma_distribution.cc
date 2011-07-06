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
int main (int argc, char** argv)
{
#if 0
    int T = 1000;

    std::vector<real> x(T);    
    if (1) {
        GammaDistribution source(2.0, 2.0);
        printf("Data is Gamma!\n");
        for (int t=0; t<T; ++t) {
            x[t] = source.generate();
        }
    } else {
        NormalDistribution source(2.0, 1.0);
        printf("Data is Normal!\n");
        for (int t=0; t<T; ++t) {
            x[t] = exp(source.generate());
        }
    }
    int max_T = T;
#else
    Matrix data;
    ReadFloatDataASCII(data, "./weights_scaled.dat");
    int T = data.Rows();
    std::vector<real> x(T);    
    for (int t=0; t<T; ++t) {
        x[t] = data(t, 1);
    }
    int max_T = data.Rows();
#endif

    for (T = 1; T<=max_T + 1; T+=1) {
        if (T > max_T) {
            T = max_T;
        }
    MultivariateNormalUnknownMeanPrecision normal_prior;
    real log_norm_pdf = 0;
    for (int t=0; t<T; ++t) {
        Vector z(1);
        z[0] = log(x[t]) + 3.0;
        real log_p_z = normal_prior.log_pdf(z);
        log_norm_pdf += log_p_z;
        normal_prior.Observe(z);
    }
    //printf("# %f # normal likelihood\n", log_norm_pdf);
    real a = 1;
    int K = 10e2;
    GammaDistributionUnknownShapeScale gamma_prior(a, a, a);
    real log_gamma_pdf =  gamma_prior.LogLikelihood(x, K);
    //printf ("%f %d %f\n", a, K, log_gamma_pdf);

    real prior_gamma = 0.5;
    real log_prior_gamma = log(prior_gamma);
    real log_prior_norm = log(1.0 - prior_gamma);

    real log_posterior_gamma = log_prior_gamma + log_gamma_pdf
        - logAdd(log_prior_gamma + log_gamma_pdf, log_prior_norm + log_norm_pdf);
    printf ("%d %f %f %f %f\n", T, log_gamma_pdf, log_norm_pdf, log_posterior_gamma, exp(log_posterior_gamma));
        if (T > max_T) {
            break;
        }
    }
};

#endif
