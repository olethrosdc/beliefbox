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
#include "ContextTreeKDTree.h"

Vector Test(std::vector<real>& x, int T, int K);

int main (int argc, char** argv)
{

    if (argc != 4) {
        fprintf(stderr, "Usage: gamma_distribution scaling N_samples filename\n");
        exit(-1);
    }
    real scaling = atof(argv[1]);
    int number_of_samples = atoi(argv[2]);
    char* fname = argv[3];
        

    Matrix data;
    if (!ReadFloatDataASCII(data, fname)) {
        Serror("No data read from file %s\n", fname);
        exit(-1);
    }
    int T = data.Rows();
    std::vector<real> x(T);    
    for (int t=0; t<T; ++t) {
        x[t] = data(t, 0) * scaling; // + 0.1*urandom(); // maybe add a bit of noise
    }

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

    return 0;
}

Vector Test(std::vector<real>& x_orig, int T, int K)
{

    real alpha = 1.0;

    std::vector<real> x(T); // x is the set of truncated data
    for (int t=0; t<T; ++t) {
        x[t] = x_orig[t];
    }
    //printf ("Running test for T = %d\n", T);
    //printf ("=======================\n\n");
    NormalUnknownMeanPrecision normal_prior(alpha, alpha);
    //real log_norm_pdf = 0;
    std::vector<real> z(T); // x is the set of logarithmic data
    for (int t=0; t<T; ++t) {
        z[t] = log(x[t]);
    }
    std::vector<real> y(T); // y is the same as x (!)
    for (int t=0; t<T; ++t) {
        y[t] = x[t];
    }
    real min_x = Min(x);
    real max_x = Max(x);
    real log_norm_pdf = normal_prior.LogLikelihoodLogNormal(z, K);

    GammaDistributionUnknownShapeScale gamma_prior(alpha, alpha, alpha);
    real log_gamma_pdf =  gamma_prior.LogLikelihood(x, K);

    GammaExponentialPrior exponential_prior(alpha, alpha);
    real log_exp_pdf = exponential_prior.LogLikelihood(x, K);

    BetaDistributionMCPrior beta_prior(alpha, alpha);
    real log_beta_pdf = beta_prior.LogLikelihood(x, K);
    
    Vector lower_bound(1);
    Vector upper_bound(1);
    lower_bound(0) = min_x - 0.5;
    upper_bound(0) = max_x + 0.5;
    ContextTreeKDTree kd_tree (2, 128, lower_bound, upper_bound);
    real log_kd_tree = 0;
    for (int t=0; t<T; ++t) {
        Vector v(1);
        v(0) = x[t];
        real log_p = log(kd_tree.Observe(v));
        if (log_p < -40) {
            log_p = -40;
        }
        log_kd_tree += log_p;
        //printf ("%f %f %f\n", x[t], log_p);
    }
    #if 0
    for (real X = lower_bound(0); X<=upper_bound(0); X+= 0.001) {
        Vector v(1);
        v(0) = X;
        real p_tree = kd_tree.pdf(v);
        
        printf ("%f %f # PDF\n", X, p_tree);
    }
#endif
    int n_c = 5;
    Vector log_P (n_c);
    log_P(0) = log_norm_pdf;
    log_P(1) = log_exp_pdf;
    log_P(2) = log_gamma_pdf;
    log_P(3) = log_beta_pdf;
    log_P(4) = log_kd_tree;
    Vector log_prior (n_c);
    for (int i=0; i<n_c; ++i) {
        log_prior(i) = log(1.0/(real) n_c);
    }
    Vector log_posterior = log_prior + log_P;
    log_posterior -= log_posterior.logSum();
    
    Vector posterior = exp(log_posterior);
    log_posterior.print(stdout); printf("# Posterior (norm, exp, gamma, beta, tree)\n");

    return log_posterior;
    NormalDistribution ML_normal;
    GammaDistribution ML_gamma;
    ExponentialDistribution ML_exp;
    BetaDistribution ML_beta;
    
    real log_ML_norm = ML_normal.setMaximumLikelihoodParametersLogNormal(z);
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
