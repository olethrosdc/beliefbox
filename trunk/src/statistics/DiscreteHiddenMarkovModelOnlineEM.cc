/* -*- Mode: c++;  -*- */
// copyright (c) 2009 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "DiscreteHiddenMarkovModelOnlineEM.h"
#include "RandomNumberGenerator.h"
#include <cassert>
#include <cmath>

#if 0
#define dbg_printf (void)
#else
#define dbg_printf printf
#endif

DiscreteHiddenMarkovModelOnlineEM::DiscreteHiddenMarkovModelOnlineEM(int n_states_, int n_observations_)  : DiscreteHiddenMarkovModel(n_states_, n_observations_), q(n_states)
{
    InitPhi();
    T = 0;
    Reset();
}


DiscreteHiddenMarkovModelOnlineEM::DiscreteHiddenMarkovModelOnlineEM(Matrix& Pr_S, Matrix& Pr_X) : DiscreteHiddenMarkovModel(Pr_S, Pr_X), q(n_states)
{
    InitPhi();
    T = 0;
    Reset();
}

void DiscreteHiddenMarkovModelOnlineEM::InitPhi()
{
    std::vector<int> x(4);
    x[0] = n_states;
    x[1] = n_states;
    x[2] = n_states;
    x[3] = n_observations;
    Phi = new Tensor(x);
    q[0] = 1;
}

void DiscreteHiddenMarkovModelOnlineEM::Reset()
{
    DiscreteHiddenMarkovModel::Reset();
    Phi->Reset();
    real p = 1.0 / (real) n_states;
    for (int i=0; i<n_states; ++i) {
        q[i] = p;
    }
     for (int i=0; i<n_states; ++i) {
         for (int j=0; j<n_states; ++j) {
             PrS(i,j) = 1.0 / (real) n_states;
         }
         PrS(i,i)+=1;
         for (int j=0; j<n_observations; ++j) {
             PrX(i,j) = 1.0 / n_observations;
         }
         PrX(i,i)+=1;
         
         real sum = 0;
         for (int j=0; j<n_states; ++j) {
             sum += PrS(i,j);
         }
         real invsum = 1.0 / sum;
         for (int j=0; j<n_states; ++j) {
             PrS(i,j) *= invsum;
         }

         sum = 0;
         for (int j=0; j<n_observations; ++j) {
             sum += PrX(i,j);
         }         
         invsum = 1.0 / sum;
         for (int j=0; j<n_observations; ++j) {
             PrX(i,j) *= invsum;
         }
    }

    T = 0;
}

DiscreteHiddenMarkovModelOnlineEM::~DiscreteHiddenMarkovModelOnlineEM()
{
    // nothing to do
    delete Phi;
}

/** Update EM estimate
   
Based on: G. Mongillo and S. Deneve, "Online learning with hidden
Markov models", Neural computation 20, 1706-1716 (2008).
*/
real DiscreteHiddenMarkovModelOnlineEM::Observe(int x)
{
    assert(x >= 0 && x < n_observations);
    T++;
    real eta = 1.0 / (T);
    real fudge = 0.00000001;
    // Calculate gamma
    Matrix gamma(n_states, n_states);
    Matrix AB(n_states, n_states);
    for (int l=0; l<n_states; ++l) {
        for (int h=0; h<n_states; ++h) {
            AB(l,h) = PrS(l, h) * PrX(h, x);
            //dbg_printf ("P_S(%d %d) = %f, P_X(%d %d)= %f\n", l, h, PrS(l,h), h, x, PrX(h, x));
		}
    }

    real Z_gamma = 0;
    for (int m=0; m<n_states; ++m) {
        for (int n=0; n<n_states; ++n) {
            Z_gamma += q[m] * AB(m,n);
        }
    }

    for (int l=0; l<n_states; ++l) {
        for (int h=0; h<n_states; ++h) {
            gamma(l, h) = AB(l,h) / Z_gamma;
            dbg_printf("gamma(%d %d) = %f\n", l, h, gamma(l, h));
        }
    }

    // update phi
    {
        std::vector<int> index(4);
        int& h = index[0];
        int& i = index[1];
        int& j = index[2];
        int& k = index[3];
        for (i=0; i<n_states; ++i) {
            for (j=0; j<n_states; ++j) {
                for (k=0; k<n_observations; ++k) {
                    // for each ijk combination, no other is used.
                    // copy previous values
                    Vector phi_prev(n_states);
                    for (h=0; h<n_states; ++h) {
                        phi_prev(h) = Phi->Y(index);
                        dbg_printf ("phi_{%d %d %d}^%d(T-1) = %f\n",i, j, k, h, phi_prev(h));
                    }
                    
                    // caclulate new set of values for phi_{ijk}^h
                    for (h=0; h<n_states; ++h) {
                        real sum = 0;
                        for (int l=0; l<n_states; ++l)  {
                            real hq = 0;
                            if (x == k && i == l && j == h) {
                                hq = q(l);
                            }
                            sum += gamma(l, h) * (phi_prev(l) + eta  * (hq - phi_prev(l)));
                        } // l
                        Phi->Y(index) = sum;
                        dbg_printf ("phi_{%d %d %d}^%d(T) = %f\n", i, j, k, h, Phi->Y(index));
                    } // k
                } // h
            } // j
        } // i
    }

    // update q
    Vector q_old = q;
    real log_p_x = 0;
    for (int l=0; l<n_states; ++l) {
        real sum = 0.0;
        for (int m=0; m<n_states; ++m) {
            sum += gamma(m, l) * q_old[m];
        }
        q(l) = sum;
        log_p_x += log(q(l)) + log(PrX(l, x));
        dbg_printf ("q(%d) = %f\n", l, q(l));
    }
    
    // update A
    {
        std::vector<int> index(4);
        int& h = index[0];
        int& i = index[1];
        int& j = index[2];
        int& k = index[3];
        for (i=0; i<n_states; ++i) {
            real sum_ai = 0.0;
            for (j=0; j<n_states; ++j) {
                real sum = 0.0;
                for (k=0; k<n_observations; ++k) {
                    for (h=0; h<n_states; ++h) {
                        sum += Phi->Y(index) + fudge;
                    }
                }
                PrS(i, j) = sum;
                sum_ai += sum;
            }
            dbg_printf ("sum_%d: %f\n", i, sum_ai);
            real inv_sum = 1.0 / sum_ai;
            for (j=0; j<n_states; ++j) {
                PrS(i, j) *= inv_sum;
                dbg_printf ("%f ", PrS(i, j));
            }
            dbg_printf("# P_S\n");
        }
    }
    // update B
    {
        std::vector<int> index(4);
        int& h = index[0];
        int& i = index[1];
        int& j = index[2];
        int& k = index[3];
        for (j=0; j<n_states; ++j) {
            real sum_bj = 0.0;
            for (k=0; k<n_observations; ++k) {
                real sum = 0.0;
                for (i=0; i<n_states; ++i) {
                    for (h=0; h<n_states; ++h) {
                        sum += Phi->Y(index) + fudge;
                    } // h
                } // i
                PrX(j, k) = sum;
                sum_bj += sum;
            } // k
            real inv_sum = 1.0 / sum_bj;
            for (k=0; k<n_observations; ++k) {
                PrX(j, k) *= inv_sum;
                dbg_printf ("%f ", PrX(j, k));
            } // k
            dbg_printf("# P_X\n");
        } // j
    }
    return log_p_x;
}

Vector DiscreteHiddenMarkovModelOnlineEM::getPrediction()
{

    Vector p(n_states);
    for (int i=0; i<n_states; ++i) {
        real sum = 0;
        for (int j=0; j<n_states; ++j) {
            sum += q(i) * PrS(i, j);
        }
        p(i) = sum;
    }

    Vector v(n_observations);
    for (int i=0; i<n_states; ++i) {
        real sum = 0;
        for (int k=0; k<n_observations; ++k) {
            sum += p(i) * PrX(i, k);
        }
        v(i) = sum;
    }

    return v;
}
