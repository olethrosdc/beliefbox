/* -*- Mode: c++ -*- */
/* VER: $Id: MathFunctions.h,v 1.2 2006/11/06 15:48:53 cdimitrakakis Exp cdimitrakakis $ */
// copyright (c) 2004 by Christos Dimitrakakis <dimitrak@idiap.ch>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "SpecialFunctions.h"

real Beta(real x, real y)
{
    return exp(logGamma(x) + logGamma(y) - logGamma(x+y));
}

real logBeta(real x, real y)
{
    return logGamma(x) + logGamma(y) - logGamma(x+y);
}


real Beta(Vector& x)
{
    return exp(logBeta(x));
}

real logBeta(Vector& x)
{
	int N = x.Size();
	real sum = 0.0;
	real logsum = 0.0;
	for (int i=0; i<N; i++) {
		sum += x[i];
		logsum += logGamma(x[i]);
	}
    return exp(logsum - sum);
}



float betacf(float a, float b, float x)
//Used by betai: Evaluates continued fraction for incomplete beta function by modified Lentz's method (§5.2).
{
    const int MAXIT = 100;
    const float EPS = 3.0e-7;
    const float FPMIN = 1.0e-30;

  int m,m2;
  float aa,c,d,del,h,qab,qam,qap;
  // These q's will be used in factors that occur
  // in the coefficients (6.4.6).
  qab=a+b;
                                           
  qap=a+1.0;
  qam=a-1.0;
  // First step of Lentz's method.
  c=1.0;
  d=1.0-qab*x/qap;
  if (fabs(d) < FPMIN) d=FPMIN;
  d=1.0/d;
  h=d;
  for (m=1;m<=MAXIT;m++) {
      m2=2*m;
      aa=m*(b-m)*x/((qam+m2)*(a+m2));
      // One step (the even one) of the recurrence.
      d=1.0+aa*d;
      if (fabs(d) < FPMIN) d=FPMIN;
      c=1.0+aa/c;
      if (fabs(c) < FPMIN) c=FPMIN;
      d=1.0/d;
      h *= d*c;
      aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
      // Next step of the recurrence (the odd one).
      d=1.0+aa*d;
      if (fabs(d) < FPMIN) d=FPMIN;
      c=1.0+aa/c;
      if (fabs(c) < FPMIN) c=FPMIN;
      d=1.0/d;
      del=d*c;
      h *= del;
      // Are we done?
      if (fabs(del-1.0) < EPS) break;
  }
  if (m > MAXIT) fprintf(stderr, "a or b too big, or MAXIT too small in betacf");
  return h;
}

//Returns the incomplete beta function Ix (a, b).
float BetaInc(float x, float a, float b)

{
    //float betacf(float a, float b, float x);
    //float gammln(float xx);

    float bt;
    if (x < 0.0 || x > 1.0) fprintf(stderr, "Bad x in routine betai");
    if (x == 0.0 || x == 1.0) {
        bt=0.0;
    }  else {
        // Factors in front of the continued fraction.
        bt=exp(logGamma(a+b) - logGamma(a) - logGamma(b) + a*log(x)+b*log(1.0-x));
    }

    if (x < (a+1.0)/(a+b+2.0)) {
        // Use continued fraction directly.
        return bt*betacf(a,b,x)/a;
    } else {
        // Use continued fraction after making the symmetry transformation
        return 1.0-bt*betacf(b,a,1.0-x)/b;
    }
}

unsigned long binomial (int n, int k)
{
    if (n <= 0 || k < 0 || k > n)
        return 0;

    if (k > n/2)
        k = n-k; // faster

    //long long n_minus_k = n-k;
    long double accum = 1;
    int k2 = k;
    for (int i = 1; i++ < k2;) {
        //accum = accum * (n_minus_k+i) / i;
        accum *= n / k;
        n--;
        k--;
    }

    return (unsigned long) round(accum + 0.5); // avoid rounding error
}



