// -*- Mode: c++ -*-
#ifndef BOUNDS_H
#define BOUNDS_H

#include "real.h"

real WeissmanBound(int n_outcomes, int n_samples, real delta);
real HoeffdingBound(int n_outcomes, int n_samples, real delta);


#endif
