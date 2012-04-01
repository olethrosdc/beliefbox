// -*- Mode: c++ -*-
#include "Bounds.h"
#include <cmath>

real WeissmanBound(int n_outcomes, int n_samples, real delta)
{
	if (n_samples > 0) {
		return sqrt(2.0 * ((real) (n_outcomes - 1) * log(2.0) - log(delta)) / (real) n_samples);
	}
	return 1.0;
}

real HoeffdingBound(int n_outcomes, int n_samples, real delta)
{
	if (n_samples > 0) {
		return sqrt(log((real) n_outcomes / delta) / (2.0 * (real) n_samples));
	}
	return 1.0;
}
