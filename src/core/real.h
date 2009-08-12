// -*- Mode: c++ -*-

#ifndef REAL_H
#define REAL_H

#include <limits>

#ifdef real

#error "Real already defined in another header! Binaries might not link properly!"

#endif /* real */


#ifdef USE_FIXED_POINT
typedef fixed real;

#else /* USE_FIXED_POINT */

#ifdef USE_DOUBLE
typedef double real;
#else
typedef float real;
#endif

#endif /* USE_FIXED_POINT */

#define LOG_ONE 0.0
#define LOG_2_PI 1.83787706640934548355
#define LOG_ZERO -INF

#define INF std::numeric_limits<real>::infinity()


#endif /* REAL_H */
