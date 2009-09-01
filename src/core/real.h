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

#define INF std::numeric_limits<real>::infinity()

#define LOG_ONE 0.0
#define LOG_2_PI 1.83787706640934548355
#define LOG_ZERO -INF


#ifndef INT_MAX
#define INT_MAX std::numeric_limits<int>::max()
#endif

#ifndef uint
typedef unsigned int uint;
#endif

#ifndef ulong
typedef unsigned long ulong;
#endif


#endif /* REAL_H */
