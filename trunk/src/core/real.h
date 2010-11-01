// -*- Mode: c++ -*-
/** \file Defines the real type. 
    
    It defaults to single precision.

    There is a possibility to use double precision via USE_DOUBLE.

    The USE_FIXED flag, that uses fixed point, is unimplemented.
    
    Possibly we could look at http://gmplib.org/ for a variable
    precision implementation.

 */
#ifndef REAL_H
#define REAL_H

#include <limits>
#include <cmath>

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
#define MIN_PRECISION std::numeric_limits<real>::min()
#define REAL_RANGE std::numeric_limits<real>::max()

#define LOG_ONE 0.0
#define LOG_2_PI 1.83787706640934548355
#define LOG_ZERO -INF



#ifndef uint
typedef unsigned int uint;
#endif

#ifndef ulong
typedef unsigned long ulong;
#endif


#endif /* REAL_H */
