/* -*- Mode: C++; -*- */
/* "VER: $Id: Vector.h,v 1.1 2006/11/07 18:03:35 cdimitrakakis Exp cdimitrakakis $" */

#ifndef VECTOR_H
#define VECTOR_H

#include "debug.h"
#include "real.h"
#include "MathFunctions.h"
#include "Object.h"
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cmath>

/**
   \ingroup MathGroup
*/
/*@{*/

/** 
    \file Vector.h
	
    \brief Vector and matrix computations.
*/

/// An n-dimensional vector.
class Vector : public Object
{
public:
    enum BoundsCheckingStatus {NO_CHECK_BOUNDS=0, CHECK_BOUNDS=1};
    real* x;
    int n;
    Vector ();
    Vector (int N_, real* y, enum BoundsCheckingStatus check = NO_CHECK_BOUNDS);
    Vector (int N_, enum BoundsCheckingStatus check = NO_CHECK_BOUNDS);
    Vector (const Vector& rhs);
    ~Vector ();
    Vector& operator= (const Vector& rhs);
    void Resize(int N_);
    void CheckBounds(enum BoundsCheckingStatus check = CHECK_BOUNDS)
    {
        checking_bounds = check;
    }
    real& operator[] (int index); ///< return element for read-write
    const real& operator[] (int index) const; ///< return element for read
    int Size() const { return n;}
    real Sum() const;
    real Sum(int start, int end) const;
    Vector operator+ (const Vector& rhs);
    Vector operator- (const Vector& rhs);
    Vector operator* (const Vector& rhs);
    Vector operator/ (const Vector& rhs);
    Vector& operator+= (const Vector& rhs);
    Vector& operator-= (const Vector& rhs);
    Vector& operator*= (const Vector& rhs);
    Vector& operator/= (const Vector& rhs);
    Vector operator+ (const real& rhs);
    Vector operator- (const real& rhs);
    Vector operator* (const real& rhs);
    Vector operator/ (const real& rhs);
    Vector& operator+= (const real& rhs);
    Vector& operator-= (const real& rhs);
    Vector& operator*= (const real& rhs);
    Vector& operator/= (const real& rhs);
private:
    int maxN;
    enum BoundsCheckingStatus checking_bounds;
};


/// Use this to define an n_1 x n_2 x ... x n_N lattice
typedef struct Lattice_ {
    int dim; ///< number of dimensions
    int* n; ///< number of dimensions per axis
    real* x; ///< data
} Lattice;


Vector* NewVector (int n);///< make a new vector of length n
void CopyVector (Vector* const lhs, const Vector * const rhs); ///< Copy one vector to another.
int DeleteVector (Vector* vector); ///< Delete vector

/// Get i-th component of vector v
inline real GetVal (Vector* v, int i) {
    assert ((i>=0)&&(i<v->n));
    return v->x[i];
}

/// Set i-th component of vector v
inline void PutVal (Vector* v, int i, real x) {
    assert ((i>=0)&&(i<v->n));
    v->x[i] = x;
}


class Matrix;

void Add (const Vector* lhs, const Vector* rhs, Vector* res); 
void Sub (const Vector* lhs, const Vector* rhs, Vector* res);
void Mul (const Vector* lhs, const Vector* rhs, Vector* res);
void Div (const Vector* lhs, const Vector* rhs, Vector* res);
real Product (const Vector* const lhs, const Vector* const rhs);
void Product (const Vector* lhs, const Vector* rhs, Matrix* res);
//real EuclideanNorm (const Vector* lhs, const Vector* rhs);
//real SquareNorm (const Vector* lhs, const Vector* rhs);
/// Return \f$\|a-b\|_1\f$, the L1 norm between two vectors.
inline real L1Norm (const Vector* lhs, const Vector* rhs)
{
    assert (lhs->n==rhs->n);
    return L1Norm (lhs->x, rhs->x, lhs->n);
}

/// Return \f$\|a-b\|\f$, the euclidean norm between two vectors.
inline real EuclideanNorm (const Vector* lhs, const Vector* rhs)
{
    assert (lhs->n==rhs->n);
    return EuclideanNorm (lhs->x, rhs->x, lhs->n);
}


/// Return \f$\|a-b\|^2\f$, the square of the euclidean norm between two
/// vectors.
inline real SquareNorm (const Vector* lhs, const Vector* rhs)
{
    assert (lhs->n==rhs->n);
    return SquareNorm (lhs->x, rhs->x, lhs->n);
}

/// Get maximum element
inline real Max(const Vector* v)
{
    return Max(v->Size(), v->x);
}

/// Get minimum element
inline real Min(const Vector* v)
{
    return Min(v->Size(), v->x);
}

/// Get maximum element
inline int ArgMax(const Vector* v)
{
    return ArgMax(v->Size(), v->x);
}

/// Get minimum element
inline int ArgMin(const Vector* v)
{
    return ArgMin(v->Size(), v->x);
}

inline real Span(const Vector* v)
{
    return Max(v) - Min(v);
}

/// Get maximum element
inline real Max(const Vector& v)
{
    return Max(v.Size(), v.x);
}

/// Get minimum element
inline real Min(const Vector& v)
{
    return Min(v.Size(), v.x);
}

/// Get maximum element
inline int ArgMax(const Vector& v)
{
    return ArgMax(v.Size(), v.x);
}

/// Get minimum element
inline int ArgMin(const Vector& v)
{
    return ArgMin(v.Size(), v.x);
}

inline real Span(const Vector& v)
{
    return Max(v) - Min(v);
}

Vector exp (const Vector& rhs)
{
    int n = rhs.Size();
    Vector lhs (n);
    for (int i=0; i<n; i++) {
        lhs.x[i] = exp(rhs[i]);
    }
    return lhs;
}

/// Logarithmication
Vector log (const Vector& rhs)
{
    int n = rhs.Size();
    Vector lhs (n);
    for (int i=0; i<n; i++) {
        lhs.x[i] = log(rhs[i]);
    }
    return lhs;
}
inline void exp(const Vector& v, Vector& res)
{
    assert(v.Size() == res.Size());
    for (int i=0; i<v.Size(); ++i) {
        res[i] = exp(v[i]);
    }
}

inline void log(const Vector& v, Vector& res)
{
    assert(v.Size() == res.Size());
    for (int i=0; i<v.Size(); ++i) {
        res[i] = log(v[i]);
    }
}
#if 0
real L1Norm (const Vector* lhs, const Vector* rhs);
//real L1Norm (const Vector* lhs);
real Max(const Vector* v);
real Min(const Vector* v);
int ArgMax(const Vector* v);
int ArgMin(const Vector* v);
real Span(const Vector* v);
#endif

/*@}*/

#endif /* VECTOR_H */
