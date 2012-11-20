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

#ifdef NDEBUG
#define DEFAULT_CHECK_BOUNDS NO_CHECK_BOUNDS
#else
#define DEFAULT_CHECK_BOUNDS CHECK_BOUNDS
#endif

/// An n-dimensional vector.
class Vector : public Object
{
public:
    enum BoundsCheckingStatus {NO_CHECK_BOUNDS=0, CHECK_BOUNDS=1};
    real* x;
    int n;
    Vector ();
    Vector (int N_, real* y, enum BoundsCheckingStatus check = DEFAULT_CHECK_BOUNDS);
    explicit Vector (int N_, enum BoundsCheckingStatus check = DEFAULT_CHECK_BOUNDS);
    explicit Vector (uint N_, enum BoundsCheckingStatus check = DEFAULT_CHECK_BOUNDS);
    static Vector Unity(int N_, enum BoundsCheckingStatus check = DEFAULT_CHECK_BOUNDS);
    explicit Vector (real x, enum BoundsCheckingStatus check = DEFAULT_CHECK_BOUNDS);

    Vector (const Vector& rhs);
    Vector (const std::vector<real>& rhs);
    ~Vector ();
    Vector& operator= (const real& rhs);
    Vector& operator= (const Vector& rhs);
    void Clear();
    void Resize(int N_);
    void CheckBounds(enum BoundsCheckingStatus check = CHECK_BOUNDS)
    {
        checking_bounds = check;
    }
    real& operator[] (int index); ///< return element for read-write
    const real& operator[] (int index) const; ///< return element for read
    /// return element for read-write
    inline real& operator() (int index) 
    {
        return (*this)[index];
    }
    /// return element for read
    inline const real& operator() (int index) const 
    {
        return (*this)[index];
    }
    int Size() const { return n;}
    real logSum() const;
    real Sum() const;
    real Sum(int start, int end) const;
    real Norm(real p) const;
    real L1Norm() const;
    real L2Norm() const;
    real SquareNorm() const;
	const bool operator< (const real& rhs) const;
    const bool operator< (const Vector& rhs) const;
	const bool operator> (const real& rhs) const;
    const bool operator> (const Vector& rhs) const;
    const bool operator== (const Vector& rhs) const;
    const Vector operator+ (const Vector& rhs) const;
    const Vector operator- (const Vector& rhs) const;
    const Vector operator* (const Vector& rhs) const;
    const Vector operator/ (const Vector& rhs) const;
    Vector& operator+= (const Vector& rhs);
    Vector& operator-= (const Vector& rhs);
    Vector& operator*= (const Vector& rhs);
    Vector& operator/= (const Vector& rhs);
    const Vector operator+ (const real& rhs) const;
    const Vector operator- (const real& rhs) const;
	const Vector operator- () const;
    const Vector operator* (const real& rhs) const;
    const Vector operator/ (const real& rhs) const;
    Vector& operator+= (const real& rhs);
    Vector& operator-= (const real& rhs);
    Vector& operator*= (const real& rhs);
    Vector& operator/= (const real& rhs);
	
    void print(FILE* f) const;
    void printf(FILE* f) const;
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

template<typename T>
std::vector<T> convert(const Vector& x)
{
    std::vector<T> y(x.Size());
    for (int i=0; i<x.Size(); ++i) {
        y[i] = x(i);
    }
    return y;
}

Vector* NewVector (int n);///< make a new vector of length n
void CopyVector (Vector* const lhs, const Vector * const rhs); ///< Copy one vector to another.
int DeleteVector (Vector* vector); ///< Delete vector


class Matrix;

inline void SoftMax (Vector& src, Vector& dst, real beta)
{                                            
    int n = src.Size();
    if (dst.Size() != n) {
        dst.Resize(n);
    }                                           
    SoftMax(n, &src[0], &dst[0], beta);
}

inline void SoftMin (Vector& src, Vector& dst, real beta)
{                                            
    int n = src.Size();
    if (dst.Size() != n) {
        dst.Resize(n);
    }                                           
    SoftMin(n, &src[0], &dst[0], beta);
}

void Add (const Vector* lhs, const Vector* rhs, Vector* res); 
void Sub (const Vector* lhs, const Vector* rhs, Vector* res);
void Mul (const Vector* lhs, const Vector* rhs, Vector* res);
void Div (const Vector* lhs, const Vector* rhs, Vector* res);
real Product (const Vector* const lhs, const Vector* const rhs);
void Product (const Vector* lhs, const Vector* rhs, Matrix* res);
real Product (const Vector& lhs, const Vector& rhs);
void Product (const Vector& lhs, const Vector& rhs, Matrix& res);

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


inline real L1Norm (const Vector& lhs, const Vector& rhs)
{
    assert (lhs.n==rhs.n);
    return L1Norm (lhs.x, rhs.x, lhs.n);
}

/// Return \f$\|a-b\|\f$, the euclidean norm between two vectors.
inline real EuclideanNorm (const Vector& lhs, const Vector& rhs)
{
    assert (lhs.n==rhs.n);
    return EuclideanNorm (lhs.x, rhs.x, lhs.n);
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

inline real Volume(const Vector& x)
{
	real V = 1;
	for (int i=0; i<x.Size(); ++i) {
		V *= x[i];
	}
	return V;
}
/// Exponentiation
inline const Vector exp (const Vector& rhs)
{
    int n = rhs.Size();
    Vector lhs (n);
    for (int i=0; i<n; i++) {
        lhs.x[i] = exp(rhs[i]);
    }
    return lhs;
}

/// Power, by element
inline const Vector pow (const Vector& rhs, const real p)
{
    int n = rhs.Size();
    Vector lhs (n);
    for (int i=0; i<n; i++) {
        lhs.x[i] = pow((double) rhs[i], (double) p);
    }
    return lhs;
}


/// Hypertangentification
inline const Vector tanh (const Vector& rhs)
{
    int n = rhs.Size();
    Vector lhs (n);
    for (int i=0; i<n; i++) {
        lhs.x[i] = tanh(rhs[i]);
    }
    return lhs;
}


/// Logarithmication
inline const Vector log (const Vector& rhs)
{
    int n = rhs.Size();
    Vector lhs (n);
    for (int i=0; i<n; i++) {
        lhs.x[i] = log(rhs[i]);
    }
    return lhs;
}

/// Absolute value
inline const Vector abs (const Vector& rhs)
{
    int n = rhs.Size();
    Vector lhs (n);
    for (int i=0; i<n; i++) {
        lhs.x[i] = fabs(rhs[i]);
    }
    return lhs;
}

/// logAdd
inline const Vector logAdd (const Vector& x, const Vector& y)
{
	int n = x.Size();
	assert (x.Size() == y.Size());
	Vector z(n);
	for (int i=0; i<n; ++i) {
		z(i) = logAdd(x(i), y(i));
	}
	return z;
}

/// exponentiate to a target vector
void exp(const Vector& v, Vector& res);

/// logarthimicate to a target vector
void log(const Vector& v, Vector& res);


/**@}*/

#endif /* VECTOR_H */
